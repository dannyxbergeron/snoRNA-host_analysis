import os
from collections import Counter
from itertools import permutations
import pandas as pd
import numpy as np


from pybedtools import BedTool as bt


in_file = snakemake.input.single_id
gene_bed_biotype = snakemake.input.gene_bed_biotype

out_file = snakemake.output.merged_windows


def load_df(file):

    df = pd.read_csv(file)
    df['merge_ids'] = [
        '_'.join(sorted([id1, id2]))
        for (id1, id2) in df[['single_id1', 'single_id2']].values
    ]
    return df


def get_snoRNA_list(file):

    df = pd.read_csv(file, sep='\t')
    sno_df = df.loc[df.gene_biotype == 'snoRNA']
    snolist = list(sno_df.gene_id)

    return snolist


def snoRNA_side(df_int_, snolist):

    df_int = df_int_.copy(deep=True)

    # create all three df types
    df_snoleft = df_int[(df_int.single_id1.isin(snolist))
                        & ~(df_int.single_id2.isin(snolist))]
    df_both_sno = df_int[(df_int.single_id1.isin(snolist))
                         & (df_int.single_id2.isin(snolist))]
    df_snoright = df_int[(df_int.single_id2.isin(snolist))
                         & ~(df_int.single_id1.isin(snolist))]
    del df_int

    def switch_side(df):
        # Put all snoRNA on the same side to simplify further analysis
        df.columns = df.columns.str.replace('2', 'old2')
        df.columns = df.columns.str.replace('1', '2')
        df.columns = df.columns.str.replace('old2', '1')
        df.columns = df.columns.str.replace(
            'left', 'oldleft')  # nb of reads aligned to each side
        df.columns = df.columns.str.replace('right', 'left')
        df.columns = df.columns.str.replace('oldleft', 'right')

        return df

    def alpha_order(row):
        id_list = sorted([row.single_id1, row.single_id2])
        if id_list[0] == row.single_id1:
            return False
        return True

    # process df with sno on the right side
    df_snoright = switch_side(df_snoright)

    # process df with sno on both side
    df_both_sno['alpha_order'] = df_both_sno.apply(alpha_order, axis=1)

    df_both_left = df_both_sno.loc[~(df_both_sno.alpha_order)]
    df_both_left = df_both_left.drop(['alpha_order'], axis=1)

    df_both_right = df_both_sno.loc[df_both_sno.alpha_order]
    df_both_right = df_both_right.drop(['alpha_order'], axis=1)
    df_both_right = switch_side(df_both_right)

    df_int = pd.concat([df_snoleft, df_snoright, df_both_left, df_both_right],
                       ignore_index=True, sort=False)

    return df_int


def find_multiple_merge_ids(updated_df):

    counter = Counter(updated_df.merge_ids)
    multiples = [mi for (mi, num) in counter.most_common() if num > 1]
    multiples_df = updated_df.loc[updated_df.merge_ids.isin(multiples)]
    return multiples_df


def create_merged_list(test_list):

    columns = list(set(np.array(test_list).ravel()))

    # Create the matrix
    matrix = np.zeros((len(test_list), len(columns)), dtype=np.bool)
    for idx, row in enumerate(test_list):
        for DG in row:
            matrix[idx, columns.index(DG)] = True

    # Combine the rows
    for i in range(len(columns)):
        col = matrix[:, i]
        if sum(col) > 1:
            matrix[col, :] = np.any(matrix[col, :], axis=0)

    matrix = np.unique(matrix, axis=0)

    # Transform back in list
    master_list = []
    for row in matrix:
        tmp = []
        for idx, col in enumerate(row):
            if col:
                tmp.append(columns[idx])
        master_list.append(tmp)

    return master_list, columns


def find_rows_to_merge(multiples_df):

    def pbt(df):

        first = bt.from_dataframe(df)
        intersect = first.intersect(first, wo=True, s=True, sorted=False)
        new_cols = ['chr1', 'start1', 'end1', 'merge_ids1', 'DG1', 'strand',
                    'chr2', 'start2', 'end2', 'merge_ids2', 'DG2', 'strand2', 'overlap']
        intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                              dtype={'chr': str, 'chr2': str})
        intersect_df = intersect_df.loc[intersect_df.DG1 != intersect_df.DG2]
        intersect_df = intersect_df.loc[(intersect_df.merge_ids1 ==
                                         intersect_df.merge_ids2)]

        intersect_df['tuple'] = [
            (dg1, dg2) # CHANGED !!
            # tuple(sorted(dg1, dg2))
            for (dg1, dg2) in intersect_df[['DG1', 'DG2']].values
        ]
        int_groupby = intersect_df[[
            'merge_ids1', 'tuple'
        ]].groupby('merge_ids1')['tuple'].apply(set).reset_index()
        merged_dg_dict = dict(zip(int_groupby.merge_ids1,
                                  int_groupby['tuple']))
        return merged_dg_dict

    def extra(merged_dg_dict):

        merged_dict = {}
        for ids, the_set in merged_dg_dict.items():
            if len(the_set) == 0:
                continue

            merged_list, _ = create_merged_list(list(the_set))

            tmp_list = []
            for dg_list in merged_list:
                if len(dg_list) == 2:
                    tmp_list.append(tuple(dg_list))
                    tmp_list.append((dg_list[1], dg_list[0]))
                else:
                    perm = permutations(dg_list, 2)
                    for p in perm:
                        tmp_list.append(p)
            merged_dict[ids] = set(tmp_list)
        return merged_dict

    def validate(extra_dict, initial_dict):

        master_dict = {}
        for ids, couples in extra_dict.items():
            if ids in initial_dict:
                tmp = couples.intersection(initial_dict[ids])
                tmp = set([tuple(sorted([x,y])) for (x,y) in tmp])
                master_dict[ids] = tmp
        final_dict = extra(master_dict)

        return final_dict

    left = multiples_df[['chr1', 'start1', 'end1', 'merge_ids',
                         'DG', 'strand1']]
    right = multiples_df[['chr2', 'start2', 'end2', 'merge_ids',
                          'DG', 'strand2']]

    left_dg_dict = pbt(left)
    right_dg_dict = pbt(right)

    extrapol_left_dg_dict = extra(left_dg_dict)
    extrapol_right_dg_dict = extra(right_dg_dict)

    validated_lis = validate(left_dg_dict, extrapol_right_dg_dict)
    validated_ris = validate(right_dg_dict, extrapol_left_dg_dict)

    both_intersection = {}
    for ids in validated_lis.keys():
        if ids in validated_ris:
            both_intersection[ids] = validated_lis[ids].intersection(validated_ris[ids])

    return both_intersection


def more_than_two_merge(final_merge, df):
    """ viz function """
    test = [(x, len(x)) for x in final_merge if len(x) > 2]
    test = sorted(test, key=lambda x: x[1], reverse=True)
    print(len(final_merge) - len(test), len(test))
    for dg, num in test:
        print(num, dg)


def merge_and_return(df, to_merge, multiple_dg_list):

    def merge(df, merge_list):

        master_list = []
        for dg_list in merge_list:
            tmp_df = df.loc[df.DG.isin(dg_list)]
            harm_score = sum(tmp_df.support * tmp_df.harmonic_score) / sum(tmp_df.support)
            geo_score = sum(tmp_df.support * tmp_df.geometric_score) / sum(tmp_df.support)
            dg = '|'.join(sorted(list(set(tmp_df.DG))))
            exp = '|'.join(sorted(list(set(tmp_df.exp))))

            if len(set(tmp_df.single_id1)) != 1 or len(set(tmp_df.single_id2)) != 1:
                print(tmp_df[['start1', 'end1', 'gene_name1', 'start2', 'end2', 'gene_name2']])

            tmp_list = [
                tmp_df.chr1.iloc[0],
                min(tmp_df.start1),
                max(tmp_df.end1),
                tmp_df.chr2.iloc[0],
                min(tmp_df.start2),
                max(tmp_df.end2),
                sum(tmp_df.support),
                sum(tmp_df.left),
                sum(tmp_df.right),
                harm_score,
                geo_score,
                tmp_df.strand1.iloc[0],
                tmp_df.strand2.iloc[0],
                dg,
                tmp_df.gene_id1.iloc[0],
                tmp_df.gene_name1.iloc[0],
                tmp_df.gene_id2.iloc[0],
                tmp_df.gene_name2.iloc[0],
                tmp_df.single_id1.iloc[0],
                tmp_df.single_id2.iloc[0],
                tmp_df.gene_biotype1.iloc[0],
                tmp_df.gene_biotype2.iloc[0],
                exp,
                tmp_df.merge_ids.iloc[0],
            ]

            master_list.append(tmp_list)

        return master_list

    to_merge_df = df.loc[df.DG.isin(multiple_dg_list)]
    original = df.loc[~(df.DG.isin(multiple_dg_list))]
    col_names = df.columns.values

    merged_list = merge(to_merge_df, to_merge)
    merged_df = pd.DataFrame(data=merged_list, columns=col_names)

    final_df = pd.concat([original, merged_df], ignore_index=True, sort=False)
    final_df.sort_values(['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2'],
                          inplace=True)

    return final_df


# def get_single_names(df, file):
#
#     ref_df = pd.read_csv(file, sep='\t')
#
#     df['name1'] = df.single_id1.map(dict(zip(ref_df.gene_id, ref_df.gene_name)))
#     df['name2'] = df.single_id2.map(dict(zip(ref_df.gene_id, ref_df.gene_name)))
#
#     df['biotype1'] = df.single_id1.map(dict(zip(ref_df.gene_id, ref_df.gene_biotype)))
#     df['biotype2'] = df.single_id2.map(dict(zip(ref_df.gene_id, ref_df.gene_biotype)))
#
#     return df


def write_df(df):

    df.to_csv(out_file, sep='\t', index=None)


def main():

    data_df = load_df(in_file)
    snolist = get_snoRNA_list(gene_bed_biotype)

    updated_df = snoRNA_side(data_df, snolist)

    multiples_df = find_multiple_merge_ids(updated_df)
    to_merge_raw = find_rows_to_merge(multiples_df)

    final_to_merge = []
    more_than_one = []
    for ids, the_set in to_merge_raw.items():
        if len(the_set) != 0:
            tmp_to_merge, tmp_more_than_one = create_merged_list(list(the_set))
            final_to_merge += tmp_to_merge
            more_than_one += tmp_more_than_one

    # Just to get some viz on the DG to merge
    # more_than_two_merge(final_to_merge, updated_df)

    final_df = merge_and_return(updated_df, final_to_merge, more_than_one)
    final_df = final_df[[
        'chr1',
        'start1',
        'end1',
        'strand1',
        'chr2',
        'start2',
        'end2',
        'strand2',
        'support',
        'DG',
        'exp',
        'single_id1',
        'single_id2',
        'gene_name1',
        'gene_name2',
        'gene_biotype1',
        'gene_biotype2',
    ]]
    final_df.rename(columns={'gene_biotype1': 'biotype1',
                             'gene_biotype2': 'biotype2',
                             'gene_name1': 'name1',
                             'gene_name2': 'name2'}, inplace=True)

    # final_df = get_single_names(final_df, gene_bed_biotype)

    write_df(final_df)


if __name__ == '__main__':

    main()
