import os

from collections import defaultdict

import numpy as np
import pandas as pd


multiDG_file = snakemake.input.single_id

out_file = snakemake.output.single_id_and_dg

MIN_LENGTH = snakemake.params.min_length
HOST_MAX_OFF = snakemake.params.host_max_offset

MASTER_COUNT = 0


def get_multiple_line_DG(dataframe):

    DG_line_counts = dataframe.groupby('DG')[['chr1']].count().sort_values(
        by='chr1', axis=0, ascending=False)

    DG_all = list(zip(
        DG_line_counts.index,
        DG_line_counts.chr1
    ))
    DG_multiple_list = [(x, y) for (x,y) in DG_all if y > 1]
    DG_single_list = [(x, y) for (x,y) in DG_all if y == 1]
    print('Number of multiple DG: {}'.format(len(DG_multiple_list)))
    print('Number of single DG: {}'.format(len(DG_single_list)))

    return np.array(DG_multiple_list), np.array(DG_single_list)


def deal_with_two(tmp_df, DG):

    inverse = {1: 2, 2: 1}

    biotypes = tmp_df[['gene_biotype1', 'gene_biotype2']].values.ravel()
    ids = tmp_df[['single_id1', 'single_id2']].values

    # Deal with miRNA embeded in snoRNAs
    if np.count_nonzero(biotypes == 'miRNA') == 1: # 82
        if np.count_nonzero(biotypes.reshape(2, 2)[:,0] == 'miRNA') == 1: # 30
            if 'snoRNA' in list(tmp_df.gene_biotype1):
                return list(tmp_df.loc[tmp_df.gene_biotype1 == 'snoRNA'].values[0])
            else:
                return None
        else: # 52
            return list(tmp_df.loc[tmp_df.gene_biotype2 == 'snoRNA'].values[0])

    def deal_side(side):

        tmp = tmp_df.copy(deep=True)
        lengths = list(tmp[f'end{side}'] - tmp[f'start{side}'])
        biotypes_side = list(tmp[f'gene_biotype{side}'])
        biotypes_other_side = set(tmp['gene_biotype{}'.format(inverse[side])])

        # If the two overlapping genes are the targets
        if np.count_nonzero(tmp_df[f'gene_biotype{side}'].values == 'snoRNA') == 0: # 29
            # deal with the two small gene overlap
            if lengths[0] < MIN_LENGTH: return list(tmp.iloc[1]) # 6
            elif lengths[1] < MIN_LENGTH: return list(tmp.iloc[0]) # 13
            # deal with the other and take the bigger interval and add
            # the offset to the dictionnaty for further use
            if lengths[0] > lengths[1] and not pd.isnull(tmp.at[0, f'gene_name{side}']):
                idx = 0
            else: idx = 1
            tmp[f'start{side}'] = tmp[f'start{side}'].min()
            tmp[f'end{side}'] = tmp[f'end{side}'].max()
            return list(tmp.iloc[idx])

        # If the snoRNA is in the overlapping.
        elif np.count_nonzero(tmp_df[f'gene_biotype{side}'].values == 'snoRNA') == 1: # 964
            # I still included the ones that are really
            # longer than the snoRNA because otherwise the information is lost since there
            # is no snoRNA on one of the side. Could be filter out later.
            if biotypes_other_side != {'snoRNA',}: # 556
                if biotypes_side[0] == 'snoRNA' and lengths[0] > MIN_LENGTH: idx = 0
                elif biotypes_side[1] == 'snoRNA' and lengths[1] > MIN_LENGTH: idx = 1
                else:
                    return None
                tmp[f'start{side}'] = tmp[f'start{side}'].min()
                tmp[f'end{side}'] = tmp[f'end{side}'].max()
                return list(tmp.iloc[idx])

            # The ones that are overlapping in the target and is a sno-other
            # But that on the other side it's all snoRNA
            else: # 408
                # Remove the too small intervals
                if max(lengths) < MIN_LENGTH:
                    return None

                if biotypes_side[0] == 'snoRNA': idx = 0
                else: idx = 1
                other_idx = 0 if idx == 1 else 0

                # If the host gene is significantly greater than the snoRNA
                if lengths[idx] < lengths[other_idx]: # 17
                    new_idx = other_idx
                    return list(tmp.iloc[new_idx])
                # If the snoRNA interval is pretty much the same as the host one
                else: # 384
                    tmp[f'start{side}'] = tmp[f'start{side}'].min()
                    tmp[f'end{side}'] = tmp[f'end{side}'].max()
                    return list(tmp.iloc[idx])


    if ids[0][0] == ids[1][0]: # 581
        return deal_side(2)
    elif ids[0][1] == ids[1][1]: # 412
        return deal_side(1)


def deal_with_more(tmp_df_, DG):

    tmp_df = tmp_df_.copy(deep=True)

    inverse = {1: 2, 2: 1}

    def deal_side(side):
        tmp = tmp_df.copy(deep=True)
        lengths = list(tmp[f'end{side}'] - tmp[f'start{side}'])
        max_length = max(lengths)
        biotypes_side = list(tmp[f'gene_biotype{side}'])
        biotype_other_side = set(tmp['gene_biotype{}'.format(inverse[side])])

        tmp[f'start{side}'] = tmp[f'start{side}'].min()
        tmp[f'end{side}'] = tmp[f'end{side}'].max()

        if biotype_other_side == {'snoRNA',}:
            if biotypes_side.count('snoRNA') > 1:
                idx = -1
                for i, b in enumerate(biotypes_side):
                    if b != 'snoRNA':
                        idx = i
                        break
                return list(tmp.iloc[idx])
            elif biotypes_side.count('snoRNA') == 1: # 0
                idx = biotypes_side.index('snoRNA')
                return list(tmp.iloc[idx])
            else:
                idx = -1
                for i, l in enumerate(lengths):
                    if l == max_length:
                        idx = i
                return list(tmp.iloc[idx])
        else:
            if biotypes_side.count('snoRNA') > 1:
                return None
            idx = biotypes_side.index('snoRNA')
            return list(tmp.iloc[idx])

    def helper_four(tmp, side):
        df = tmp.copy(deep=True)
        df = df[[f'start{side}', f'end{side}', f'single_id{side}', f'gene_biotype{side}']]
        df.drop_duplicates(inplace=True)
        lengths = list(df[f'end{side}'] - df[f'start{side}'])
        sno_idx = list(df[f'gene_biotype{side}']).index('snoRNA')
        other_idx = 1 if sno_idx == 0 else 0
        if lengths[sno_idx] >= lengths[other_idx]:
            tmp_ = tmp.loc[tmp[f'gene_biotype{side}'] == 'snoRNA']
            return list(tmp_.index)
        else:
            tmp_ = tmp.loc[tmp[f'gene_biotype{side}'] != 'snoRNA']
            return list(tmp_.index)

    def deal_four():
        tmp = tmp_df.copy(deep=True)
        biotypes = tmp[['gene_biotype1', 'gene_biotype2']].values
        biotypes_left = biotypes[:,0]
        biotypes_right = biotypes[:,1]


        if np.count_nonzero(biotypes_left == 'snoRNA') == 2:
            if np.count_nonzero(biotypes_right == 'snoRNA') == 2:
                left_indices = helper_four(tmp, 1)
                tmp[f'start1'] = tmp[f'start1'].min()
                tmp[f'end1'] = tmp[f'end1'].max()
                tmp = tmp.loc[tmp.index.isin(left_indices)]

                right_indices = helper_four(tmp, 2)
                right_index = right_indices[0]

                tmp[f'start2'] = tmp[f'start2'].min()
                tmp[f'end2'] = tmp[f'end2'].max()
                tmp = tmp.loc[tmp.index == right_index]
                return list(tmp.iloc[0])
            else:
                left_indices = helper_four(tmp, 1)
                tmp[f'start1'] = tmp[f'start1'].min()
                tmp[f'end1'] = tmp[f'end1'].max()
                tmp = tmp.loc[tmp.index.isin(left_indices)]

                lengths = list(tmp[f'end2'] - tmp[f'start2'])
                idx = -1
                max_l = 0
                for i, l in enumerate(lengths):
                    if l > max_l:
                        idx = i
                        max_l = l
                return list(tmp.values[idx])
        elif np.count_nonzero(biotypes_right == 'snoRNA') == 2:
            right_indices = helper_four(tmp, 2)
            tmp[f'start2'] = tmp[f'start2'].min()
            tmp[f'end2'] = tmp[f'end2'].max()
            tmp = tmp.loc[tmp.index.isin(right_indices)]

            lengths = list(tmp[f'end1'] - tmp[f'start1'])

            biotypes = list(tmp[f'gene_biotype1'])
            if 'tRNA' in biotypes or 'scaRNA' in biotypes:
                for bio in ['tRNA', 'scaRNA']:
                    if bio in biotypes:
                        idx = biotypes.index(bio)
                tmp[f'start1'] = tmp[f'start1'].min()
                tmp[f'end1'] = tmp[f'end1'].max()
                return list(tmp.values[idx])

            idx = -1
            max_l = 0
            for i, l in enumerate(lengths):
                if l > max_l:
                    idx = i
                    max_l = l
            return list(tmp.values[idx])


    ids = tmp_df[['single_id1', 'single_id2']].values

    # Deal with the one that has only one id on one side
    if len(set(ids[:,0])) == 1: # 3
        return deal_side(2)
    elif len(set(ids[:,1])) == 1: # 7
        return deal_side(1)

    # Deal with the ones that have 4 entries per DG
    if len(tmp_df) == 4: # 44
        return deal_four()

    # Deal with the nonsens rest... # 6
    def deal_nonsense(tmp, side):
        lengths = list(tmp_df[f'end{side}'] - tmp_df[f'start{side}'])
        max_length = max(lengths)
        indices = set([i for (i,x) in enumerate(lengths) if x == max_length])
        return indices

    idx_left = deal_nonsense(tmp_df, 1)
    idx_right = deal_nonsense(tmp_df, 2)

    idx = list(idx_left.intersection(idx_right))[0]

    tmp_df[f'start1'] = tmp_df[f'start1'].min()
    tmp_df[f'end1'] = tmp_df[f'end1'].max()

    tmp_df[f'start2'] = tmp_df[f'start2'].min()
    tmp_df[f'end2'] = tmp_df[f'end2'].max()

    if 'snoRNA' in list(tmp_df.values[idx]):
        return list(tmp_df.values[idx])
    else:
        return None



def keep_best_overlap(DG_mult, DG_single, data_df_):

    main_df = data_df_.copy(deep=True)
    main_df = main_df.loc[main_df.DG.isin(DG_single[:,0])]

    data_df = data_df_.copy(deep=True)
    data_df = data_df.loc[~(data_df.DG.isin(DG_single[:,0]))]

    best_overlap = []
    for dg, num in DG_mult:
        tmp_df = data_df.loc[data_df.DG == dg]
        tmp_df.reset_index(drop=True, inplace=True)

        # deal with 2 row per DG
        if len(tmp_df) == 2: # 1075 in total
            DG_info = deal_with_two(tmp_df, dg)
            if DG_info:
                best_overlap.append(DG_info)

        # deal with more than 2 row per DG
        else: # 60
            DG_info = deal_with_more(tmp_df, dg)
            if DG_info:
                best_overlap.append(DG_info)


    best_overlap_df = pd.DataFrame(best_overlap, columns=data_df.columns)
    result_df = pd.concat([main_df, best_overlap_df], axis=0)

    return result_df


def write_filered_data(filtered_df):

    filtered_df.to_csv(out_file, index=None)


def main():

    data_df = pd.read_csv(multiDG_file)

    DG_multiple_list, DG_single_list = get_multiple_line_DG(data_df)

    filtered_df = keep_best_overlap(DG_multiple_list,
                                       DG_single_list,
                                       data_df)


    filtered_df['exp'] = [x.split('_')[1] for x in list(filtered_df.DG)]

    print(MASTER_COUNT)

    write_filered_data(filtered_df)


if __name__ == '__main__':
    main()
