import os

import numpy as np
import pandas as pd

raw_file = snakemake.input.initial_file
gene_bed_biotype = snakemake.input.gene_bed_biotype
curated_multimap_file = snakemake.input.man_curated

multiId_file = snakemake.params.multi
out_file = snakemake.output.single_id


# header of the raw file
"""['chr1' 'start1' 'end1' 'chr2' 'start2' 'end2' 'support' 'left' 'right'
 'harmonic_score' 'geometric_score' 'strand1' 'strand2' 'DG' 'gene_id1'
 'gene_name1' 'gene_id2' 'gene_name2' 'exp']"""


def get_multimap_genes(df, ref_df):


    biotype_dict = dict(zip(ref_df.gene_id, ref_df.gene_biotype))

    cols = ['gene_id1', 'gene_name1', 'gene_id2',  'gene_name2']
    multiId = set()
    multi_to_single = {}
    for matrix in [df[cols].values[:,:2], df[cols].values[:,2:]]:
        for id, name in matrix:
            if ';' in id:
                multi_to_single[id] = ''
                biotypes = []

                # Check if there is rRNA
                for single_id in id.split(';'):
                    single_biotype = biotype_dict[single_id]
                    biotypes.append(single_biotype)
                    if single_biotype == 'rRNA':
                        multi_to_single[id] = single_id
                        break;

                # Check if the most likely biotypes are in biotypes
                for b in ['snoRNA', 'tRNA', 'snRNA']:
                    if biotypes.count(b) == 1:
                        idx = biotypes.index(b)
                        multi_to_single[id] = id.split(';')[idx]
                        break
                else:
                    if biotypes.count('lncRNA') == 2 and 'SNHG' in name:
                        for idx, single_name in enumerate(name.split(';')):
                            if 'SNGH' in single_name:
                                multi_to_single[id] = id.split(';')[idx]
                        multi_to_single[id] = id.split(';')[idx]
                    elif biotypes.count('protein_coding') == 1:
                        idx = biotypes.index('protein_coding')
                        multi_to_single[id] = id.split(';')[idx]

                # Just for visualization
                if multi_to_single[id] == '' and (id, name, ';'.join(biotypes)) not in multiId:
                    print((id, name, ';'.join(biotypes)), '->' + multi_to_single[id])

                multiId.add((id, name, ';'.join(biotypes)))

    print(len(multiId))
    print(len(set(multi_to_single.values())))

    # with open(multiId_file, 'w') as f:
    #     for multi, single in multi_to_single.items():
    #         f.write('{},{}\n'.format(multi, single))



def create_single_id(df, ref_df):

    curated_ids_df = pd.read_csv(curated_multimap_file,
                                names=['mmap', 'uniq_id'])

    currated_dict = dict(zip(curated_ids_df.mmap,
                           curated_ids_df.uniq_id))

    def map_single(row, gene_num):
        gene_id = row['gene_id{}'.format(gene_num)]
        if gene_id in currated_dict:
            return currated_dict[gene_id]
        elif 'intergenic' in gene_id:
            chr = row['chr{}'.format(gene_num)]
            start = str(row['start{}'.format(gene_num)])
            return 'intergenic_{}_{}'.format(chr, start[:3])
        else:
            return gene_id

    new_df = df.copy(deep=True)
    new_df['single_id1'] = new_df.apply(map_single, gene_num='1', axis=1)
    new_df['single_id2'] = new_df.apply(map_single, gene_num='2', axis=1)

    id_name_dict = dict(zip(ref_df.gene_id, ref_df.gene_name))
    id_biotype_dict = dict(zip(ref_df.gene_id, ref_df.gene_biotype))

    # Adding the single gene_name and gene_biotype
    new_df['gene_name1'] =  new_df['single_id1'].map(id_name_dict)
    new_df['gene_name2'] =  new_df['single_id2'].map(id_name_dict)
    new_df['gene_biotype1'] =  new_df['single_id1'].map(id_biotype_dict)
    new_df['gene_biotype2'] =  new_df['single_id2'].map(id_biotype_dict)
    return new_df


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


def keep_best_overlap(DG_mult, DG_single, data_df_):

    def get_length(df, side):
        length = [
            df.at[x, f'end{side}'] - df.at[x, f'start{side}']
            for x in range(len(df))
        ]
        return length

    main_df = data_df_.copy(deep=True)
    main_df = main_df.loc[main_df.DG.isin(DG_single[:,0])]

    data_df = data_df_.copy(deep=True)
    data_df = data_df.loc[~(data_df.DG.isin(DG_single[:,0]))]

    counts = 0
    for dg, num in DG_mult:
        tmp_df = data_df.loc[data_df.DG == dg]
        tmp_df.reset_index(drop=True, inplace=True)

        lengths = list(zip(get_length(tmp_df, 1),
                            get_length(tmp_df, 2)))
        biotypes = np.array([
            (tmp_df.at[x, 'gene_biotype1'], tmp_df.at[x, 'gene_biotype2'])
            for x in range(len(tmp_df))
        ])

        if len(tmp_df) == 2:
            if abs(lengths[0][0] - lengths[1][0]) > 10:
                if 'snoRNA' in biotypes[:,0] and 'snoRNA' not in biotypes[:,1]:
                    counts+=1
                    print(f'DG: {dg}', end=' ---------\n')
                    print(tmp_df[['start1', 'end1', 'gene_name1',
                                'start2', 'end2', 'gene_name2']])
                    print(biotypes)

            # elif abs(lengths[0][1] - lengths[1][1]) > 10:


    print('total counts:', counts)







        # max_sum_interval = 0
        # for idx, i in enumerate(tmp_df.index):
        #     length1 = tmp_df.at[i, 'end1'] - tmp_df.at[i, 'start1']
        #     length2 = tmp_df.at[i, 'end2'] - tmp_df.at[i, 'start2']
        #
        #     sum_interval = length1 + length2
        #
        #     if not max_sum_interval or sum_interval > max_sum_interval:
        #         max_sum_interval = sum_interval
        #         good_index = idx
        #
        # main_df = main_df.append(tmp_df.iloc[good_index], ignore_index=True)

    return main_df


def write_filered_data(filtered_df):

    filtered_df.to_csv(out_file, index=None)


def main():

    data_df = pd.read_csv(raw_file)
    data_df = data_df.fillna('intergenic')

    annotation_df = pd.read_csv(gene_bed_biotype, sep='\t')

    # get a list of all coordinates assigned to muliples genes
    # get_multimap_genes(data_df, annotation_df)

    single_id_data = create_single_id(data_df, annotation_df)

    DG_multiple_list, DG_single_list = get_multiple_line_DG(single_id_data)

    filtered_df = keep_best_overlap(DG_multiple_list,
                                       DG_single_list,
                                       single_id_data)

    # write_filered_data(filtered_df)


if __name__ == '__main__':
    main()
