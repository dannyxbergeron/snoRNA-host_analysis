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

def write_data(single_id_data):

    single_id_data.to_csv(out_file, index=None)


def main():

    data_df = pd.read_csv(raw_file)
    data_df = data_df.fillna('intergenic')

    annotation_df = pd.read_csv(gene_bed_biotype, sep='\t')

    # get a list of all coordinates assigned to muliples genes
    # get_multimap_genes(data_df, annotation_df)

    single_id_data = create_single_id(data_df, annotation_df)
    write_data(single_id_data)


if __name__ == '__main__':
    main()
