import os
from collections import defaultdict
import ast

import numpy as np
import pandas as pd

data_file = snakemake.input.full_merge
ref_file = snakemake.input.gene_bed_biotype
snoDB_file = snakemake.input.snodb

out_interactions = snakemake.output.interactions
out_nodes = snakemake.output.nodes
out_edges = snakemake.output.edges
out_clip_only = snakemake.output.mapped_clip_only


def load_ref_dict():

    df = pd.read_csv(ref_file, sep='\t')
    dict_ = dict(zip(df.gene_name, df.gene_id))
    dict_['AARS'] = 'ENSG00000090861'
    dict_['TROVE2'] = 'ENSG00000116747'
    return dict_

def load_snodb():

    snodb_df = pd.read_csv(snoDB_file, sep='\t')
    snodb_df.rename(columns={'gene_id_annot2020': 'gene_id',
                             'gene_name_annot2020': 'gene_name'},
                             inplace=True)
    return snodb_df

def load_df():

    df = pd.read_csv(data_file, sep='\t')
    df.sort_values(['single_id1', 'single_id2'], inplace=True)

    return df


def remove_dups(df_):
    """ Simple function to remove sno-sno duplicates keeping only the ones
        that are in alphabetical order """
    df = df_.copy(deep=True)
    df['merge_ids'] = [
        '_'.join(sorted([df.at[i, 'single_id1'], df.at[i, 'single_id2']]))
        for i in df.index
    ]
    df.drop_duplicates(subset=['merge_ids', 'DG'], inplace=True)

    return df


def write_mapped_clip_data(full_clip, df_cols):

    full_clip_df = pd.DataFrame(full_clip, columns=['prot', 'sno_score', 'sno_pValue',
                                                   'target_score', 'target_pValue'] +
                                                   df_cols)

    full_clip_df.sort_values(by=['sno_pValue', 'target_pValue'], ascending=False, inplace=True)
    full_clip_df = full_clip_df.loc[full_clip_df.biotype2 == 'protein_coding']
    full_clip_df.drop(columns=['single_id1', 'single_id2', 'prot_info', 'merge_ids'], inplace=True)

    full_clip_df.to_csv(out_clip_only, sep='\t', index=False)


def create_sif(df_, id_name_dict):
    """ Function to create a siff file and to extract at the same time
        the clip data """

    df = df_.copy(deep=True)
    df = df[[
        'chr2', 'start2', 'end2','support', 'DG', 'single_id1', 'single_id2',
        'name1', 'name2', 'biotype1', 'biotype2', 'host_interaction',
        'E', 'prot_info', 'merge_ids'
    ]]
    couter_dict = defaultdict(int)

    siff_list = []
    interaction_list = []
    single_int_clip = []
    full_clip = []

    for col in df.values:
        chr, start, end, sup, DG, id1, id2, name1, name2, bio1, bio2, h_i, E, p_i, merged_id = col
        pro_edge = 'pr'

        if h_i:
            edge = 'host'
        else:
            edge = 'rr'

        if couter_dict[merged_id]:
            edge += str(couter_dict[merged_id]+1)

        couter_dict[merged_id] += 1
        siff_list.append((id1, edge, id2))
        interaction_list.append(edge)

        # data for the clip
        prot_info = str(p_i)
        if prot_info != 'nan':
            for prot_list in prot_info.replace(('['), '').replace(']', '').replace('),', ')),').split('), '):
                tup = ast.literal_eval(prot_list)
                prot, sno_score, sno_pValue, target_score, target_pValue = tup
                full_clip.append([*tup, *col])

                # Add value for the snoRNA
                tmp = (id_name_dict[prot] + '_prot', prot, pro_edge, id1,
                       name1, bio1, sno_score, sno_pValue, DG)
                single_int_clip.append(tmp)

                # Add value fot the target
                tmp = (id_name_dict[prot] + '_prot', prot, pro_edge, id2,
                       name2, bio2, target_score, target_pValue, DG)
                single_int_clip.append(tmp)


    clip_df = pd.DataFrame(single_int_clip, columns=['prot_id', 'prot_name', 'interaction',
                                               'target_id', 'target_name', 'target_biotype',
                                               'target_score', 'target_pValue', 'DG'])

    # Add the prot data to the siff list
    siff_list += set(zip(clip_df.prot_id, clip_df.interaction, clip_df.target_id))

    siff_df = pd.DataFrame(siff_list, columns=['source', 'interaction', 'target'])
    siff_df.to_csv(out_interactions, sep='\t', index=False)

    write_mapped_clip_data(full_clip, list(df.columns))

    df['interaction_type'] = interaction_list

    return df, clip_df


def create_node_info(df, clip_df, snodb_df):

    new_cols = ['gene_id', 'gene_name', 'biotype']

    def transform_df(df, col1, col2, col3):
        tmp_df = df.copy(deep=True)[[col1, col2, col3]]
        tmp_df.drop_duplicates(inplace=True)
        tmp_df.columns = new_cols
        return tmp_df

    node_df_left = transform_df(df, 'single_id1', 'name1', 'biotype1')
    node_df_right = transform_df(df, 'single_id2', 'name2', 'biotype2')
    clip_prot = clip_df.copy(deep=True)
    clip_prot['biotype'] = 'protein'
    clip_prot = transform_df(clip_prot, 'prot_id', 'prot_name', 'biotype')

    node_df = pd.concat([node_df_left, node_df_right, clip_prot],
                        axis=0, ignore_index=True)
    node_df.drop_duplicates(inplace=True)

    box_type_dict = dict(zip(snodb_df.gene_id, snodb_df['box type']))
    node_df['box_type'] = node_df.gene_id.map(box_type_dict)

    node_df.to_csv(out_nodes , sep='\t', index=False)


def create_edges_info(df_, clip_df_):

    df = df_.copy(deep=True)
    clip_df = clip_df_.copy(deep=True)
    clip_df.columns = [
        'single_id1', 'name1', 'interaction_type', 'single_id2',
        'name2', 'biotype2', 'target_score', 'target_pValue', 'DG'
    ]
    df = df.append(clip_df, sort=False)

    df['shared_name'] = df.single_id1 + ' (' + df.interaction_type + \
                            ') ' + df.single_id2
    df = df[['shared_name', 'support', 'name1', 'name2', 'biotype2',
            'host_interaction', 'E', 'target_score', 'target_pValue']]
    df.columns = ['shared_name', 'support', 'gene_name1', 'gene_name2',
                 'target_biotype', 'host_interaction', 'E', 'target_score',
                 'target_pValue']

    df.to_csv(out_edges, sep='\t', index=False)


def main():

    id_name_dict = load_ref_dict()

    snodb_df = load_snodb()

    base_df = load_df()
    base_df = remove_dups(base_df)

    int_df, clip_df = create_sif(base_df, id_name_dict)

    create_node_info(base_df, clip_df, snodb_df)

    create_edges_info(int_df, clip_df)



if __name__ == '__main__':
    main()
