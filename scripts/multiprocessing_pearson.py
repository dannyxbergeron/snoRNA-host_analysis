from scipy import stats
import pandas as pd
import numpy as np

import multiprocessing
from itertools import repeat

in_file = snakemake.input.tpm
ref_file = snakemake.input.parsed_gtf

out_file = snakemake.output.pearson_corr

CORES = snakemake.threads


def load_df(file, ref_file):

    ref_df = pd.read_csv(ref_file, sep='\t')
    id_biotype_dict = dict(zip(ref_df.gene_id, ref_df.gene_biotype))

    df = pd.read_csv(file, sep='\t')
    df = df[[
        'gene_id', 'gene_name',
        'HCT116_1', 'HCT116_2',
        'INOF_1',
        'MCF7_1', 'MCF7_2',
        'PC3_1', 'PC3_2',
        'SKOV_frg_1', 'SKOV_frg_2',
        # 'SKOV_nf_1', 'SKOV_nf_2',
        'TOV112D_1', 'TOV112D_2',
        'Liver_1', 'Liver_2', 'Liver_3',
        'Breast_1', 'Breast_2', 'Breast_3', 'Breast_4',
        'Testis_1', 'Testis_2', 'Testis_3',
        'Prostate_1', 'Prostate_2', 'Prostate_3',
        'Ovary_1', 'Ovary_2', 'Ovary_3',
        'SkeletalMuscle_1', 'SkeletalMuscle_2', 'SkeletalMuscle_3',
        'Brain_1', 'Brain_2', 'Brain_3',
        'humanRef_1', 'humanRef_2', 'humanRef_3',
        'brainLam_1', 'brainLam_2', 'brainLam_3'
    ]]
    df['gene_biotype'] = df.gene_id.map(id_biotype_dict)

    return df


def process(df):

    new_df = df.copy(deep=True)
    new_df['avg'] = new_df[new_df.columns[2:-1]].mean(axis=1)
    new_df = new_df.loc[new_df.avg > 1]
    new_df.index = new_df.gene_id
    new_df = new_df.drop(['gene_id', 'gene_name', 'gene_biotype', 'avg'], axis=1)
    return new_df


def pearson(sno_cols, sno_names, other_cols, other_names):
    res = []
    for id, sno_col in zip(sno_names, sno_cols):
        for o_id, o_col in zip(other_names, other_cols):
            if id != o_id:
                coef, p_value = stats.pearsonr(sno_col, o_col)
                # CHANGED !!!!!
                res.append([id, o_id, coef, p_value])
                # if p_value < 0.05:
                #     res.append([id, o_id, coef, p_value])

    return res


def main():

    df = load_df(in_file, ref_file)

    sno_df = df.loc[df.gene_biotype == 'snoRNA']
    sno_df = process(sno_df)

    others = process(df)

    with multiprocessing.Pool(processes=CORES) as pool:
        res = pool.starmap(pearson, zip(np.array_split(sno_df.values, CORES),
                                        np.array_split(sno_df.index, CORES),
                                        repeat(others.values),
                                        repeat(others.index)))

    master_res = []
    for mat in res:
        master_res += mat

    results = pd.DataFrame(master_res, columns=['sno_id', 'other_id', 'coef', 'p_value'])

    results['sno_name'] = results.sno_id.map(dict(zip(df.gene_id, df.gene_name)))
    results['other_name'] = results.other_id.map(dict(zip(df.gene_id, df.gene_name)))
    results['other_biotype'] = results.other_id.map(dict(zip(df.gene_id, df.gene_biotype)))

    results = results.sort_values('p_value')
    results.to_csv(out_file, sep='\t', index=None)


if __name__ == '__main__':

    main()
