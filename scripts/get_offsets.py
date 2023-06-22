from collections import defaultdict

import pandas as pd
import numpy as np

from pybedtools import BedTool as bt

in_file = snakemake.input.merged_windows
gene_bed_biotype_file = snakemake.input.gene_bed_biotype

out_file = snakemake.output.merged_with_offset

def get_ncRNA_dict(df_, ncRNA_df):

    df = df_.copy(deep=True)

    ncRNA_dict = dict(zip(ncRNA_df.gene_id, zip(ncRNA_df.start, ncRNA_df.end)))

    ncRNA_off_dict = {}
    sno_data_df = df[['start1', 'end1', 'single_id1', 'DG']]
    sno_data_df['side'] = 1

    other_side_sno = df.loc[df.biotype2.isin(['snoRNA', 'tRNA', 'scaRNA'])]
    other_side_sno = other_side_sno[['start2', 'end2', 'single_id2', 'DG']]
    other_side_sno['side'] = 2

    other_side_sno.columns = sno_data_df.columns

    merged_df = pd.concat([sno_data_df, other_side_sno], axis=0)

    master_dict = defaultdict(dict)
    for start,end,id,DG,side in merged_df.values:
        sno_start,sno_end = ncRNA_dict[id]

        off_left = sno_start - start
        if off_left <= 0: off_left = 0

        off_right = end - sno_end
        if off_right <= 0: off_right = 0

        total_offset = off_left + off_right
        master_dict[side][DG] = total_offset

    return master_dict



def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_name1', 'DG', 'strand',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df



def get_other_dict(df_, ncRNA_df, ncRNA_off_dict):

    df = df_.copy(deep=True)
    other_df = df.loc[~(df.biotype2.isin(['snoRNA', 'tRNA', 'scaRNA']))]
    other_df = other_df[['chr2', 'start2', 'end2', 'name2', 'DG', 'strand2']]

    ncRNA_df = ncRNA_df.drop(columns=['gene_biotype'])

    intersect_df = bedtools(other_df, ncRNA_df)
    print(intersect_df)

    gb_intersect_df = intersect_df[['DG', 'overlap']].groupby(['DG']).sum().reset_index()
    for dg, overlap in gb_intersect_df.values:
        ncRNA_off_dict[2][dg] = overlap

    return ncRNA_off_dict


def main():

    data_df = pd.read_csv(in_file, sep='\t')
    ref_df = pd.read_csv(gene_bed_biotype_file, sep='\t')
    ncRNA_df = ref_df.loc[ref_df.gene_biotype.isin(['snoRNA', 'tRNA', 'scaRNA'])]

    ncRNA_off_dict = get_ncRNA_dict(data_df, ncRNA_df)

    full_off_dict = get_other_dict(data_df, ncRNA_df, ncRNA_off_dict)

    data_df['offset1'] = data_df.DG.map(full_off_dict[1])
    data_df['offset2'] = data_df.DG.map(full_off_dict[2])

    data_df['offset2'] = data_df['offset2'].fillna(value=0)

    # for i in data_df.index:
    #     d1 = data_df.loc[data_df.index == i][['chr1', 'start1', 'end1', 'name1',
    #                                             'single_id1', 'strand1']]
    #     d2 = data_df.loc[data_df.index == i][['chr2', 'start2', 'end2', 'name2',
    #                                             'single_id2', 'strand2']]
    #     inter = bedtools(d1, d2)
    #     if len(inter) > 0 and inter.gene_name1.iloc[0] != inter.gene_id2.iloc[0]:
    #         print(inter)

    data_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
