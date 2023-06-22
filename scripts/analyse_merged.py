import os

import numpy as np
import pandas as pd


in_file = snakemake.input.merged_windows
gene_bed_biotype_file = snakemake.input.gene_bed_biotype
sno_df_file = snakemake.input.snoDB

out_file = snakemake.output.merged_double

TRESH = snakemake.params.min_length


def load_and_filter_df(file):

    df = pd.read_csv(file, sep='\t')
    df = df.fillna('intergenic')

    # ------------------Filters-------------------
    # remove less than 9 nt interactions
    f1 = (df.end1 - df.start1 >= TRESH)
    f2 = (df.end2 - df.start2 >= TRESH)
    df = df.loc[f1 & f2]

    # remove everything intergenic
    df = df.loc[df.biotype2 != 'intergenic']
    # remove the intramolecular interactions
    df = df.loc[df.single_id1 != df.single_id2]
    # ------------------Filters-------------------

    return df

def load_refs(snoDB_file, ref_file):

    snodb_df = pd.read_csv(snoDB_file, sep='\t')
    snodb_df.rename(columns={'gene_id_annot2020': 'gene_id',
                             'gene_name_annot2020': 'gene_name'},
                             inplace=True)
    ref_df = pd.read_csv(ref_file, sep='\t')
    return snodb_df, ref_df



def flag_host_interactions(df, snodb_df, ref_df):
    """Create a positive booleen columns for row with sno-host interaction."""
    snodb_df = snodb_df[['gene_id', 'host gene id']]
    snodb_df = snodb_df.dropna()
    snoDB_dict = dict(zip(snodb_df['gene_id'], snodb_df['host gene id']))

    host_col = []
    for i in df.index:
        sno_id = df.at[i, 'single_id1']
        other_id = df.at[i, 'single_id2']
        if sno_id in snoDB_dict and other_id == snoDB_dict[sno_id]:
            host_col.append(True)
        else:
            tmp = ref_df.loc[ref_df.gene_id == other_id].values
            chr, start, end, gene_id, gene_name, strand, gene_biotype = tmp[0]

            sno_chr = df.at[i, 'chr1']
            sno_start = df.at[i, 'start1']
            sno_end = df.at[i, 'end1']

            if chr == sno_chr and start < sno_start and end > sno_end:
                host_col.append(True)
            else:
                host_col.append(False)

    df['host_interaction'] = host_col

    return df


def find_substring(df, substrings):

    def find(row, sub=substrings):
        for s in substrings:
            if s not in row.exp:
                return False
        return True

    new_df = df.copy(deep=True)
    new_df['search'] = new_df.apply(find, sub=substrings, axis=1)
    new_df = new_df.loc[new_df.search]

    length = len(new_df.loc[new_df.search])

    sommation = sum(list(new_df.support))
    return length, sommation


def method_ovelap(df):

    print('number of different interactions: {}'.format(len(df)))
    print('in all three: {}'.format(find_substring(df, ['P', 'L', 'S'])))
    all_df = df.loc[(df.exp.str.contains('L'))
                    & (df.exp.str.contains('S'))]
    print(all_df)
    print()
    print('in Paris-Ligrt-seq: {}'.format(find_substring(df, ['P', 'L'])))
    print('in Paris-Splash: {}'.format(find_substring(df, ['P', 'S'])))
    print('in ligrt-seq-Splash: {}'.format(find_substring(df, ['L', 'S'])))


def double_sno_sno(df_):

    df = df_.copy(deep=True)

    sno_sno_df = df.loc[(df.biotype1 == 'snoRNA') &
                          (df.biotype2 == 'snoRNA')]

    rev_cols = []
    for col in df.columns:
        if '1' in col:
            rev_cols.append(col.replace('1', '2'))
        elif '2' in col:
            rev_cols.append(col.replace('2', '1'))
        else:
            rev_cols.append(col)

    rev_df = sno_sno_df.copy(deep=True)
    rev_df.columns = rev_cols

    final_df = pd.concat([df, rev_df], ignore_index=True, sort=False)

    return final_df


def write_df(df, out_file):

    df.to_csv(out_file, sep='\t', index=None)


def main():

    # Load the data
    df = load_and_filter_df(in_file)

    # Load the references files
    snodb_df, ref_df = load_refs(sno_df_file, gene_bed_biotype_file)

    # Find the sno-host interactions
    df_with_host = flag_host_interactions(df, snodb_df, ref_df)

    # find the overlap between 3 methods
    method_ovelap(df_with_host)

    # double and inverse the sno_sno interactions
    double_df = double_sno_sno(df_with_host)

    write_df(double_df, out_file)


if __name__ == '__main__':
    main()
