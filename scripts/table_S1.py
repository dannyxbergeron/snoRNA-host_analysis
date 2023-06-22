import numpy as np
import pandas as pd

in_file = snakemake.input.ext_ratio
out_file = snakemake.output.tableS1


def beautify(df_):

    df = df_.copy(deep=True)
    df = df[[
        'DG', 'chr1',
        'start1', 'end1',
        'start2', 'end2',
        'strand1',
        'single_id1', 'single_id2',
        'name1', 'name2',
        'intron_start', 'intron_end',
        'sno_tpm', 'host_tpm',
        'support', 'E', 'other_mean_cons', 'splice_dist', 'ext_ratio'
    ]]
    df.rename(columns={
        'chr1': 'chr',
        'strand1': 'strand',
        'single_id1': 'snoRNA_id', 'single_id2': 'host_id',
        'name1': 'snoRNA_name', 'name2': 'host_name',
        'sno_tpm': 'snoRNA_mean_tpm', 'host_tpm': 'host_mean_tpm',
        'other_mean_cons': 'target_mean_cons',
        'ext_ratio': 'extension_ratio'

    }, inplace=True)
    df.sort_values(by=['target_mean_cons', 'support'],
                   ascending=[False, False],
                   inplace=True)
    return df

def main():

    df = pd.read_csv(in_file, sep='\t')
    df = beautify(df)

    df.to_csv(out_file, index=False)


if __name__ == '__main__':
    main()
