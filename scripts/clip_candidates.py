from collections import defaultdict
import ast

import numpy as np
import pandas as pd

file = snakemake.input.full_merged

out_file = snakemake.output.tmp_candidates
out_bed = snakemake.output.bed_file

def load_df():
    df = pd.read_csv(file, sep='\t', dtype={'prot_info': str})
    df = df.dropna(subset=['prot_info'])
    df = df.loc[df.biotype2 == 'protein_coding']

    rbp_count = defaultdict(int)
    for i in df.index:
        prot_info = df.at[i, 'prot_info']
        rbp = []
        for prot_list in prot_info.replace(('['), '').replace(']', '').replace('),', ')),').split('), '):
            tup = ast.literal_eval(prot_list)
            rbp.append(tup[0])
            rbp_count[tup[0]] += 1
        df.at[i, 'prot_info'] = ','.join(rbp)

    sorted_rbp_count = {
        k: v
        for k, v in sorted(rbp_count.items(), key=lambda item: item[1], reverse=True)
    }
    return df, sorted_rbp_count


def main():

    df_, prot_dict = load_df()

    # To create the bed to vizualize the interactions
    bed_df = df_.copy(deep=True)
    bed_df = bed_df[['chr2', 'start2', 'end2', 'name1', 'E', 'strand2']]
    bed_df.sort_values(['chr2', 'start2', 'end2'], inplace=True)
    bed_df.to_csv(out_bed, sep='\t', index=False, header=False)


    df = df_.copy(deep=True)
    df = df[['DG', 'chr2', 'start2', 'end2', 'name1',
             'name2', 'host_interaction', 'E', 'prot_info']]

    df.columns = ['DG', 'chr', 'start', 'end', 'snoRNA',
                    'gene_name', 'host_interaction', 'E', 'prot_info']
    df.sort_values(['host_interaction', 'E'], inplace=True)

    df.to_csv(out_file, index=False)




if __name__ == '__main__':
    main()
