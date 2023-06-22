import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

file = snakemake.input.merged_file
gab_host_file = '/home/danx/Documents/projects/projectarticlefamille/dan_reanalysis/data/ref/final_host_gene_101.csv'

def get_num_sno_per_host(host_df):

    print(host_df.columns)
    gb_df = host_df[[
        'gene_id', 'sno_id'
    ]].groupby('gene_id').count()

    gb_df = gb_df.reset_index()
    # print(gb_df.loc[gb_df.gene_id == 'ENSG00000156976'])

    return dict(zip(gb_df.gene_id, gb_df.sno_id))

def main():

    df = pd.read_csv(file, sep='\t')
    host_df = pd.read_csv(gab_host_file)
    host_df = host_df.loc[host_df.gene_biotype == 'protein_coding']

    df = df[[
        'single_id1', 'name1',
        'single_id2', 'name2',
        'host_interaction'
    ]]

    host_interacting_df = df.loc[df.host_interaction]

    sno_per_host_dict = get_num_sno_per_host(host_df)

    # ---------------------------------------------------------
    all_sno_ids = set(df.single_id1)
    host_int_ids = set(host_interacting_df.single_id1)
    all_sno_ids = all_sno_ids - host_int_ids # Check all the sno not interacting with their host !!! CHANGED !

    intronic_all_sno_ids = all_sno_ids.intersection(set(host_df.sno_id))
    intronic_host_int_ids_with_host = host_int_ids.intersection(set(host_df.sno_id))

    print(f'Total number of snoRNAs: {len(all_sno_ids)}')
    print(f'Number of snoRNA interacting with their host: {len(host_int_ids)}')
    print('----------------------------------------\n')

    print(f'Total number of intronic snoRNAs: {len(intronic_all_sno_ids)}')
    print(f'Number of intronic snoRNA interacting with their host: {len(intronic_host_int_ids_with_host)}')
    print('----------------------------------------\n')


    sno_host_dict = dict(zip(host_df.sno_id, host_df.gene_id))
    all_num_sno_per_gene = [
        sno_per_host_dict[sno_host_dict[sno_id]]
        for sno_id, host_id in sno_host_dict.items()
        if sno_id in intronic_all_sno_ids
    ]

    print(all_num_sno_per_gene)
    print(np.mean(all_num_sno_per_gene))


    sno_host_num_sno_per_gene = [
        sno_per_host_dict[sno_host_dict[sno_id]]
        for sno_id, host_id in sno_host_dict.items()
        if sno_id in intronic_host_int_ids_with_host
    ]

    print(sno_host_num_sno_per_gene)
    print(np.mean(sno_host_num_sno_per_gene))


    sns.kdeplot(data=all_num_sno_per_gene, shade=True, linewidth=1, alpha=.3,
                    label='All intronic network snoRNAs not interacting', color='green')
    sns.kdeplot(data=sno_host_num_sno_per_gene, shade=True, linewidth=1, alpha=.3,
                    label='Intronic network snoRNAs interacting with their host', color='red')

    plt.xlabel('Number of snoRNA in the host gene')
    plt.ylabel('Density of snoRNA')

    plt.show()






if __name__ == '__main__':
    main()
