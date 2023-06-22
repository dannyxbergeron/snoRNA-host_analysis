import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']
# sns.set_theme()

from pybedtools import BedTool as bt

snodb_file = snakemake.input.snodb
gene_bed_file = snakemake.input.gene_bed
sno_host_file = snakemake.input.sno_host
sno_host_loc_file = snakemake.input.sno_host_loc
full_network_file = snakemake.input.full_network_interactions
network_host_file = snakemake.input.net_sno_host

out_file = snakemake.output.svg


def load_df(file):
    return pd.read_csv(file, sep='\t')

def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=False, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id1', 'gene_name1', 'strand1', 'biotype1',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'biotype2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df

def process(ref_df_, snodb_df, sno_host_df, sno_host_loc_df, data_df, network_df):

    id_name_dict = dict(zip(ref_df_.gene_id, ref_df_.gene_name))

    sno_df = ref_df_.copy(deep=True)
    sno_df = sno_df.loc[sno_df.gene_biotype == 'snoRNA']
    TOTAL_SNORNA = set(sno_df.gene_id)

    other_ref_df = ref_df_.copy(deep=True)
    other_ref_df = other_ref_df.loc[~(other_ref_df.gene_biotype.isin([
        'snoRNA', 'miRNA', 'snRNA', 'rRNA', 'tRNA', 'scaRNA'
    ]))]

    intersect_df = bedtools(sno_df, other_ref_df)
    inter_set = set(intersect_df.gene_id1)

    INTERGENIC_SNORNA = TOTAL_SNORNA - set(intersect_df.gene_id1)
    EMBEDDED = TOTAL_SNORNA - INTERGENIC_SNORNA

    embeded_snoRNA = sno_df.copy(deep=True)

    snodb_df_dict = dict(zip(snodb_df.gene_id_annot2020, snodb_df['host gene id']))
    sno_host_dict = dict(zip(sno_host_df.sno_id, sno_host_df.host_id))

    PROT_CODING_HOST = set([
        x for x in sno_host_df.sno_id
        if not pd.isna(x)
        and x != 'ENSG00000239039'
    ])
    NON_CODING_HOST = EMBEDDED - PROT_CODING_HOST

    FILTERED_P_CODING_HOST = set(sno_host_loc_df.gene_id)
    OUT_P_CODING_HOST = PROT_CODING_HOST - FILTERED_P_CODING_HOST

    NETWORK_SNO = set(network_df.single_id1).intersection(FILTERED_P_CODING_HOST)
    NON_NETWORK_SNO = FILTERED_P_CODING_HOST - NETWORK_SNO
    
    NETWORK_SNO_HOST_INTERACTION = set(data_df.single_id1)
    NETWORK_SNO_NO_HOST_INTERACTION = NETWORK_SNO - NETWORK_SNO_HOST_INTERACTION

    INTRA_NET_SNO = set(data_df.loc[data_df.interaction_type == 'intra'].single_id1)
    OTHER_NET_SNO = NETWORK_SNO_HOST_INTERACTION - INTRA_NET_SNO


    print(f'Total number of snoRNA: {len(TOTAL_SNORNA)}')
    print('--------------------------------------------')
    print(f'Intergenic snoRNA: {len(INTERGENIC_SNORNA)}')
    print(f'Embedded snoRNA: {len(EMBEDDED)}')
    print('--------------------------------------------')
    print(f'snoRNA in non coding genes : {len(NON_CODING_HOST)}')
    print(f'snoRNA in protein coding genes : {len(PROT_CODING_HOST)}')
    print('--------------------------------------------')
    print(f'snoRNA in filtered out protein coding genes : {len(OUT_P_CODING_HOST)}')
    print(f'snoRNA in filtered protein coding genes : {len(FILTERED_P_CODING_HOST)}')
    print('--------------------------------------------')
    print(f'Filtered snoRNA Not in network : {len(NON_NETWORK_SNO)}')
    print(f'Filtered snoRNA in network : {len(NETWORK_SNO)}')
    print('--------------------------------------------')
    print(f'In network and having host interaction : {len(NETWORK_SNO_HOST_INTERACTION)}')
    print(f'In network, no host interaction : {len(NETWORK_SNO_NO_HOST_INTERACTION)}')
    print('--------------------------------------------')
    print(f'snoRNA interacting elsewhere : {len(OTHER_NET_SNO)}')
    print(f'snoRNA interacting in the same intron : {len(INTRA_NET_SNO)}')
    print('--------------------------------------------')

    values = [
        len(TOTAL_SNORNA),
        len(INTERGENIC_SNORNA),
        len(EMBEDDED),
        len(NON_CODING_HOST) + len(OUT_P_CODING_HOST),
        len(PROT_CODING_HOST) - len(OUT_P_CODING_HOST),
        # len(OUT_P_CODING_HOST),
        # len(FILTERED_P_CODING_HOST),
        len(NON_NETWORK_SNO),
        len(NETWORK_SNO),
        len(NETWORK_SNO_NO_HOST_INTERACTION),
        len(NETWORK_SNO_HOST_INTERACTION),
        len(OTHER_NET_SNO),
        len(INTRA_NET_SNO),
    ]

    groups = [
        'Total',
        'Intergenic',
        'Embedded',
        'In non-coding HG or transcripts, exonic or in opposite strand as the HG',
        'In protein coding',
        # 'Filtered out',
        # 'Filtered',
        'Not detected in PARIS/LIGR-seq/SPLASH',
        'Detected in PARIS/LIGR-seq/SPLASH',
        'No host interaction',
        'Host interaction',
        'Interaction elsewhere',
        'Same intron interaction'
    ]

    groups = [f'{g} ({v})' for g,v in zip(groups, values)]
    values = [x/len(TOTAL_SNORNA) for x in values]



    return groups, values


def sankey_graph(groups, values):

    fig, ax = plt.subplots(figsize=(16, 8))

    sankey = Sankey(ax=ax, unit=None, shoulder=.075,)
    sankey.add(
        flows=[values[0], -values[2], -values[1]],
        orientations=[0, 0, -1],
        labels=[groups[0], groups[2], groups[1]],
        facecolor='#808080',
        edgecolor='white',
    )
    sankey.add(
        flows=[values[2], -values[4], -values[3]],
        orientations=[0, 0, -1],
        pathlengths=[0, 0, 0.3],
        labels=['', groups[4], groups[3]],
        prior=0,
        connect=(1, 0),
        facecolor='#333333',
        edgecolor='white',
    )
    # sankey.add(
    #     flows=[values[4], -values[6], -values[5]],
    #     orientations=[0, 0, -1],
    #     pathlengths=[0, 0, 0.3],
    #     labels=['', groups[6], groups[5]],
    #     prior=1,
    #     connect=(1, 0),
    #     facecolor='#fdbf6f',
    #     edgecolor='white',
    # )
    sankey.add(
        flows=[values[4], -values[6], -values[5]],
        orientations=[0, 0, -1],
        pathlengths=[0, 0, 0.3],
        labels=['', groups[6], groups[5]],
        prior=1,
        connect=(1, 0),
        facecolor='#e41a1c',
        edgecolor='white',
    )
    
    sankey.add(
        flows=[values[6], -values[8], -values[7]],
        orientations=[0, 0, -1],
        pathlengths=[0, 0.5, 0.2],
        labels=['', groups[8], groups[7]],
        prior=2,
        connect=(1, 0),
        facecolor='#756bb1',
        edgecolor='white',
    )
    
    sankey.add(
        flows=[values[8], -values[10], -values[9]*1.2],
        orientations=[0, 0, -1],
        pathlengths=[-.3, 0, 0.2],
        labels=['', groups[10], groups[9]],
        prior=3,
        connect=(1, 0),
        facecolor='#1f78b4',
        edgecolor='white',
    )
    sankey.finish()

    plt.xticks([])
    plt.yticks([])

    left, right = plt.xlim()
    plt.xlim(left, right*1.05)

    plt.savefig(out_file, format='svg')
    # plt.show()


def main():

    ref_df = load_df(gene_bed_file)
    snodb_df = load_df(snodb_file)
    sno_host_df = load_df(sno_host_file)
    sno_host_loc_df = load_df(sno_host_loc_file)
    full_network_df = load_df(full_network_file)
    data_df = load_df(network_host_file)

    groups, values = process(ref_df, snodb_df, sno_host_df, sno_host_loc_df, data_df, full_network_df)

    sankey_graph(groups, values)

if __name__ == '__main__':
    main()
