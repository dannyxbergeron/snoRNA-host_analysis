import random
from collections import Counter
from statistics import stdev

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu

import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']
# sns.set_theme()

from pybedtools import BedTool as bt

sno_host_file = snakemake.input.sno_host_loc
bedgraph_file = snakemake.input.bedgraph
data_file = snakemake.input.cons

# svg = snakemake.output.svg

THRESH = snakemake.params.threshold


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df

def load_bedgraph():
    df = pd.read_csv(bedgraph_file, sep='\t', names=['chr', 'start', 'end', 'score'])
    return df

def get_stats(df_):

    df = df_.copy(deep=True)

    side = []
    dist = []
    for s1, e1, s2, e2, strand in df[['sno_start', 'sno_end', 'start2', 'end2', 'strand1']].values:
        if strand == '-':
            max_val = max(e1, e2)
            s1, e1, s2, e2 = [(v - max_val) * -1 for v in [e1, s1, e2, s2]]

        if s2 < s1:
            side.append('upstream')
            d = s1 - e2 if s1 - e2 > 0 else 0
            dist.append(d)
        else:
            side.append('downstream')
            d = s2 - e1 if s2 - e1 > 0 else 0
            dist.append(d)

    df['side'] = side
    df['dist'] = dist

    df['min'] = df[['start2', 'end2', 'sno_start', 'sno_end']].min(axis=1)
    df['start2'] = df['start2'] - df['min']
    df['end2'] = df['end2'] - df['min']
    df['sno_start'] = df['sno_start'] - df['min']
    df['sno_end'] = df['sno_end'] - df['min']

    print(df[['DG', 'start2', 'end2', 'sno_start', 'sno_end', 'strand1', 'side', 'dist']].sort_values('dist', ascending=False))


    # print(df)

    print('--------------------- STATS ---------------------')

    counter = Counter(list(df.side))
    downstream_ratio = counter['downstream'] / (counter['downstream'] + counter['upstream'])
    median_dist = df.dist.median()
    int_length_mean = df.other_len_cons.mean()
    stdev_length = stdev(df.other_len_cons)

    print(f'downstream ratio: {downstream_ratio:.2f}')
    print(f'median distance of the interaction {median_dist}')
    print(f'average lenth of the interacting region: {int_length_mean:.2f}, stdev: {stdev_length:.2f}')
    print(f'------> Average cons {df.other_mean_cons.mean()}')

    return downstream_ratio, median_dist, int_length_mean, stdev_length


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, sorted=True)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id', 'gene_name', 'strand',
                'chr2', 'start2', 'end2', 'score', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df

def get_cons(data_df, df, bedgraph_df):

    CONSERVATION_OFFSET = 2

    d_ratio, med_dist, len_inter, stdev_length = get_stats(data_df)
    inverse_side ={0: 1, 1: 0}

    sides = {}
    bt_list = []
    for i in df.index:
        chr = df.at[i, 'chr']
        start = df.at[i, 'start']
        end = df.at[i, 'end']
        intron_start = df.at[i, 'intron_start']
        intron_end = df.at[i, 'intron_end']
        strand = df.at[i, 'strand']
        gene_id = df.at[i, 'gene_id']
        gene_name = df.at[i, 'gene_name']

        len_offset = random.randint(0, int(stdev_length))
        if random.choice([0, 1]):
            int_len = len_inter + len_offset
        else:
            int_len = len_inter - len_offset
        side = random.choice([0, 1])
        if strand == '-': side = inverse_side[side]

        if side and intron_end - end > 8:
            int_len = int_len if intron_end - end > int_len else intron_end - end
            new_start, new_end = end + CONSERVATION_OFFSET, end + int_len
        elif side == 0 and start - intron_start > 8:
            int_len = int_len if start - intron_start > int_len else start - intron_start
            new_start, new_end = start - int_len, start - CONSERVATION_OFFSET
        elif side and start - intron_start > 8:
            int_len = int_len if start - intron_start > int_len else start - intron_start
            side = inverse_side[side]
            new_start, new_end = start - int_len, start - CONSERVATION_OFFSET
        elif side == 0 and intron_end - end > 8:
            int_len = int_len if intron_end - end > int_len else intron_end - end
            side = inverse_side[side]
            new_start, new_end = end + CONSERVATION_OFFSET, end + int_len
        else:
            print('ERROR !!!!!!!')

        if (side and strand == '+') or (not side and strand == '-'): s = 'upstream'
        else: s = 'downstream'
        sides[gene_id] = s

        bt_list.append(['chr'+chr, int(new_start), int(new_end), gene_id, gene_name, strand])

    colnames = ['chr', 'start', 'end', 'gene_id', 'gene_name', 'strand']
    bt_df = pd.DataFrame(bt_list, columns=colnames)
    bt_df.sort_values(['chr', 'start', 'end'], inplace=True)
    print(df)
    print(bt_df)

    intersect_df = bedtools(bt_df, bedgraph_df)

    intersect_df['product'] = intersect_df.overlap * intersect_df.score
    wanted_cols = ['gene_id', 'start1', 'end1', 'product']

    tmp_intersect = intersect_df[wanted_cols].groupby('gene_id').agg({'start1': 'min',
                                                                      'end1':'max',
                                                                      'product':'sum'}).reset_index()
    tmp_intersect['length'] = tmp_intersect['end1'] - tmp_intersect['start1']
    tmp_intersect['mean'] = tmp_intersect['product'] / tmp_intersect['length']

    df['length_cons'] = df.gene_id.map(dict(zip(tmp_intersect.gene_id, tmp_intersect['length'])))
    df['mean_cons'] = df.gene_id.map(dict(zip(tmp_intersect.gene_id, tmp_intersect['mean'])))
    df['side'] = df.gene_id.map(sides)

    df.fillna(0, inplace=True)
    df.sort_values(['mean_cons'], inplace=True)
    # print(df.describe())

    # print(df)

    print(f'------> Average cons for the other: {df.mean_cons.mean()}')

    return df


def graph(df, ref_df_cons):

    print('\n=========== STATS - mann-whitney u test ==========')
    net_stats, net_pval = mwu(df['other_mean_cons'],
                              ref_df_cons['mean_cons'],
                              alternative='two-sided')
    print(f'For network host interacting vs others p_value: {net_pval}')
    print('===================================================\n')


    def get_vals(list_):
        pos, neg = 0, 0
        for v in list_:
            if v >= THRESH: pos += 1
            else: neg +=1
        return pos, neg

    network_values = get_vals(df['other_mean_cons'])
    other_values = get_vals(ref_df_cons['mean_cons'])
    print(f'Network vals: {network_values}, other vals: {other_values}')


    groups = [
        'snoRNA in the network',
        'others'
    ]

    fig, ax = plt.subplots()
    fig.canvas.draw()

    sns.kdeplot(data=df['other_mean_cons'], shade=True, linewidth=1, alpha=.3,
                label=groups[0], ax=ax, bw_adjust=1,
                color='#377eb8')
    sns.kdeplot(ref_df_cons['mean_cons'], shade=True, linewidth=1, alpha=.3,
                label=groups[1], ax=ax, bw_adjust=1,
                color='#e41a1c')

    plt.title('Distribution of conservation around snoRNAs')
    plt.xlabel('Average conservation of the regions')
    plt.legend()
    # plt.savefig(svg, format='svg')
    plt.show()


def main():

    df = load_df(data_file)
    df = df.loc[df.interaction_type == 'intra']
    # df.drop_duplicates(subset=['merged_name'], inplace=True)

    ref_df = load_df(sno_host_file)
    ref_df = ref_df.loc[~ref_df.gene_id.isin(df.single_id1)]

    print(ref_df.columns)
    bedgraph_df = load_bedgraph()

    ref_df_cons = get_cons(df, ref_df, bedgraph_df)

    graph(df, ref_df_cons)



if __name__ == '__main__':
    main()
