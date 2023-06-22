import sys
import random
from statistics import stdev

import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

import multiprocessing as mp

sno_host_file = snakemake.input.sno_host_loc
bedgraph_file = snakemake.input.bedgraph
data_file = snakemake.input.cons


THRESH = snakemake.params.threshold

# sno_host_file = 'prot_cod_sno_host_loc.tsv'
# bedgraph_file = 'simplified_host_bedgraph.bedgraph'
# data_file = 'sno_host_with_cons.tsv'


def load_df(file):
    return pd.read_csv(file, sep='\t')

def load_bedgraph():
    colnames = ['chr', 'start', 'end', 'score']
    return pd.read_csv(bedgraph_file, sep='\t', names=colnames)

def get_stats(df_):

    df = df_.copy(deep=True)

    dist = []
    side = []
    locs = ['sno_start', 'sno_end', 'start2', 'end2', 'strand1']
    for s1, e1, s2, e2, strand in df[locs].values:
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

    df['side'], df['dist'] = side, dist

    downstream_ratio = list(df.side).count('downstream') / len(df.side)
    median_dist = df.dist.median()
    int_length_mean = df.other_len_cons.mean()
    stdev_length = stdev(df.other_len_cons)

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

def get_cons(data_df, df, bedgraph_df, net_stats):

    CONSERVATION_OFFSET = 2

    d_ratio, med_dist, len_inter, stdev_length = net_stats
    inverse_side ={0: 1, 1: 0}

    sides = {}
    bt_list = []
    cols = ['chr', 'start', 'end', 'intron_start',
            'intron_end', 'strand', 'gene_id', 'gene_name']
    for chr, start, end, intron_start, intron_end, strand, \
                gene_id, gene_name in df[cols].values:

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

    return df


def stats(df, ref_df_cons):

    def get_vals_wrapper(original_list):
        def get_vals(list_):
            pos, neg = 0, 0
            for v in list_:
                if v >= THRESH: pos += 1
                else: neg +=1
            return pos, neg
        vals = get_vals(original_list)
        return vals[0] / sum(vals)

    network_ratio = get_vals_wrapper(df['other_mean_cons'])
    other_ratio = get_vals_wrapper(ref_df_cons['mean_cons'])

    # print((other_ratio > network_ratio), network_ratio, other_ratio)
    return [other_ratio > network_ratio, network_ratio, other_ratio]


def run(args):
    func_args, n = args
    return [
        stats(func_args[0], get_cons(*func_args))
        for _ in range(n)
    ]


def main():

    # load and process data
    df = load_df(data_file)
    df = df.loc[df.interaction_type == 'intra']
    # df.drop_duplicates(subset=['merged_name'], inplace=True)

    # load and process the reference
    ref_df = load_df(sno_host_file)
    ref_df = ref_df.loc[~ref_df.gene_id.isin(df.single_id1)]

    # load the conservation bedgraph
    bedgraph_df = load_bedgraph()

    # get the stats for the data
    net_stats = get_stats(df)

    # ========================= MULTIPROCESSING ==========================
    N = 16
    CORES = 6
    pool = mp.Pool(processes = CORES)

    # Map apply_func function with data and my_function to cores
    results = pool.map(run, [
        ((df, ref_df, bedgraph_df, net_stats), len(it))
        for it in np.array_split(range(N), CORES)
    ])
    # Wait for everything to close and close cores
    pool.close()
    # Results is a list of lists, flatten
    results = [item for sublist in results for item in sublist]
    print('\n'.join([
        '\t'.join([str(val) for val in vals])
        for vals in results
    ]))

    # print('---------------------------')
    # if not results.count(False):
    #     print(f'pvalue: 0')
    # else:
    #     print(f'pvalue: {results.count(True) / len(results)}')



if __name__ == '__main__':
    main()
