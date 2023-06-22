import random
from collections import Counter
from statistics import stdev

from collections import defaultdict
from statistics import mean, median

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu as mwu

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']

from pybedtools import BedTool as bt

data_file = snakemake.input.alt_splice
sno_host_file = snakemake.input.sno_host_loc
bg_files = snakemake.input.bg

bed_viz = snakemake.output.bed_viz
out_file = snakemake.output.ext_ratio


def load_df(file, header=True):
    if header:
        df = pd.read_csv(file, sep='\t')
    else:
        df = pd.read_csv(file, sep='\t',
                         names=['chr', 'start', 'end', 'value'])
    return df

def load_beds(files):

    dfs = {}
    for file in files:
        name = file.split('/')[-1].split('.')[0].replace('intron_', '')
        dfs[name] = load_df(file, header=False)

    return dfs


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

    print('--------------------- STATS ---------------------')
    counter = Counter(list(df.side))
    downstream_ratio = counter['downstream'] / (counter['downstream'] + counter['upstream'])
    median_dist = df.dist.median()
    int_length_mean = df.other_len_cons.mean()
    stdev_length = stdev(df.other_len_cons)

    print(f'downstream ratio: {downstream_ratio:.2f}')
    print(f'median distance of the interaction {median_dist}')
    print(f'average lenth of the interacting region: {int_length_mean:.2f}, stdev: {stdev_length:.2f}')
    print('-------------------------------------------------')

    return downstream_ratio, median_dist, int_length_mean, stdev_length

def create_fake(data_df, df):

    d_ratio, med_dist, len_inter, stdev_length = get_stats(data_df)
    inverse_side ={0: 1, 1: 0}

    sides = {}
    fake_list = []
    for i in df.index:
        chr = df.at[i, 'chr']
        start = df.at[i, 'start']
        end = df.at[i, 'end']
        intron_start = df.at[i, 'intron_start']
        intron_end = df.at[i, 'intron_end']
        strand = df.at[i, 'strand']
        gene_id = df.at[i, 'gene_id']
        gene_name = df.at[i, 'gene_name']
        host_name = df.at[i, 'host_name']
        host_id = df.at[i, 'host_id']

        len_offset = random.randint(0, int(stdev_length))
        if random.choice([0, 1]):
            int_len = len_inter + len_offset
        else:
            int_len = len_inter - len_offset
        side = random.choice([0, 1])
        if strand == '-': side = inverse_side[side]

        if side and intron_end - end > 8:
            int_len = int_len if intron_end - end > int_len else intron_end - end
            new_start, new_end = end, end + int_len
        elif side == 0 and start - intron_start > 8:
            int_len = int_len if start - intron_start > int_len else start - intron_start
            new_start, new_end = start - int_len, start
        elif side and start - intron_start > 8:
            int_len = int_len if start - intron_start > int_len else start - intron_start
            side = inverse_side[side]
            new_start, new_end = start - int_len, start
        elif side == 0 and intron_end - end > 8:
            int_len = int_len if intron_end - end > int_len else intron_end - end
            side = inverse_side[side]
            new_start, new_end = end, end + int_len
        else:
            print('ERROR !!!!!!!')

        if (side and strand == '+') or (not side and strand == '-'): s = 'upstream'
        else: s = 'downstream'
        sides[gene_id] = s

        DG = f'FAKE_{gene_id}-{host_id}'
        fake_info = [
            chr, '-'.join([str(int(new_start)), str(int(new_end))]), intron_start,
            intron_end, start, end, strand, DG, gene_name, host_name
        ]
        fake_list.append(fake_info)

    colnames = [
        'chr1', 'int_portion_start2', 'intron_start', 'intron_end',
        'sno_start', 'sno_end', 'strand1', 'DG', 'name1', 'name2'
    ]
    fake_df = pd.DataFrame(fake_list, columns=colnames)

    return fake_df


def get_regions(df_):

    THRESH = 5
    OFFSET = 2

    df = df_.copy(deep=True)
    df = df.loc[~pd.isnull(df.int_portion_start2)]
    df_val = df[[
        'chr1', 'int_portion_start2', 'intron_start', 'intron_end',
        'sno_start', 'sno_end', 'strand1', 'DG', 'name1', 'name2'
    ]].values

    master_list = []
    for chr, portion2, int_start, int_end, sno_start, sno_end, strand, DG, name1, name2 in df_val:
        p_start, p_end = [int(x) for x in portion2.split('-')]
        sno_start -= (OFFSET + 1)
        sno_end += OFFSET
        if p_start < sno_start:
            ext_start, ext_end = p_start, sno_start
            snoExt_start, snoExt_end = p_start, sno_end
        else:
            ext_start, ext_end = sno_end, p_end
            snoExt_start, snoExt_end = sno_start, p_end

        ext = [chr, ext_start, ext_end, DG, 'ext', strand]
        master_list.append(ext)

        left = [chr, int_start, snoExt_start, DG, 'left', strand]
        right = [chr, snoExt_end, int_end - 1, DG, 'right', strand]

        if snoExt_start - int_start > THRESH:
            master_list.append(left)
        if int_end - snoExt_end > THRESH:
            master_list.append(right)

    master_df = pd.DataFrame(master_list, columns=['chr', 'start', 'end', 'DG', 'side', 'strand'])
    master_df['chr'] = 'chr' + master_df['chr']
    master_df['start'] = master_df.start.map(int)
    master_df['end'] = master_df.end.map(int)

    viz_df = master_df.copy(deep=True)
    viz_df = viz_df.loc[~(viz_df.DG.str.startswith('FAKE'))]
    viz_df.to_csv(bed_viz, sep='\t', index=False, header=False)

    return master_df


def bedtools(df1, df2):

    df1 = df1.sort_values(['chr', 'start', 'end'])
    df2 = df2.sort_values(['chr', 'start', 'end'])

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=False, sorted=True)
    new_cols = ['chr', 'start', 'end', 'DG', 'side', 'strand',
                'chr2', 'start2', 'end2', 'score', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})

    intersect_df['sum_score'] = intersect_df['score'] * intersect_df['overlap']

    return intersect_df


def get_values(df_, beds):

    df = df_.copy(deep=True)
    # Just for testing with SNORD2...
    # df = df.loc[df.DG == '273061_L0|273063_L0|488949_L1|488953_L1|692461_P0']

    means_ratio = defaultdict(list)
    for name, bed in beds.items():
        print(name)

        intersect_df = bedtools(df, bed)

        ext = intersect_df.loc[intersect_df.side == 'ext'].copy(deep=True)
        rest = intersect_df.loc[~(intersect_df.side == 'ext')].copy(deep=True)

        for DG in set(df.DG):
            tmp_ext = ext.loc[ext.DG == DG]
            tmp_rest = rest.loc[rest.DG == DG]

            if len(tmp_ext) == 0:
                ext_val = 0.5
            else:
                ext_values = tmp_ext.values[0]
                ext_sum = ext_values[2] - ext_values[1]
                ext_val = tmp_ext.sum_score.sum() / ext_sum

            if len(tmp_rest) == 0:
                rest_val = 0.5
            else:
                rest_gb = tmp_rest.groupby('side').median().reset_index()
                rest_sum = sum([
                    rest_gb.at[x, 'end'] - rest_gb.at[x, 'start']
                    for x in rest_gb.index
                ])
                rest_val = tmp_rest.sum_score.sum() / rest_sum

            ratio = ext_val / rest_val
            means_ratio[DG].append(ratio)


    final_dict = {}
    for DG in set(df.DG):
        final_dict[DG] = median(means_ratio[DG]) # CHANGED, no big difference

    return final_dict

def graph(combined_df, ratio_dict):

    MAX_VAL = 20

    df = combined_df.copy(deep=True)
    df['ext_ratio'] = combined_df.DG.map(ratio_dict)

    in_net_df = df.loc[~(df.DG.str.startswith('FAKE'))]
    others_df = df.loc[df.DG.str.startswith('FAKE')]

    print('\n=========== STATS - mann-whitney u test ==========')
    net_stats, net_pval = mwu(in_net_df['ext_ratio'],
                              others_df['ext_ratio'],
                              alternative='two-sided')
    print(f'For network host interacting vs others p_value: {net_pval}')
    print('===================================================\n')

    in_net = np.where(in_net_df['ext_ratio'] <= MAX_VAL,
                      in_net_df['ext_ratio'], MAX_VAL)
    others = np.where(others_df['ext_ratio'] <= MAX_VAL,
                      others_df['ext_ratio'], MAX_VAL)

    groups = [
        'snoRNA in the network',
        'others'
    ]

    fig, ax = plt.subplots()
    fig.canvas.draw()

    sns.kdeplot(data=in_net, shade=True, linewidth=1, alpha=.3,
                label=groups[0], ax=ax, bw_adjust=.7,
                color='#377eb8')
    sns.kdeplot(others, shade=True, linewidth=1, alpha=.3,
                label=groups[1], ax=ax, bw_adjust=.7,
                color='#e41a1c')

    print(in_net_df['ext_ratio'].mean(), in_net_df['ext_ratio'].median())
    print(others_df['ext_ratio'].mean(), others_df['ext_ratio'].median())

    print(in_net_df[['ext_ratio']].describe())
    print(others_df[['ext_ratio']].describe())

    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-2] = str(tick_labels[-2]) + '+'
    ax.set_xticklabels(tick_labels)

    left, right = plt.xlim()
    plt.xlim(left, 29)

    plt.title('Distribution of ratio of reads in extension compared to the intron of the snoRNA')
    plt.xlabel('ratio of reads in extension (extension/rest of intron)')
    plt.legend()
    plt.show()


def main():

    df = load_df(data_file)
    ref_df = load_df(sno_host_file)
    beds = load_beds(bg_files)


    # TEST ---------------------------------
    # print(f'original len ref: {len(ref_df)}')
    # print(f'original len data: {len(df)}')
    # ref_df = ref_df.loc[ref_df.sno_tpm > 10]
    # df = df.loc[df.single_id1.isin(ref_df.gene_id)]
    # print(f'filtered len ref: {len(ref_df)}')
    # print(f'filtered len data: {len(df)}')
    # TEST ---------------------------------

    # Added to create false regions in other snoRNAs
    ref_df_ = ref_df.loc[~(ref_df.gene_id.isin(df.single_id1))]

    fake_df = create_fake(df, ref_df_)

    combined_df = pd.concat([df, fake_df], axis=0, sort=False)

    df_regions = get_regions(combined_df)
    # df_regions = get_regions(df)

    ratio_dict = get_values(df_regions, beds)

    graph(combined_df, ratio_dict)

    df['ext_ratio'] = df['DG'].map(ratio_dict)
    df.sort_values('ext_ratio', ascending=False, inplace=True)

    df['sno_tpm'] = df.single_id1.map(dict(zip(ref_df.gene_id, ref_df.sno_tpm)))
    df['host_tpm'] = df.single_id2.map(dict(zip(ref_df.host_id, ref_df.target_tpm)))

    df.to_csv(out_file, sep='\t', index=False)



if __name__ == '__main__':
    main()
