import numpy as np
import pandas as pd
import scipy.stats as stats

import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']
# sns.set_theme()

from pybedtools import BedTool as bt

parsed_file = snakemake.input.parsed
sno_host_loc_file = snakemake.input.sno_host_loc
sno_host_file = snakemake.input.cons

out_file = snakemake.output.alt_splice
# graph1 = snakemake.output.graph1
graph2 = snakemake.output.graph2


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    print(file)
    print(df.columns)
    print('-----------------------------')
    return df


def get_sno_and_host(gtf_df, sno_host_loc_df):

    sno_df = gtf_df.loc[(gtf_df.gene_biotype == 'snoRNA') & (gtf_df.feature == 'gene')]
    prot_cod_df = gtf_df.loc[gtf_df.gene_biotype == 'protein_coding']

    snoRNA_ids = set(sno_host_loc_df.gene_id)
    protein_coding_set = set(prot_cod_df.gene_id)

    sno_host_dict = dict(zip(sno_host_loc_df.gene_id, sno_host_loc_df.host_id))

    colnames = ['seqname', 'start', 'end', 'gene_id', 'gene_name', 'strand']
    sno_df = sno_df.loc[sno_df.gene_id.isin(sno_host_dict.keys())][colnames]
    prot_cod_df = prot_cod_df.loc[prot_cod_df.gene_id.isin(sno_host_dict.values())]

    return sno_host_dict, prot_cod_df, sno_df


def create_introns(df):
    master_list = []
    df_values = df.values
    for i, row in enumerate(df_values):
        chr, start, end, exon_number, exon_id, strand = row
        master_list.append(row)
        if i != len(df) - 1:
            if exon_number < df_values[i+1][3]:
                intron_id = f'intron_{exon_id}'
                if strand == '+':
                    int_start = end
                    int_end = df_values[i + 1][1]
                else:
                    int_start = df_values[i + 1][2]
                    int_end = start
                intron_row = [chr, int_start, int_end,
                              exon_number, intron_id, strand]
                master_list.append(intron_row)
    return pd.DataFrame(master_list, columns=df.columns)


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['seqname1', 'start1', 'end1', 'gene_id', 'gene_name', 'strand1',
                'seqname2', 'start2', 'end2', 'exon_number', 'exon_id', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def get_modulation(int_start, int_end, full_intron_exon_df, bt1_):

    bt1 = bt1_.copy(deep=True)
    sno_start = bt1['start'].values[0]
    sno_end = bt1['end'].values[0]

    bt1['start'] = int_start + 1
    bt1['end'] = int_end - 1

    intersect_df = bedtools(bt1, full_intron_exon_df)

    intersect_df = intersect_df.loc[intersect_df.exon_id.str.contains('intron')]
    min_distance = 1000000
    for start, end, gene_name in intersect_df[['start2', 'end2', 'gene_name']].values:
        if (start == int_start and end != int_end):
            canonical_dist = int_end - sno_end
            alt_dist = end - sno_end
            if alt_dist < 0:
                alt_dist = 0 if end - sno_start > 0 else sno_start - end
            min_distance = min(min_distance, canonical_dist, alt_dist)

            # if gene_name == 'SNORD133':
            #     print('---------------------------------')
            #     print('start == int_start', min_distance)
            #     print(start, end)
            #     print('=================================')

        elif (end == int_end and start != int_start):
            canonical_dist = sno_start - int_start
            alt_dist = sno_start - start
            if alt_dist < 0:
                alt_dist = 0 if sno_end - start > 0 else start - sno_end
            min_distance = min(min_distance, canonical_dist, alt_dist)

            # if gene_name == 'SNORD133':
            #     print('---------------------------------')
            #     print('end == int_end', min_distance)
            #     print(start, end)
            #     print('=================================')
    if min_distance != 1000000:
        return True, min_distance
    return False, np.nan


def get_sno_intron(sno_host_dict, prot_cod_df, sno_df_, sno_host_loc_df):

    sno_df = sno_df_.copy(deep=True)
    ref_df = prot_cod_df.copy(deep=True)

    intron_start = dict(zip(sno_host_loc_df.gene_id, sno_host_loc_df.intron_start))
    intron_end = dict(zip(sno_host_loc_df.gene_id, sno_host_loc_df.intron_end))
    snoId_transId_dict = dict(zip(sno_host_loc_df.gene_id, sno_host_loc_df.host_transcript_id))
    snoId_transName_dict = dict(zip(sno_host_loc_df.gene_id, sno_host_loc_df.host_transcript_name))

    splicing_modulation = []
    distances = []
    for idx in sno_df.index:
        sno_id = sno_df.at[idx, 'gene_id']
        sno_name = sno_df.at[idx, 'gene_name']
        host_id = sno_host_dict[sno_id]
        host_df = ref_df.loc[ref_df.gene_id == host_id]
        sno_data_start = sno_df.at[idx, 'start']
        sno_data_end = sno_df.at[idx, 'end']

        int_start = intron_start[sno_id]
        int_end = intron_end[sno_id]

        bt1 = sno_df.loc[sno_df.index == idx]

        full_exon_df = host_df.loc[host_df.feature == 'exon']
        full_exon_df = full_exon_df[['seqname', 'start', 'end', 'exon_number', 'exon_id', 'strand']]
        full_intron_exon_df = create_introns(full_exon_df)

        s_module, distance = get_modulation(int_start, int_end, full_intron_exon_df, bt1)
        distances.append(distance)
        splicing_modulation.append(s_module)


    filt_sno_df = sno_df.copy(deep=True)
    filt_sno_df['host'] = filt_sno_df.gene_id.map(sno_host_dict)
    filt_sno_df['intron_start'] = filt_sno_df.gene_id.map(intron_start)
    filt_sno_df['intron_end'] = filt_sno_df.gene_id.map(intron_end)
    filt_sno_df['splicing_hypothesis'] = splicing_modulation
    filt_sno_df['distance'] = distances

    print('==================================')
    print('len original: {}, len filtered: {}'.format(len(sno_df), len(filt_sno_df)))
    original = sno_df.copy(deep=True)
    original['host'] = original.gene_id.map(sno_host_dict)
    print('Original nb of hosts:{}, filtered nb: {}'.format(len(set(original.host)),
                                                            len(set(filt_sno_df.host))))
    print('==================================')

    return filt_sno_df


def get_stats(sno_df, sno_host_df):

    def get_ratio(df, tresh):
        tmp = df.copy(deep=True)
        tmp_pos = tmp.loc[(tmp.distance <= tresh) & (tmp.splicing_hypothesis)]
        all_true = list(tmp_pos["splicing_hypothesis"]).count(True)
        all = len(tmp)
        ratio = all_true / all
        return all, all_true, ratio

    def fisher(net_true, nb_net, other_true, nb_other):
        net_neg = nb_net - net_true
        other_neg = nb_other - other_true

        obs = np.array([[net_true, net_neg], [other_true, other_neg]])
        return stats.fisher_exact(obs)

    sno_df_in_data = sno_df.loc[sno_df.gene_id.isin(sno_host_df.single_id1)]
    sno_df_not_in_data = sno_df.loc[~sno_df.gene_id.isin(sno_host_df.single_id1)]

    master_list = []
    threshold = [75, 100, 150, 300, 500, 1000, 10000]
    print(sno_df)
    for thres in threshold:

        nb_all, all_true, all_ratio = get_ratio(sno_df, thres)
        nb_net, net_true, net_ratio = get_ratio(sno_df_in_data, thres)
        nb_other, other_true, other_ratio = get_ratio(sno_df_not_in_data, thres)

        print(f'------------------ {thres} -------------------')
        print(f'All snoRNA ratio: {all_ratio:.2f} ({all_true}/{nb_all})')
        print(f'snoRNA inteacting with their host ratio: {net_ratio:.2f} ({net_true}/{nb_net})')
        print(f'Other snoRNA ratio: {other_ratio:.2f} ({other_true}/{nb_other})\n')
        print('------> Fisher exact test:', end=' ')
        odds, p_val = fisher(net_true, nb_net, other_true, nb_other)
        print(f'Odds: {odds:.2f}, pValue: {p_val:.5f}')

        master_list.append([net_true, nb_net, other_true, nb_other, thres])

    return pd.DataFrame(master_list,
                        columns=['net_true', 'nb_net', 'other_true', 'nb_other', 'thres'])


def graph(df):

    # Data
    r = list(range(len(df)))
    bar_width = 0.30

    # fig, axes = plt.subplots(1, 2, figsize=(12, 8))
    fig, ax = plt.subplots(figsize=(12, 8))
    barWidths = [bar_width, -bar_width]
    names = [str(x) for x in df.thres]
    colnames = [('net_true', 'nb_net'), ('other_true', 'nb_other')]
    colors = ['#377eb8', '#e41a1c']
    titles = ['snoRNA interacting with their host intron\nwith possible alternative splicing',
              'Other snoRNA with a possible alternative splicing']
    labels = ['snoRNA-host interacting', 'snoRNA-host not interacting']
    for (pos, nb), barWidth, color, label, title in zip(colnames, barWidths, colors, labels, titles):

        pos = df[pos] / df[nb] * 100
        rest = 100 - pos

        ax.bar(r, pos, color=color, edgecolor='white', align='edge', width=barWidth, label=label)
        # ax.bar(r, rest, bottom=pos, color='grey',
        #        edgecolor='white', width=barWidth, align='edge')
        ax.set_title(title)

        # Custom x axis
        ax.set_xticks(r)
        ax.set_xticklabels(names)

        ax.set_xlabel("Maximum distance of the splicing event from the snoRNA (pb)")
        ax.set_ylabel("% of snoRNA")

        ax.legend()

    # Show graphic
    plt.legend()
    # plt.savefig(graph1, format='svg')
    # plt.show()
    plt.close()


def stats_km(net_, other_):

    print(len(net_), len(other_))
    net = sorted(list([x for x in net_ if not pd.isnull(x)]))
    other = sorted(list([x for x in other_ if not pd.isnull(x)]))
    # print(len(net), len(other))

    o, p = stats.mannwhitneyu(net, other)
    print('mannwhitneyu: ', o, p)

    maximum = max(max(net), max(other))
    net_idx = 0
    other_idx = 0
    dists = []
    x = []
    for i in range(2000):
        if net_idx < len(net) -1:
            while net[net_idx] == i:
                net_idx += 1
        if other_idx < len(other) - 1:
            while other[other_idx] == i:
                other_idx += 1

        if not dists or (net_idx, other_idx) != dists[-1]:
            dists.append((net_idx, other_idx))
            x.append(i)

    dists = np.array(dists)
    net_dist = dists[:, 0] / len(net_) * 100
    other_dist = dists[:, 1] / len(other_) * 100

    # print(net_dist)
    # print(other_dist)

    # ==================== STATS ==========================
    stat, pvalue = stats.kstest(net_dist, other_dist, alternative='two-sided')
    print(f'kolmogorov-smirnov test, stat: {stat}, pvalue: {pvalue}')
    # ==================== STATS ==========================

    plt.plot(x, net_dist, color='blue', label='network', linewidth=3)
    plt.plot(x, other_dist, color='red', label='other', linewidth=3)
    plt.show()


def graph_cumsum(df_, intra_net_sno, not_intra_net_sno):

    df_copy = df_.copy(deep=True)
    df_copy.sort_values('distance', inplace=True)
    intra_net_sno = df_copy.loc[df_copy.gene_id.isin(intra_net_sno)]
    not_intra_net_sno = df_copy.loc[df_copy.gene_id.isin(not_intra_net_sno)]
    other_df = df_copy.loc[~df_copy.gene_id.isin(intra_net_sno + not_intra_net_sno)]

    # For stats----------------------
    # stats_km(net_df.distance.values, other_df.distance.values)

    # print(len(net_df), len(other_df))

    fig, ax = plt.subplots(figsize=(12, 8))
    dfs = [intra_net_sno, not_intra_net_sno, other_df]
    colors = ['#377eb8', 'purple', '#e41a1c']
    labels = ['intra snoRNA-host interacting', 'not intra snoRNA-host interacting', 'snoRNA-host not interacting']

    for df, color, label in zip(dfs, colors, labels):
        print(df)
        size = len(df)
        value = 100 / size

        x = [0]
        y = [0]
        for dist in df.distance.values:
            if dist == -1:
                continue
            if dist > 2000:
                break

            y.append(y[-1])
            x.append(dist)
            x.append(dist)
            y.append(y[-1] + value)

        plt.plot(x, y, color=color, label=label, linewidth=3)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)


    # plt.grid(b=True, which='major', color='lightgray', linestyle='-')
    plt.title('snoRNA splicing modulation potential', fontsize=25)
    plt.xlabel('Closest alternative splicing event (bp)', fontsize=20)
    plt.ylabel('Cumulative % of snoRNA', fontsize=20)
    plt.legend(fontsize=16)
    plt.savefig(graph2, format='svg')
    # plt.show()


def main():

    gtf_df = load_df(parsed_file)
    sno_host_loc_df = load_df(sno_host_loc_file)
    sno_host_df = load_df(sno_host_file)
    # better but weird if removed... #CHANGED....
    intra_sno_host_df = sno_host_df.loc[sno_host_df.interaction_type == 'intra']
    not_intra_sno_host_df = sno_host_df.loc[sno_host_df.interaction_type != 'intra']

    sno_host_dict, prot_cod_df, sno_df = get_sno_and_host(gtf_df, sno_host_loc_df)

    sno_df = get_sno_intron(sno_host_dict, prot_cod_df, sno_df, sno_host_loc_df)

    # Get stats to make a graph
    stats_df = get_stats(sno_df, intra_sno_host_df)
    graph(stats_df)

    print('===============================================================')

    graph_cumsum(sno_df, list(intra_sno_host_df.single_id1), list(not_intra_sno_host_df.single_id1))

    print('===============================================================')
    name_id_dict = dict(zip(prot_cod_df.gene_id, prot_cod_df.gene_name))
    sno_df['host_name'] = sno_df['host'].map(name_id_dict)
    sno_df['in_net'] = np.where(sno_df.gene_id.isin(intra_sno_host_df.single_id1),
                                True,
                                False)
    # print(sno_df[['seqname', 'start', 'end', 'gene_name',
    #               'host_name', 'splicing_hypothesis', 'distance', 'in_net']])

    splicing_dict = dict(zip(sno_df.gene_id, sno_df.splicing_hypothesis))
    distance_dict = dict(zip(sno_df.gene_id, sno_df.distance))
    int_start_dict = dict(zip(sno_df.gene_id, sno_df.intron_start))
    int_end_dict = dict(zip(sno_df.gene_id, sno_df.intron_end))

    intra_sno_host_df['splicing_hypothesis'] = intra_sno_host_df.single_id1.map(splicing_dict)
    intra_sno_host_df['splice_dist'] = intra_sno_host_df.single_id1.map(distance_dict)
    intra_sno_host_df['intron_start'] = intra_sno_host_df.single_id1.map(int_start_dict)
    intra_sno_host_df['intron_end'] = intra_sno_host_df.single_id1.map(int_end_dict)

    intra_sno_host_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
