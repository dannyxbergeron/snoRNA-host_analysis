from io import StringIO
from os.path import join

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import subprocess
import shlex

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']


input_file = snakemake.input.coords
fasta_dir = snakemake.params.fasta_dir
svg = snakemake.output.svg
cumsum_svg = snakemake.output.cumsum_svg


def intaRNA(df):

    def get_pred(cmd):
        cmd = shlex.split(cmd)
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        stdout = process.communicate()[0].decode('utf-8')
        data = StringIO(stdout)
        tmp = pd.read_csv(data, sep=',')
        if tmp.empty:
            row = [['dummy1', 0, 0, 'dummy2', 0, 0, np.nan, np.nan, np.nan]]
            return pd.DataFrame(row, columns=tmp.columns)
        return tmp

    cols = ['id1', 'start1', 'end1', 'id2', 'start2', 'end2', 'subseqDP', 'hybridDP', 'E', 'DG', 'type']
    master_df = pd.DataFrame([], columns=cols)
    for dg in df.DG.values:
        dg_path = join(fasta_dir, dg)
        for type in ['target', 'simulated']:
            cmd = f'IntaRNA -q {dg_path}-snoRNA.fa -t {dg_path}-{type}.fa --outMode C --outSep=,'
            tmp = get_pred(cmd)
            tmp['DG'] = dg
            tmp['type'] = type
            master_df = pd.concat([master_df, tmp])

    return master_df


def process_for_graph(df_, nan_process=''):

    df = df_.copy(deep=True)
    if nan_process == 'drop':
        df = df.dropna()
    elif nan_process == 'fill':
        df = df.fillna(0)
    else:
        pass

    target = df.loc[df.type == 'target'].E
    simulated = df.loc[df.type == 'simulated'].E
    datas = [target, simulated]
    labels = [
        'SnoRNA target',
        'Matched negative',
    ]
    # colors = ['#4daf4a', '#e41a1c', '#377eb8']
    colors = ['#F5895E', '#6CDEEB', '#F55D69']

    return datas, labels, colors


def stats(df, datas, labels, colors):

    def gb(tmp):
        type = tmp.type.values[0]
        null = tmp.loc[tmp.E.isnull()]
        print(f'Number of {type} null E values: {len(null)}')
        not_null = tmp.loc[~(tmp.E.isnull())]
        print(f'Number of {type} non-null E values: {len(not_null)}')

    target = df.loc[df.type == 'target']
    simulated = df.loc[df.type == 'simulated']

    print('\n=================== STATS =====================')
    gb(target)
    gb(simulated)
    print('-------------------- END ---------------------\n')

    ############# STATS ##########################
    from scipy.stats import mannwhitneyu as mwu
    print('\n=========== STATS - mann-whitney u test ==========')
    odd_ratio, pval = mwu(datas[1], datas[0], alternative='two-sided')
    print(f'For snoRNA targets vs matched negative p_value: {pval}, odd_ratio: {odd_ratio}')
    print('===================================================\n')

    print('\n=========== STATS - Kolmogorov-Smirnov test ==========')
    from scipy.stats import ks_2samp
    odds, pval = ks_2samp(datas[1], datas[0])
    print(f'For snoRNA targets vs matched negative p_value: {pval}, odd_ratio: {odd_ratio}')
    print('===================================================\n')
    ############# STATS ##########################


def graph_kde(datas, labels, colors):

    fig, ax = plt.subplots(figsize=(12,8))
    fig.canvas.draw()

    for label, data, color in zip(labels, datas, colors):
        sns.kdeplot(data=data, shade=True, linewidth=1, alpha=.3,
                    label=label, ax=ax, bw=2.5,
                    color=color)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)

    # plt.title('IntaRNA prediction\n(NaN removed, 60 for target and 63 for simulated)', fontsize=25)
    plt.title('IntaRNA prediction', fontsize=25)
    plt.xlabel('Molecular free energy (mfe)', fontsize=20)
    plt.ylabel('Distribution', fontsize=20)
    plt.legend(fontsize=16)
    plt.savefig(svg, format='svg')
    # plt.show()
    plt.close()


def graph_cumsum(datas, labels, colors):

    fig, ax = plt.subplots(figsize=(12, 8))

    for data, color, label in zip(datas, colors, labels):

        size = len(data)
        value = 100 / size
        data = sorted([
            x
            for x in data
            if not pd.isnull(x)
        ])

        x = [data[0]]
        y = [0]
        for mfe in data:
            y.append(y[-1])
            x.append(mfe)
            x.append(mfe)
            y.append(y[-1] + value)

        plt.plot(x, y, color=color, label=label, linewidth=3)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(16)

    plt.ylim(0, 100)
    plt.xlim(xmax=0)
    plt.title('IntaRNA prediction cummulative graph', fontsize=25)
    plt.xlabel('Molecular free energy (mfe)', fontsize=20)
    plt.ylabel('Cumulative % of snoRNA', fontsize=20)
    plt.legend(fontsize=16)
    plt.savefig(cumsum_svg, format='svg')
    # plt.show()



def main():

    df = pd.read_csv(input_file, sep='\t')

    intaRNA_df = intaRNA(df)
    print(intaRNA_df.columns)
    row = intaRNA_df.loc[intaRNA_df.DG == '273061_L0'][[
        'subseqDP',
        'hybridDP',
        'E',
        'type'
    ]]
    for col, val in zip(row.columns, row.values[0]):
        if col == 'subseqDP':
            val = val.replace('U', 'T')
        print(f'{col}: {val}')

    datas, labels, colors = process_for_graph(intaRNA_df, 'drop')

    stats(intaRNA_df, datas, labels, colors)

    graph_kde(datas, labels, colors)

    graph_cumsum(datas, labels, colors)





if __name__ == '__main__':
    main()
