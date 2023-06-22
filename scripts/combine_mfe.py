import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

sno_intron_coord_file = snakemake.input.sno_intron_coordinates
inputs = snakemake.input.mfe_files
sno_host_loc_file = snakemake.input.sno_host_loc

out_file = snakemake.output.merged_mfe
svg_out = snakemake.output.svg


def load_files(files):

    data = []
    for file in files:
        file_ = file.replace('NR_', 'NR#').replace('cluster_', 'cluster#')
        file_name = file_.split('/')[-1]
        structure = file_name.split('.')[0]
        id, group, side = structure.split('_')
        id = id.replace('#', '_')

        with open(file, 'r') as f:
            mfe_sentence = f.read().splitlines()[0]
            mfe = float(mfe_sentence.split()[-2])

        data.append([id, group, side, mfe])

    cols = ['sno_id', 'group', 'side', 'mfe']
    df = pd. DataFrame(data, columns=cols)

    df['unique_id'] = df.sno_id + '-' + df.side
    return df


def graph(df_):

    df = df_.copy(deep=True)

    intra_good = df.loc[(df.group == 'intra')
                        & (df.side == 'good')]
    intra_bad = df.loc[(df.group == 'intra')
                       & (df.side == 'bad')]
    others_up = df.loc[(df.group == 'others')
                       & (df.side == 'upstream')]
    others_down = df.loc[(df.group == 'others')
                         & (df.side == 'downstream')]

    groups = [
        'SnoRNA intron',
        'Matched negative',
    ]
    colors = [
        '#4CB5D4',
        '#DE523E',
    ]

    datas = [
        intra_good.mfe_norm,
        np.concatenate((others_up.mfe_norm.values,
                        others_down.mfe_norm.values),
                       axis=0)
    ]

    fig, ax = plt.subplots(figsize=(12, 8))

    for label, data, color in zip(groups, datas, colors):
        sns.kdeplot(data=data, shade=True, linewidth=1, alpha=.3,
                    label=label, ax=ax, bw=.025,
                    color=color)

    # tick_labels = [
    #     int(tick_label)
    #     for tick_label in ax.get_xticks().tolist()
    # ]
    #
    # tick_labels[-2] = str(tick_labels[-2]) + '+'
    # ax.set_xticklabels(tick_labels)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)

    plt.title('MFE per bp of snoRNA-intron for the different groups', fontsize=25)
    plt.xlabel('Normalize mfe (mfe/length)', fontsize=20)
    plt.ylabel('Distribution of snoRNA-intron', fontsize=20)
    legend = plt.legend(fontsize=16)

    # ============= STATS ==================
    from scipy.stats import mannwhitneyu as mwu
    print('\n=========== STATS - mann-whitney u test ==========')
    intra_good = list(intra_good.mfe_norm)
    intra_bad = list(intra_bad.mfe_norm)
    others_up = list(others_up.mfe_norm)
    others_down = list(others_down.mfe_norm)
    others = others_up + others_down
    good_bad_stats, good_bad_pval = mwu(intra_good, intra_bad, alternative='two-sided')
    print(f'For intra_good vs intra_bad p_value: {good_bad_pval}')
    good_others_stats, good_others_pval = mwu(intra_good, others, alternative='two-sided')
    print(f'For intra_good vs others p_value: {good_others_pval}')
    print('===================================================\n')
    # ============= STATS ==================

    plt.savefig(svg_out, format='svg')
    # plt.show()


def main():

    df = pd.read_csv(sno_intron_coord_file, sep='\t')
    df['unique_id'] = df.sno_id + '-' + df.side

    mfe_df = load_files(inputs)

    sno_host_df = pd.read_csv(sno_host_loc_file, sep='\t')

    df['mfe'] = df.unique_id.map(dict(zip(mfe_df.unique_id, mfe_df.mfe)))
    df['mfe_norm'] = df.mfe / df.length

    print(sno_host_df.columns)

    df['sno_name'] = df.sno_id.map(dict(zip(sno_host_df.gene_id, sno_host_df.gene_name)))
    df.sort_values('mfe_norm', inplace=True)
    df.to_csv(out_file, sep='\t', index=False)

    graph(df)


if __name__ == '__main__':
    main()
