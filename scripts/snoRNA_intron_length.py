import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

sno_host_loc_file = snakemake.input.sno_host_loc
full_sno_host_file = snakemake.input.full_sno_host

svg_out = snakemake.output.svg

def graph(data_df):

    groups = [
        'Interacting snoRNA host intron',
        'Intron of other snoRNA hosts',
    ]
    colors = [
        '#9802C2',
        '#12B203',
    ]

    THRESH = 10_000
    datas = [
        data_df.loc[data_df.intra].intron_length.values,
        data_df.loc[~(data_df.intra)].intron_length.values,
    ]
    datas = [
        np.where(data > THRESH, THRESH, data)
        for data in datas
    ]


    fig, ax = plt.subplots(figsize=(12,8))

    for label, data, color in zip(groups, datas, colors):
        sns.kdeplot(data=data, shade=True, linewidth=1, alpha=.3,
                    label=label, ax=ax, bw=500,
                    color=color)

    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-3] = str(tick_labels[-3]) + '+'
    tick_labels[-2] = str('')
    ax.set_xticklabels(tick_labels)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)

    plt.title('Lenght of introns for snoRNAs interacting with\ntheir host in the same intron vs the others', fontsize=25)
    plt.xlabel('Lenght of introns', fontsize=20)
    plt.ylabel('Distribution of introns embedding snoRNAs', fontsize=20)
    legend = plt.legend(fontsize=16)

    # ============= STATS ==================
    from scipy.stats import mannwhitneyu as mwu
    print('\n=========== STATS - mann-whitney u test ==========')
    intra = datas[0]
    others = datas[1]
    odds, pval = mwu(intra, others, alternative='two-sided')
    print(f'intra introns vs other introns: odds {odds:.2f}, p_val {pval}')
    print('===================================================\n')
    # ============= STATS ==================

    plt.savefig(svg_out, format='svg')
    # plt.show()

def main():

    sno_host_df = pd.read_csv(sno_host_loc_file, sep='\t')
    df = pd.read_csv(full_sno_host_file, sep='\t')

    print(sno_host_df.columns)
    print(df.columns)

    intra_sno = set(df.loc[df.interaction_type == 'intra'].single_id1)
    print(len(intra_sno))
    sno_host_df['intra'] = [
        x in intra_sno
        for x in sno_host_df.gene_id.values
    ]

    data = sno_host_df[[
        'gene_id', 'gene_name',
        'host_id', 'host_name',
        'intron_start', 'intron_end',
        'intra'
    ]]
    data['intron_length'] = data['intron_end'] - data['intron_start']

    graph(data)

if __name__ == '__main__':
    main()
