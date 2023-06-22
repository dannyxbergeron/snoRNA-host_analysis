import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from pybedtools import BedTool as bt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 10})
plt.rcParams['font.sans-serif'] = ['Arial']

bg_file = snakemake.input.EIF4A3_CLIP
ref_file = snakemake.input.gene_bed_biotype

TRESHOLD = 300 # TODO what they used in the article


def load_df(file):

    df = pd.read_csv(file, sep='\t', names=['chr', 'start', 'end', 'score'])
    df['protein'] = 'EIF4A3'
    df['strand'] = np.where(df['score'] >= 0, '+', '-')
    df['score'] = np.where(df['score'] >= 0, df['score'], df['score'] * -1)
    df = df.loc[df.score > TRESHOLD]
    return df

def get_ref(file):
    ref_df = pd.read_csv(file, sep='\t')
    ref_df['chr'] = 'chr' + ref_df['chr']
    return ref_df

def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = [
        'chr1', 'start1', 'end1', 'score', 'protein', 'strand1',
        'chr2', 'start2', 'end2', 'gene_id', 'gene_name', 'strand2',
        'gene_biotype', 'overlap'
    ]
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    wanted_cols = [
        'chr1', 'start1', 'end1', 'score', 'gene_id',
        'gene_name', 'gene_biotype', 'strand1'
    ]
    new_names = [
        'chr', 'start', 'end', 'score', 'gene_id',
        'gene_name', 'gene_biotype', 'strand'
    ]
    intersect_df = intersect_df[wanted_cols]
    intersect_df.columns = new_names
    return intersect_df

def process(df_):

    df = df_.copy(deep=True)
    df = df.loc[df.gene_biotype == 'snoRNA']
    df_gb = df[['gene_id', 'score']].groupby('gene_id').max().reset_index()

    snord2_val = df_gb.loc[df_gb.gene_id == 'ENSG00000238942'].values[0][1]

    return df_gb, snord2_val


def graph(all_data, unfiltered_data, val):

    MAX_VAL = 10000
    all = [
        x if x <= MAX_VAL else MAX_VAL
        for x in all_data
    ]
    data = [
        x if x <= MAX_VAL else MAX_VAL
        for x in unfiltered_data
    ]

    fig, ax = plt.subplots()
    fig.canvas.draw()
    plt.plot([], label='SNORD2', color='red')
    sns.kdeplot(data=data, shade=True, linewidth=1, alpha=.3,
                label='snoRNA', ax=ax, bw_adjust=1,
                color='purple')
    sns.kdeplot(data=all, shade=True, linewidth=1, alpha=.3,
                label='all', ax=ax, bw_adjust=1,
                color='orange')

    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-3] = str(tick_labels[-3]) + '+'
    ax.set_xticklabels(tick_labels)

    plt.vlines(val, 0, 0.0006, color='red')

    plt.title('Distribution of max peak score for the best 165 snoRNAs\n(peak max > 75)')
    plt.xlabel('Max score per snoRNA')
    plt.legend()
    # plt.savefig('/data/labmeetings/host_interactions/EIF4A3_density.svg', format='svg')
    plt.show()


def main():

    df = load_df(bg_file)
    ref_df = get_ref(ref_file)

    intersect_df = bedtools(df, ref_df)

    processed_df, snord2_val = process(intersect_df)

    # snoRNA profile
    graph(df.score, processed_df.score, snord2_val)

    # all profile
    # graph(df.score, snord2_val)




if __name__ == "__main__":
    main()
