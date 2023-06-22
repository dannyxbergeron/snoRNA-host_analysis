import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

from pybedtools import BedTool as bt


intra_file = snakemake.input.intra_with_bp

out_file = snakemake.output.svg_branch_point


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    print(df.columns)
    return df

def prepare_bt(df_):

    df = df_.copy(deep=True)
    # Test to see if some were close...
    # df['bp_start'] = df['bp_start'] - 10
    # df['bp_end'] = df['bp_end'] + 10
    df_bt = df[[
        'chr2', 'start2', 'end2', 'single_id1', 'DG', 'strand2'
    ]]
    bp_bt = df[[
        'chr2', 'bp_start', 'bp_end', 'name1', 'name2', 'strand2'
    ]]
    return df_bt, bp_bt


def bedtools(df1, df2):

    df1 = df1.sort_values(['chr2', 'start2', 'end2'])
    df2 = df2.sort_values(['chr2', 'bp_start', 'bp_end'])

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=True)
    new_cols = ['chr', 'start', 'end', 'sno_id1', 'DG1', 'strand',
                'chr2', 'start2', 'end2', 'name1', 'name2', 'strand' 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})

    return intersect_df


def prepare_graph(df, intersect):

    print('\n=============== STATS ====================\n')
    tmp = df.drop_duplicates(subset=['single_id1', 'single_id2'])
    print(f'Number of snoRNA-host interactions: {len(tmp)}')
    print(f'Number of snoRNA interacting in the branch point region: {len(intersect)}')
    print('')

    return {
        'Interacting': len(intersect),
        'Not interacting': len(tmp) - len(intersect)
    }

def plot_piechart(data):
    fig, ax= plt.subplots(figsize=(8,8))
    ax.axis('equal')
    labels = data.keys()
    values = data.values()
    colors = ['#ef8a62', '#67a9cf']
    explode = (0.15, 0)

    def autopct(val):
        num = (val / 100) * sum(values)
        return f'{val:.1f}% ({num:.0f})'

    ax.pie(values, labels=labels, explode=explode,
           autopct=autopct, colors=colors,
           shadow=True, startangle=120,
           textprops={'color':'black', 'fontsize': 15})
    ax.set_title('Proportion of snoRNA interacting with their host\nthat also interact with the branch point region', fontsize=18)
    # plt.show()
    plt.savefig(out_file, format='svg')


def main():

    df = load_df(intra_file)

    df_bt, bp_bt = prepare_bt(df)

    intersect = bedtools(df_bt, bp_bt)
    intersect.drop_duplicates(subset=['sno_id1'], inplace=True)

    data = prepare_graph(df, intersect)

    plot_piechart(data)



if __name__ == '__main__':
    main()
