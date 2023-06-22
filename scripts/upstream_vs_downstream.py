import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

in_file = snakemake.input.intra_sno_hosts

out_file = snakemake.output.svg


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    print(df.columns)
    df = df[[
        'chr1', 'start1', 'end1', 'strand1', 'chr2', 'start2', 'end2',
        'strand2', 'support', 'DG', 'single_id1', 'single_id2', 'name1',
        'name2', 'biotype2', 'offset1', 'offset2', 'E', 'loc1', 'ex_int_num1',
        'ex_int_id1', 'ext_pb1', 'loc2', 'ex_int_num2', 'ex_int_id2', 'ext_pb2',
        'target_trans_id', 'target_trans_name', 'interaction_type',
        'merged_name', 'sno_start', 'sno_end', 'sno_length'
    ]]
    return df

def process_data(df_):

    df = df_.copy(deep=True)

    cols = ['start2', 'end2', 'sno_start', 'sno_end', 'strand2']
    rev_dict = {'upstream': 'downstream', 'downstream': 'upstream'}
    sides = []
    for i, (s2, e2, sno_start, sno_end, strand) in zip(df.index, df[cols].values):
        if s2 > sno_start: side = 'downstream'
        else: side = 'upstream'

        if strand == '-': side = rev_dict[side]

        sides.append(side)

    df['side'] = sides

    # print(df[['start2', 'end2', 'sno_start', 'sno_end', 'strand2', 'side']])
    return df


def graph_bar_chart(df_):

    df = df_.copy(deep=True)
    df['target_loc'] = np.where(df.interaction_type == 'intra', 'Same intron', 'Not same intron')

    data = df[['target_loc', 'side', 'interaction_type']].groupby(['target_loc', 'side']).count().reset_index()
    total = data.groupby('target_loc').sum().reset_index().copy(deep=True)
    upstream = data.loc[data.side == 'upstream'].copy(deep=True)

    total['interaction_count'] = [i / j * 100 for i,j in zip(total['interaction_type'], total['interaction_type'])]
    upstream['interaction_count'] = [i / j * 100 for i,j in zip(upstream['interaction_type'], total['interaction_type'])]
    print(data)
    print(total)
    print(upstream)

    fig, ax = plt.subplots(figsize=(6,8))

    # plot with support
    bar_tot = sns.barplot(x='target_loc', y='interaction_count', data=total, color='#1f78b4')
    bar_upstream = sns.barplot(x='target_loc', y='interaction_count', data=upstream, color='#b2df8a')

    plt.title('Position of the interaction from the snoRNA', size=15)
    plt.xlabel('Type of interaction', size=12)
    plt.ylabel('Percentage of interactions', size=12)

    top_bar = mpatches.Patch(color='#1f78b4', label='Downstream')
    bottom_bar = mpatches.Patch(color='#b2df8a', label='Upstream')
    plt.legend(handles=[top_bar, bottom_bar])

    # plt.show()
    plt.savefig(out_file, format='svg', transparent=True)
    plt.close()

    # Test with support
    # fig, ax = plt.subplots(figsize=(12,8))
    # boxplot = sns.boxplot(x='interaction_type', y='support', data=df)
    # sns.stripplot(x = "interaction_type",
    #               y = "support",
    #               color = 'black',
    #               size = 10,
    #               alpha = 0.3,
    #               data = df)
    # ax.set_yscale('log')
    # plt.title('Support of the interaction in each group', size=15)
    # plt.xlabel('Interaction location in host gene', size=12)
    # plt.ylabel('Support of the interaction (log10)', size=12)
    # plt.show()


def graph_pie_chart(df):

    print(df.columns)
    df_groupby = df[['single_id1', 'side']].groupby('side').count().reset_index()
    data = dict(zip(df_groupby.side, df_groupby.single_id1))

    fig, ax= plt.subplots(figsize=(8,8))
    ax.axis('equal')
    labels = [x.capitalize() for x in data.keys()]
    values = data.values()
    colors = ['#1f78b4', '#b2df8a']

    # Wedge properties
    wp = { 'linewidth' : 2, 'edgecolor' : "white" }

    def autopct(val):
        num = (val / 100) * sum(values)
        return f'{val:.1f}% ({num:.0f})'

    wedges, texts, autotexts = ax.pie(values, labels=labels,
           autopct=autopct, colors=colors,
           shadow=False, startangle=120,
           textprops={'color':'black', 'fontsize': 20},
           wedgeprops=wp)
    ax.set_title('Position of the interaction from the snoRNA\nfor snoRNA interating with their host', fontsize=18)

    ax.legend(wedges, labels,
          title ="SnoRNA target position",
          loc="center left",
          bbox_to_anchor=(1, 0, 0.5, 1),
          fontsize=16)

    # plt.show()
    plt.savefig(out_file, format='svg')



def main():

    df = load_df(in_file)

    df = process_data(df)

    # graph_bar_chart(df)

    graph_pie_chart(df)


if __name__ == '__main__':
    main()
