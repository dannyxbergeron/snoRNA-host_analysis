from collections import defaultdict

import numpy as np
import pandas as pd

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']
sns.set_theme()
plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
plt.rc('ytick', labelsize=20)    # fontsize of the tick labels

data_file = snakemake.input.data_file

BIO_COLORS = {
    'protein_coding': '#33a02c',
    'lncRNA': '#6a3d9a',
    'tRNA': '#ff7f00',
    # 'scaRNA': '#a6cee3',
    'rRNA': '#e31a1c',
    'snoRNA': '#1f78b4',
    # 'intergenic': '#999999',
    'other': '#252a2c',
    # 'miRNA': '#a65628',
    'snRNA': '#cab2d6',
}

def load_df(file):
    return pd.read_csv(file, sep='\t')

def replace_weird_biotypes(df_):
    """Replace all weird byotypes with other."""
    df = df_.copy(deep=True)
    df['biotype2'].fillna('intergenic', inplace=True)
    # df = df.loc[df.biotype2 != 'rRNA'] # to remove rRNA !!!!!!!!!!!!!!!!

    df['simple_biotype2'] = np.where(df.biotype2.isin(BIO_COLORS.keys()),
                                    df.biotype2, 'other')
    print('===================== STATS ===================')
    total_interactions_df = df[['simple_biotype2', 'support']].groupby('simple_biotype2').agg({'support': ['sum', 'count']})
    total_interactions_df = total_interactions_df.reset_index().sort_values(('support', 'count'), ascending=False)
    print(total_interactions_df)
    print('===============================================')

    return df, total_interactions_df


def graph_pie_chart(df_):

    df = df_.copy(deep=True)
    df['color'] = df.simple_biotype2.map(BIO_COLORS)

    fig1, ax1 = plt.subplots()
    ax1.pie(df[('support', 'count')], labels=df['simple_biotype2'], labeldistance=1.03, autopct='%1.1f%%',
            shadow=False, startangle=0, pctdistance=0.8, colors=df['color'])

    # draw circle
    # centre_circle = plt.Circle((0, 0), 0.70, fc='white')
    # fig = plt.gcf()
    # fig.gca().add_artist(centre_circle)

    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.8))
    plt.title("Composition of biotypes interacting with snoRNAs")
    plt.tight_layout()

    # plt.savefig('/data/labmeetings/host_interactions/target_biotypes_count.svg', format='svg', transparent=True)
    plt.show()


def barh_snoRNA_host(df_):

    NUM_TOP_SNO = 30
    FONTSIZE = 20
    SNO_SIZE = 15

    df = df_.copy(deep=True)
    gene_id_name_dict = dict(zip(df.single_id1, df.name1))
    sno_host_df = df.loc[df.host_interaction]

    sno_list = list(sno_host_df.single_id1)
    all_sno_host = df.loc[df.single_id1.isin(sno_list)]

    top = all_sno_host[['single_id1', 'support']].groupby('single_id1').sum().reset_index()
    top.sort_values('support', inplace=True, ascending=False)
    sno_top_list = top.single_id1.values[:NUM_TOP_SNO]

    viz_df = pd.DataFrame(BIO_COLORS.keys(), columns=['biotype'])
    viz_df.set_index('biotype', drop=True, inplace=True)
    print(viz_df)
    sno_host_dict = {}
    for sno in sno_top_list:
        tmp = all_sno_host.loc[all_sno_host.single_id1 == sno]

        bio_tmp = tmp[['simple_biotype2', 'support']]
        # bio_count = bio_tmp.groupby('simple_biotype2').count().reset_index()
        bio_support = bio_tmp.groupby('simple_biotype2').sum().reset_index()
        viz_df[sno] = viz_df.index.map(dict(zip(bio_support.simple_biotype2,
                                                  bio_support.support)))

        host_tmp = tmp[['host_interaction', 'support']]
        # host_count = host_tmp.groupby('host_interaction').count().reset_index()
        host_support = host_tmp.groupby('host_interaction').sum().reset_index()
        non_host, host = host_support.support
        sno_host_dict[sno] = host / (host + non_host) * 100

    viz_df = viz_df.fillna(0).T
    viz_df['support'] = viz_df.sum(axis=1)
    viz_df['host_int'] = viz_df.index.map(sno_host_dict)
    for col in viz_df.columns[:-2]:
        print(col)
        viz_df[col] = viz_df[col] / viz_df['support'] * 100
    print(viz_df)
    viz_df.sort_values('host_int', inplace=True, ascending=True) # CHANGED !!

    r = [x for x in range(NUM_TOP_SNO)]

    fig = plt.figure(constrained_layout=True, figsize=(20, 10))
    gs = fig.add_gridspec(1, 2)
    ax = fig.add_subplot(gs[0, :1])
    ax2 = fig.add_subplot(gs[0, 1:])
    offset = np.zeros(NUM_TOP_SNO)
    offset_dict = {}
    for biotype, color in BIO_COLORS.items():
        ax.barh(r, list(viz_df[biotype]), left=offset, color=color,
                 height=0.9, label=biotype)
        offset_dict[biotype] = offset.copy()
        offset += np.array(list(viz_df[biotype]))

    # create host values
    ax2.barh(r, viz_df.host_int, color='#33393B',
             height=0.9, label='host')


    # Custom axis
    sno_top_names = [gene_id_name_dict[x] for x in sno_top_list]
    plt.ylabel('snoRNA', size=FONTSIZE)
    plt.yticks(r, sno_top_names, size=SNO_SIZE)
    plt.xlabel(r"% of the interactions", size=FONTSIZE)
    plt.margins(y=0)
    ax.set_yticks(r, sno_top_names)
    ax.set_yticklabels([])
    ax.set_xlabel(r"% of the interactions", fontdict={'size': FONTSIZE})
    ax.margins(y=0)

    # add the labels
    for i, v in enumerate(list(viz_df.support)):
        plt.text(101, i, str(v), color='blue', va='center', size=SNO_SIZE)

    plt.text(94, NUM_TOP_SNO + 0.5, 'Interaction counts',
             color='black', va='center', fontdict={'size': FONTSIZE})

    # Add a legend
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1), ncol=1)
    plt.title(f'Top {NUM_TOP_SNO} snoRNA', size=FONTSIZE+2)
    plt.margins(y=0)
    plt.tight_layout()
    # Show graphic
    # plt.savefig(f'/data/labmeetings/host_interactions/barh_top_{NUM_TOP_SNO}.svg',
    #             format='svg', transparent=True)
    plt.show()




def main():

    df = load_df(data_file)
    df = df.drop_duplicates(subset='DG')
    df, total_interactions_df = replace_weird_biotypes(df)

    # graph_pie_chart(total_interactions_df)

    barh_snoRNA_host(df)



if __name__ == '__main__':
    main()
