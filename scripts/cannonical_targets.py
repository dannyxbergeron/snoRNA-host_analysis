import numpy as np
import pandas as pd

import scipy.stats as stats

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

network_file = snakemake.input.merged_double
sno_host_loc_file = snakemake.input.sno_host_loc
sno_host_file = snakemake.input.sno_host
snodb_file = snakemake.input.snodb
bio_function_file = snakemake.input.bio_function

out_can = snakemake.output.svg_can_targets
out_bio_func = snakemake.output.svg_bio_functions
out_box_type = snakemake.output.svg_box_type

TITLE_SIZE = 25
AXIS_LABEL_SIZE = 20
AXIS_TICK_SIZE = 16
LEGEND_SIZE = 16

def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df

def fisher_test(cond1_pos, cond1_neg, cond2_pos, cond2_neg):
    print('----------------- Fisher exact test ---------------------')
    obs = np.array([[cond1_pos, cond1_neg], [cond2_pos, cond2_neg]])
    odds, p_val = stats.fisher_exact(obs)
    print(f'Odds: {odds:.2f}, pValue: {p_val:.5f}')
    print('-------------------------------------------------------\n')

def get_ids(df, network_df, loc_df):

    sno_interacting = set(df.single_id1)
    intra = set(df.loc[df.interaction_type == 'intra'].single_id1)
    diff_loc = set(df.loc[~(df.interaction_type == 'intra')].single_id1)

    others = set(network_df.single_id1)
    others = others - sno_interacting
    others = others.intersection(set(loc_df.gene_id))

    print('--> Length of intra, diff_loc and others')
    print(len(intra), len(diff_loc), len(others))
    print('-------------------------\n')

    return sno_interacting, intra, diff_loc, others

def analyse_can_targets(all_sno_interating, others, snodb_df):

    def get_percent(ids, can_sno):
        ids_set = set(ids)
        can_sno_set = set(can_sno)
        intersect = ids_set.intersection(can_sno_set)
        total = len(ids_set)
        percent_can = len(intersect) / total * 100
        can_pos = len(intersect)
        can_neg = total - can_pos
        return percent_can, can_pos, can_neg

    print('================ Canonical targets ====================')
    snodb_df = snodb_df[[
        'gene_id_annot2020', 'rrna', 'snrna'
    ]].copy(deep=True)
    snodb_df.dropna(thresh=2, inplace=True)
    canonical_snoRNA = list(snodb_df.gene_id_annot2020)

    all_sno_interating, all_pos, all_neg = get_percent(all_sno_interating, canonical_snoRNA)
    others_percent, others_pos, others_neg = get_percent(others, canonical_snoRNA)

    # Fisher exact test
    fisher_test(all_pos, all_neg, others_pos, others_neg)

    print(f'All snoRNA interacting % of canonical: {all_sno_interating:.1f}%')
    print(f'Others % of canonical: {others_percent:.1f}%')
    print('=======================================================\n')

    return ['SnoRNA interacting\nwith their host', 'Others'], [all_sno_interating, others_percent]

def graph_canonical(names, data):

    total = [[100] for x in data]
    data = [[x] for x in data]

    fig, ax = plt.subplots(figsize=(8,8))

    bar_tot = sns.barplot(data=total, color='#fc8d62')
    bar_can = sns.barplot(data=data, color='#66c2a5')

    tick_labels = [
        label
        for label, tick_label in zip(names, ax.get_xticks().tolist())
    ]
    ax.set_xticklabels(tick_labels)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)

    plt.title('Number of canonical snoRNA\nin the different groups', size=TITLE_SIZE)
    plt.xlabel('Groups', size=AXIS_LABEL_SIZE)
    plt.ylabel('Percentage of snoRNAs', size=AXIS_LABEL_SIZE)

    top_bar = mpatches.Patch(color='#fc8d62', label='Orphans')
    bottom_bar = mpatches.Patch(color='#66c2a5', label='Cannonical snoRNAs')
    plt.legend(handles=[top_bar, bottom_bar], fontsize=LEGEND_SIZE)

    # plt.show()
    plt.savefig(out_can, format='svg')
    plt.close()

def prepare_data(all_sno_interating, others, loc_df, bio_func_df):

    def put_in_df(ids, name):
        tmp = pd.DataFrame(ids, columns=['gene_id']).copy(deep=True)
        tmp['group'] = name
        return tmp

    all_interacting_df = put_in_df(all_sno_interating, 'SnoRNA interacting\nwith their host')
    others_df = put_in_df(others, 'Others')
    master_df_original = pd.concat([all_interacting_df, others_df])

    master_df = master_df_original.copy(deep=True)
    master_df['host_id'] = master_df.gene_id.map(dict(zip(loc_df.gene_id,
                                                                   loc_df.host_id)))
    master_df['host_function'] = master_df.host_id.map(dict(zip(bio_func_df.host_id,
                                                                bio_func_df.host_function)))
    # print('\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
    # print(loc_df.columns)
    # name_dict = dict(zip(loc_df.host_id, loc_df.host_name))
    # id_list = []
    # name_list = []
    # for id in master_df.loc[pd.isnull(master_df.host_function)].host_id.values:
    #     id_list.append(id)
    #     name_list.append(name_dict[id])
    #
    # print('\n'.join(id_list))
    # print('-------------')
    # print('\n'.join(name_list))
    # print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n')
    # master_df['host_function'] = master_df['host_function'].fillna('Not investigated')

    master_df_gb = master_df[['group', 'host_function', 'gene_id']]\
                        .groupby(['group', 'host_function'])\
                        .count()\
                        .reset_index()
    total_pergroup = master_df_gb[['group', 'gene_id']].groupby('group').sum().reset_index()
    master_df_gb['sum_group'] = master_df_gb.group.map(dict(zip(total_pergroup.group,
                                                                total_pergroup.gene_id)))
    master_df_gb['percentage_func'] = master_df_gb['gene_id'] / master_df_gb['sum_group'] * 100
    master_df_gb = master_df_gb.round({'percentage_func': 1})

    print(master_df_gb)
    # Fisher exact test
    for function in sorted(list(set(master_df_gb.host_function.values))):
        print(function)
        tmp = master_df_gb.loc[master_df_gb.host_function == function].copy(deep=True)
        tmp['neg'] = tmp.sum_group - tmp.gene_id
        data = tmp[['gene_id', 'neg']].values.ravel()
        fisher_test(*data)

    return master_df_gb, master_df_original

def graph_bio_functions(data_df):

    groups = ['SnoRNA interacting\nwith their host', 'Others']
    host_fct = [
        'ribosomal protein',
        'Ribosome biogenesis & translation',
        'RNA binding, processing, splicing',
        'poorly characterized',
        'Other',
        # 'Not investigated'
    ]
    colors = [
        '#80b1d3',
        '#fdb462',
        '#fb8072',
        '#8dd3c7',
        '#bebada',
        '#d9d9d9',
        '#ffffb3',
    ]


    print('\n'.join(host_fct))

    fig, ax = plt.subplots(1, figsize=(10, 8))
    x = np.arange(len(groups))  # the label locations
    width = 0.9  # the width of the bars

    bottom = np.zeros(len(groups))
    for func, color in zip(host_fct, colors):
        tmp = data_df.loc[data_df.host_function == func]
        tmp_dict = dict(zip(tmp.group, tmp.percentage_func))
        tmp_data = [
            tmp_dict[group]
            if group in tmp_dict
            else 0
            for group in groups
        ]
        label = func.capitalize() if not func.startswith('RNA') else func
        ax.bar(x, tmp_data, width, label=label, bottom=bottom, color=color)
        bottom += np.array(tmp_data)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)


    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_title('Proportions of function of snoRNA host gene\naccording to their locations', size=TITLE_SIZE)
    ax.set_ylabel('Percentage of host in fuction categories', size=AXIS_LABEL_SIZE)
    ax.set_xlabel('Groups', size=AXIS_LABEL_SIZE)
    ax.set_xticks(x)
    ax.set_xticklabels(groups)

    # ax.legend()
    plt.legend(bbox_to_anchor=([1, 1, 0, 0]), fontsize=LEGEND_SIZE)

    fig.subplots_adjust(right=0.7)
    plt.savefig(out_bio_func, format='svg')
    # plt.show()
    plt.close()

def process_box_type(merged_df_):
    merged_df = merged_df_.copy(deep=True)
    merged_df.dropna(subset=['box_type'], inplace=True) # 2 scaRNA
    data_df = merged_df.groupby(['group', 'box_type']).count().reset_index()
    total = data_df.groupby('group').sum().reset_index()
    data_df['total_pergroup'] = data_df.group.map(dict(zip(total.group, total.gene_id)))
    data_df['percent_box_type'] = data_df.gene_id / data_df.total_pergroup * 100
    data_df = data_df.round({'percent_box_type': 1})

    return data_df

def graph_box_type(data_df_):

    order_dict = {'SnoRNA interacting\nwith their host': 0, 'Others': 1}
    data_df = data_df_.copy(deep=True)
    data_df['order'] = data_df.group.map(order_dict)
    data_df.sort_values(['order'], inplace=True)
    labels = ['SnoRNA interacting\nwith their host', 'Others']

    print(data_df)
    # Fisher exact test
    values = data_df.gene_id.values
    fisher_test(*values)

    total = data_df.groupby('group').sum().reset_index()
    cd = data_df.loc[data_df.box_type == 'C/D']

    fig, ax = plt.subplots(figsize=(8,8))

    bar_tot = sns.barplot(x='group', y='percent_box_type', data=total, color='#7570b3')
    bar_can = sns.barplot(x='group', y='percent_box_type', data=cd, color='#d95f02')

    tick_labels = [
        label
        for label, tick_label in zip(labels, ax.get_xticks().tolist())
    ]
    ax.set_xticklabels(tick_labels)

    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)
    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(AXIS_TICK_SIZE)

    plt.title('Proportion of snoRNA box type\nin different groups', size=TITLE_SIZE)
    plt.xlabel('Groups', size=AXIS_LABEL_SIZE)
    plt.ylabel('Percentage of snoRNAs', size=AXIS_LABEL_SIZE)

    top_bar = mpatches.Patch(color='#7570b3', label='H/ACA')
    bottom_bar = mpatches.Patch(color='#d95f02', label='C/D')
    plt.legend(handles=[top_bar, bottom_bar], fontsize=LEGEND_SIZE)

    # plt.show()
    plt.savefig(out_box_type, format='svg')
    plt.close()

def main():

    df = load_df(sno_host_file)
    loc_df = load_df(sno_host_loc_file)
    snodb_df = load_df(snodb_file)
    network_df = load_df(network_file)

    print('df', df.columns)
    print('loc_df', loc_df.columns)
    print('snodb_df', snodb_df.columns)
    print('network_df', network_df.columns)
    print('==================================\n')

    all_sno_interating, intra, diff_loc, others = get_ids(df, network_df, loc_df)

    # Analyse all snoRNA interacting with their host together
    names, data = analyse_can_targets(all_sno_interating, others, snodb_df)

    graph_canonical(names, data)

    # Bio function graph
    print('\n============================== Bio function analysis ==============================\n')
    bio_func_df = load_df(bio_function_file)
    data_df, merged_df = prepare_data(all_sno_interating, others, loc_df, bio_func_df)

    graph_bio_functions(data_df)

    # Box type analysis
    print('\n============================== Box type analysis ==============================\n')
    merged_df['box_type'] = merged_df.gene_id.map(dict(zip(snodb_df.gene_id_annot2020,
                                                           snodb_df['box type'])))
    data_df = process_box_type(merged_df)
    graph_box_type(data_df)




if __name__ == '__main__':
    main()
