import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']

tpm_file = snakemake.input.NMD_tpm
transId_geneId_file = snakemake.input.transId_geneId
significant_file = snakemake.input.sign_nmd_trans


def trans_of_genes(file):

    corresp_dict = {}
    with open(file, 'r') as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            corresp_dict[fields[0]] = fields[1]
    return corresp_dict


def load_df(file, corresp_dict, GOI):
    df = pd.read_csv(file, sep='\t',
                     dtype={'ctrl': float,
                               'XRN1': float,
                            'SMG6-XRN1': float,
                            'UPF1-XRN1': float})

    df['gene'] = df['transcript'].map(corresp_dict)
    df = df.loc[df.gene == GOI]
    df.set_index('transcript', inplace=True)
    df.drop('gene', axis=1, inplace=True)

    df.sort_index(inplace=True)
    print(df)

    return df


def load_process_colombo(file, corresp_dict, GOI, significant_t):


    def get_GOI(df):
        df = df[['ctrl', 'UPF1_KD', 'SMG6_KD', 'double_KD', 'UPF1_res',
                 'SMG6_res',  'double_resSMG6']]

        df['gene'] = df.index.map(corresp_dict)
        df = df.loc[df.gene == GOI]
        df.drop(columns=['gene'], inplace=True)
        df.sort_index(inplace=True)
        return df


    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.transcript.isin(significant_t)]

    # For only the transcript not having exon 4
    # df = df.loc[df.transcript.isin(['ENST00000429589', 'ENST00000467585'])]
    print(df)

    df.set_index('transcript', inplace=True)
    df = df.transpose()

    new_index = []
    for idx in df.index:
        if 'ctrl' in idx:
            new_index.append('ctrl')
        elif 'KD' in idx:
            new_index.append(idx[:-1])
        elif 'double' in idx:
            new_index.append(idx[:-2])
        else:
            new_index.append(idx[:-1])

    df.index = new_index
    df.reset_index(inplace=True)
    stdev_df = df.groupby('index').std().transpose().copy(deep=True)
    df = df.groupby('index').mean().transpose()

    wanted_cols = ['ctrl', 'UPF1_KD', 'SMG6_KD', 'double_KD',
                   'UPF1_res', 'SMG6_res',  'double_resSMG6']
    df = df[wanted_cols]
    stdev_df = stdev_df[wanted_cols]

    stdev_df = stdev_df.loc[stdev_df.index.isin(df.index)]

    print('Final df----------------------------')
    print(df)
    print('stdev df----------------------------')
    print(stdev_df)

    return df, stdev_df


def graph_col(df, stdev_df, trans_id_name_dict):

    colors = {
        'ctrl': '#252525',
        'UPF1_KD': '#33a02c',
        'SMG6_KD': '#1f78b4',
        'double_KD': '#e31a1c',
        'UPF1_res': '#b2df8a',
        'SMG6_res': '#a6cee3',
        'double_resSMG6': '#fb9a99'
    }

    fig, ax = plt.subplots(1, figsize=(16, 8))
    x = np.arange(len(df.index))  # the label locations
    width = 0.1  # the width of the bars
    start = (len(df.columns)/2)*width

    for i, col in enumerate(df.columns):
        ax.bar(x - start + width*i + 0.5*width, list(df[col]), width, label=col, color=colors[col],
               yerr=list(stdev_df[col]), error_kw=dict(ecolor='k', lw=1, capsize=2, capthick=1))

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('TPM')
    ax.set_title('Colombo 2017 - Average of 3 replicates')
    ax.set_xticks(x)
    ax.set_xticklabels([trans_id_name_dict[x] for x in df.index])
    if len(df) != 0:
        ax.legend()

    fig.subplots_adjust(bottom=0.25)

    # plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha="right", rotation_mode="anchor")


def graph(df_colombo, stdev_df, trans_id_name_dict):

    graph_col(df_colombo, stdev_df, trans_id_name_dict)
    plt.show()
    # plt.savefig('/data/labmeetings/host_interactions/NMD_sign_trans.svg',
    #             format='svg', transparent=True)


def main():

    GOI = 'EIF4A2'

    corresp_dict = trans_of_genes(transId_geneId_file)

    sign_df = pd.read_csv(significant_file, sep='\t')
    sign_GOI_df = sign_df.loc[sign_df.gene_name == GOI]
    significant_t = set(sign_GOI_df.transcript)
    trans_id_name_dict = dict(zip(sign_GOI_df.transcript,
                                  sign_GOI_df.transcript_name))

    colombo_trans_df, stdev_df = load_process_colombo(tpm_file, corresp_dict,
                                                      GOI, significant_t)

    graph(colombo_trans_df, stdev_df, trans_id_name_dict)


if __name__ == '__main__':
    main()
