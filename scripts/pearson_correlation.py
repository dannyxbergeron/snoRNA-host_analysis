import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import style
from statistics import mean, stdev

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

tmp_file = snakemake.input.tpm

out_svg = snakemake.output.svg


def tpm_matrix():
    """get the tmp matrix of all tissues"""
    matrix_df = pd.read_csv(tmp_file, sep='\t')

    print(matrix_df.columns)
    matrix_df = matrix_df[['gene_id', 'gene_name',
        'HCT116_1', 'HCT116_2',
        # 'INOF_1',
        'MCF7_1', 'MCF7_2',
        'PC3_1', 'PC3_2',
        'SKOV_frg_1', 'SKOV_frg_2',
        # 'SKOV_nf_1', 'SKOV_nf_2',
        'TOV112D_1', 'TOV112D_2',
        'Liver_1', 'Liver_2', 'Liver_3',
        'Breast_1', 'Breast_2', 'Breast_3', #'Breast_4',
        'Testis_1', 'Testis_2', 'Testis_3',
        'Prostate_1', 'Prostate_2', 'Prostate_3',
        'Ovary_1', 'Ovary_2', 'Ovary_3',
        'SkeletalMuscle_1', 'SkeletalMuscle_2', 'SkeletalMuscle_3',
        'Brain_1', 'Brain_2', 'Brain_3',
        'humanRef_1', 'humanRef_2', 'humanRef_3',
        'brainLam_1', 'brainLam_2', 'brainLam_3'
    ]]
    return matrix_df


def get_values(matrix, id):

    row = matrix.loc[matrix.gene_id == id].values[0]
    row_values = row[3:]
    row_name = row[1]
    return row_name, row_values


def read_test_file():

    file_test = '/home/danx/Documents/projects/snoRNA_network_V2/raw_data/' + \
        'pearson_corr/test.tsv'
    couples = []
    with open(file_test, 'r') as f:
        for line in f.read().splitlines():
            couples.append(line.split('\t'))
    return couples


def process_ids(matrix_df, couples):

    for couple in couples:

        sno_id = couple[0]
        other_id = couple[1]

        sno_name, sno = get_values(matrix_df, sno_id)
        other_name, other = get_values(matrix_df, other_id)

        if mean(sno) < 10 or mean(other) < 10:
            status = 'out'
        else:
            status = 'ok'

        sorted_coord = sorted(list(zip(sno, other)), key=lambda x: x[0])

        sorted_sno, sorted_other = zip(*sorted_coord)

        coef, p_value = stats.pearsonr(sno, other)
        s = '{}  \t{}    \t{:.2f}\t{}\t{}'.format(
            sno_name, other_name, coef, p_value, status)
        print(s)


def get_new_EIF4A1():
    file = '/data/heavy_stuff/bigwig/EIF4A2_1_bedgraphs/Estimated_EIF4A1_tpm_NORMAL_TISSUES.csv'
    df = pd.read_csv(file)
    return list(df['EIF4A1_est_TPM'])


def graph_single(matrix_df, first_id, second_id_list):
    columns = [
        x.replace('brainLam', 'HBR').replace('humanRef', 'UHR')
        for x in matrix_df.columns
    ]
    matrix_df.columns = columns

    colors = [
        '#33a02c',
        '#fb9a99',
        '#1f78b4',
        '#6a3d9a',
        '#ff7f00',
        '#b15928',
        '#e31a1c',
        '#a6cee3',
        '#b2df8a',
        '#fdbf6f',
        '#cab2d6',
        '#ffff99',
        'black',
        'grey',
        '#b3de69',
        '#ccebc5',
        '#fccde5',
        '#fb9a99'
    ]

    tissues = set([x.split('_')[0] for x in matrix_df.columns[3:]])
    sorted_tissues = sorted(list(tissues))

    # To get the wanted order
    sorted_tissues = [
        'Brain', 'Breast', 'Liver',
        'Ovary', 'Prostate', 'SkeletalMuscle',
        'Testis',
        'HCT116', 'MCF7', 'PC3',
        'SKOV', 'TOV',
        'HBR', 'UHR'
    ]

    fig, ax = plt.subplots(figsize=(12,8))

    for other_id in second_id_list:
        sno_name, sno = get_values(matrix_df, first_id)
        other_name, other = get_values(matrix_df, other_id)

        # sorted_coord = sorted(list(zip(sno, other)), key=lambda x: x[0])
        # sorted_sno, sorted_pc = zip(*sorted_coord)

        for idx, tissue in enumerate(sorted_tissues):
            X = []
            Y = []
            for i, sample in enumerate(matrix_df.columns[3:]):
                if tissue in sample:
                    X.append(sno[i])
                    Y.append(other[i])
            plt.scatter(X, Y, label=tissue, color=colors[idx], edgecolor='k', alpha=0.8, s=100)

        coef, p_value = stats.pearsonr(sno, other)
        plt.title(
            'correlation graph PEARSON TISSUES\n coef: {:.2f}\n p_value: {:.5f}'.format(coef, p_value),
            fontsize=25)

        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(16)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(16)

        # sdev = stdev(sorted_sno)

        # plt.axvline(x=mean(sorted_sno), color='r')
        # plt.axvline(x=mean(sorted_sno) + sdev/2, color='b')
        # plt.axvline(x=mean(sorted_sno) - sdev/2, color='b')

        # ------------------------------------------
        # plt the trendline
        from sklearn import linear_model
        np_sno = np.array(sno).reshape(-1, 1)
        regr = linear_model.LinearRegression()
        regr.fit(np_sno, other)
        trendline = regr.predict(np_sno)
        plt.plot(sno, trendline, color='black', linewidth=3, label='trendline', alpha=0.5, zorder=0)
        # ------------------------------------------

        # Adding the cell_lines in the graph
        # for i, cl in enumerate(matrix_df.columns[3:]):
        #     rep = cl.split('_')[1]
        #     plt.text(sno[i], other[i], rep, fontsize=6)

        plt.xlabel(f'{sno_name} expression (tpm)', fontsize=20)
        plt.ylabel(f'{other_name} expression (tpm)', fontsize=20)
        plt.legend(fontsize=16)

        # plt.gcf().set_size_inches(10, 6.5)

        # plt.tight_layout()
        plt.savefig(out_svg, format='svg')
        # plt.show()
        plt.close()




def main():

    matrix_df = tpm_matrix()

    x_axis = 'ENSG00000156976'
    y_axis = ['ENSG00000238942']

    graph_single(matrix_df, x_axis, y_axis)


if __name__ == '__main__':
    main()
