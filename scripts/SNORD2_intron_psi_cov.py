import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib import style

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

psi_ext_file = snakemake.input.psi_ext

out_SNORD2 = snakemake.output.svg_SNORD2
out_SNORD2_intron = snakemake.output.svg_SNORD2_intron



def graph(df, corr_type, region, norm, out_file):

    columns = [
        x.replace('brainLam', 'HBR').replace('humanRef', 'UHR')
        for x in df.columns
    ]
    df.columns = columns

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


    tissues = set([x.split('_')[0] for x in df['sample'].values])
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

    for idx, tissue in enumerate(sorted_tissues):
        data = df.loc[df['sample'].str.contains(tissue)]

        X = data['psi'].values
        Y = data[corr_type].values

        plt.scatter(X,
                    Y,
                    label=tissue,
                    color=colors[idx],
                    edgecolor='k',
                    alpha=0.8,
                    s=100)

    all_X = df['psi'].values
    all_Y = df[corr_type].values
    coef, p_value = stats.pearsonr(all_X, all_Y)
    supp = ' extension' if region == 'Extension' else ''
    plt.title(
        f'Correlation between {norm}SNORD2{supp} and exon 4 exclusion\ncoef: {coef:.2f}, p_value: {p_value:.5f}',
        fontsize=25)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(16)

    # ------------------------------------------
    # plt the trendline
    from sklearn import linear_model
    np_all_X = np.array(all_X).reshape(-1, 1)
    regr = linear_model.LinearRegression()
    regr.fit(np_all_X, all_Y)
    trendline = regr.predict(np_all_X)
    plt.plot(all_X, trendline, color='black', linewidth=3, label='trendline', alpha=0.5, zorder=0)


    plt.xlabel('PSI value (inclusion / total)', fontsize=20)
    plt.ylabel(f'{region} coverage {norm}(raw counts)', fontsize=20)
    plt.legend(fontsize=16)

    plt.savefig(out_file, format='svg')
    # plt.show()
    plt.close()




def main():

    df = pd.read_csv(psi_ext_file, sep='\t')

    corr_types = [
        # 'ext_cov',
        'ext_cov_norm',
        # 'snoRNA_cov',
        'snoRNA_cov_norm'
    ]
    out_files = [out_SNORD2, out_SNORD2_intron]
    for corr_type, out_file in zip(corr_types, out_files):
        region = 'Extension' if 'ext' in corr_type else 'SNORD2'
        norm = 'normalized ' if 'norm' in corr_type else ''
        graph(df, corr_type, region, norm, out_file)


if __name__ == '__main__':
    main()
