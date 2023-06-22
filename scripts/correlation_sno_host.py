import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu

import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams.update({'font.size': 15})
plt.rcParams['font.sans-serif'] = ['Arial']
sns.set_theme()

data_file = snakemake.input.cons
pearson_file = snakemake.input.pearson_corr
sno_host_file = snakemake.input.sno_host_loc


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df

def process_pearson(file):

    df = load_df(file)
    pearson_dict = {
        (sno_id, other_id): (coef, pvalue)
        for sno_id,other_id,coef,pvalue in df[[
            'sno_id', 'other_id', 'coef', 'p_value'
        ]].values
    }
    return pearson_dict

def process(df, sno_host_df):

    sno_couples = [
        (sno_id, host_id)
        for sno_id, host_id in df[['single_id1', 'single_id2']].values
    ]
    other_couples = [
        (gene_id, host_id)
        for gene_id, host_id in sno_host_df[['gene_id', 'host_id']].values
        if gene_id not in list(df.single_id1)
    ]

    return sno_couples, other_couples


def analysis(net_couples, other_couples, pearson_dict):

    def stats(couples):
        coefs = []
        pValues = []
        for couple in couples:
            if couple in pearson_dict:
                coef, pvalue = pearson_dict[couple]
                coefs.append(coef)
                pValues.append(pvalue)
            else:
                pass
                # coefs.append(0)
                # pValues.append(1)
        return np.array(list(set(zip(coefs, pValues))))

    net_stats = stats(net_couples)
    other_stats = stats(other_couples)

    print(f'number of couple in network: {net_stats.shape[0]}')
    print(f'number of couple in others: {other_stats.shape[0]}')

    # ============= STATS ==================
    print('\n=========== STATS - mann-whitney u test ==========')
    coef_stats, coef_pval = mwu(net_stats[:,0], other_stats[:,0], alternative='two-sided')
    print(f'pValue for the coefficient: {coef_pval}')
    pval_stats, pval_pval = mwu(net_stats[:,1], other_stats[:,1], alternative='two-sided')
    print(f'pValue for the pvalues: {pval_pval}')
    print('===================================================\n')
    # ============= STATS ==================


    fig, ax = plt.subplots()
    fig.canvas.draw()

    # For the coefs
    sns.kdeplot(data=net_stats[:,0], shade=True, linewidth=1, alpha=.3,
                label='network coefs', ax=ax, bw_adjust=1,
                color='#377eb8')
    sns.kdeplot(data=other_stats[:,0], shade=True, linewidth=1, alpha=.3,
                label='other coefs', ax=ax, bw_adjust=1,
                color='#e41a1c')

    plt.title('snoRNA-host coeficient correlation')
    plt.ylabel('Density')
    plt.xlabel('Pearson correlation coeficient')
    plt.legend()
    plt.show()

    # For the pvalues
    fig, ax = plt.subplots()
    fig.canvas.draw()

    sns.kdeplot(data=net_stats[:,1], shade=True, linewidth=1, alpha=.3,
                label='network pvalues', ax=ax, bw_adjust=1,
                color='cyan')

    sns.kdeplot(data=other_stats[:,1], shade=True, linewidth=1, alpha=.3,
                label='other pvalues', ax=ax, bw_adjust=1,
                color='purple')

    plt.title('snoRNA-host pvalue correlation')
    plt.ylabel('Density')
    plt.xlabel('Pearson correlation pValue')
    plt.legend()
    plt.show()


def main():

    df = load_df(data_file)
    sno_host_df = load_df(sno_host_file)
    pearson_dict = process_pearson(pearson_file)

    net_couples, other_couples = process(df, sno_host_df)

    analysis(net_couples, other_couples, pearson_dict)


if __name__ == '__main__':
    main()
