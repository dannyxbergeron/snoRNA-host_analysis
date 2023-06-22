import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

in_file = snakemake.input.tableS1

out_file = snakemake.output.svg


def load_df(file):
    df = pd.read_csv(file)
    df = df[[
        'support', 'E', 'target_mean_cons',
        'splice_dist', 'extension_ratio'
    ]]
    return df


def graph(df):

    fig, ax = plt.subplots(figsize=(12,8))
    sns.kdeplot(data=df.support, shade=True, linewidth=1, alpha=.3,
                label=None, ax=ax, color='green', bw=.85)

    plt.title('Distribution of support for snoRNA interactions', size=15)
    plt.xlabel('Support of the interaction', size=12)
    plt.ylabel('Density of snoRNA interaction', size=12)
    plt.legend([])
    # plt.show()
    plt.savefig(out_file, format='svg', transparent=True)


def main():

    df = load_df(in_file)

    graph(df)


if __name__ == '__main__':
    main()
