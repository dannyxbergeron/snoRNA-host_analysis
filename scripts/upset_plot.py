import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from upsetplot import UpSet

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


def process(df_):

    df = df_.copy(deep=True)

    thresh = {
        'support': 3,
        'E': 0,
        'target_mean_cons': .2,
        'splice_dist': 150,
        'extension_ratio': 2
    }
    above = ['support', 'target_mean_cons', 'extension_ratio']

    for col in df.columns:
        if col in above:
            df[col] = np.where(df_[col] > thresh[col], True, False)
        else:
            df[col] = np.where(df_[col] < thresh[col], True, False)

    graph_df = df_.copy(deep=True)
    graph_df = graph_df.set_index(pd.MultiIndex.from_frame(df))

    return graph_df


def graph(graph_df):

    print(graph_df)

    upset = UpSet(graph_df, subset_size='count',
                  intersection_plot_elements=5,
                  # orientation='vertical',
                  # sort_by='cardinality',
                  facecolor='#33393B',
                  show_counts='%d')
    # upset.add_catplot(value='support', kind='strip', color='blue', alpha=0.3)
    # upset.add_catplot(value='E', kind='box', color='red')
    # upset.add_catplot(value='splice_dist', kind='bar', color='purple')
    upset.plot()
    plt.show()
    # plt.savefig(out_file,
    #             format='svg', transparent=True)


def main():

    df = load_df(in_file)

    graph_df = process(df)

    graph(graph_df)


if __name__ == '__main__':
    main()
