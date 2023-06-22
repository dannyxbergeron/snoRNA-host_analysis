import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.sans-serif'] = ['Arial']

in_file = snakemake.input.intra_sno_hosts

out_file = snakemake.output.svg


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    print(df.columns)
    df = df[[
        'start1', 'end1', 'strand1',
        'start2', 'end2',
        'sno_start', 'sno_end',
        'DG'
    ]]
    return df

def get_distance_from_snoRNA(df_):

    df = df_.copy(deep=True)

    side = []
    dist = []
    for s1, e1, s2, e2, strand in df[['sno_start', 'sno_end', 'start2', 'end2', 'strand1']].values:
        if strand == '-':
            max_val = max(e1, e2)
            s1, e1, s2, e2 = [(v - max_val) * -1 for v in [e1, s1, e2, s2]]

        if s2 < s1:
            side.append('upstream')
            d = s1 - e2 if s1 - e2 > 0 else 0
            dist.append(d)
        else:
            side.append('downstream')
            d = s2 - e1 if s2 - e1 > 0 else 0
            dist.append(d)

    df['side'] = side
    df['dist'] = dist

    df_stats = df.loc[df.dist < 5000] # Hack...
    print('\n===================== STATS =====================')
    print(f'Average distance: {df_stats["dist"].mean():.0f}')
    print(f'Median distance: {df_stats["dist"].median():.0f}')
    print(f'stdev distance: {df_stats["dist"].values.std():.2f}')
    print('')


    df['min'] = df[['start2', 'end2', 'sno_start', 'sno_end']].min(axis=1)
    df['start2'] = df['start2'] - df['min']
    df['end2'] = df['end2'] - df['min']
    df['sno_start'] = df['sno_start'] - df['min']
    df['sno_end'] = df['sno_end'] - df['min']

    print(df[['DG', 'start2', 'end2', 'sno_start', 'sno_end', 'strand1', 'side', 'dist']].sort_values('dist', ascending=False))

    return df


def graph(df):

    TRESH = 1000 # Remove higher values as it gives a ugly plot...
    df['dist'] = np.where(df.dist > TRESH, TRESH, df.dist)


    fig, ax = plt.subplots(figsize=(12,8))
    sns.kdeplot(data=df.dist, shade=True, linewidth=1, alpha=.3,
                label='distance', ax=ax, color='#D45D4C', bw=25)


    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-2] = str(tick_labels[-2]) + '+'
    ax.set_xticklabels(tick_labels)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(16)
    for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(16)

    plt.title('Distribution of distance of the target region from the snoRNA', size=25)
    plt.xlabel('Distance of the target region from the snoRNA', size=20)
    plt.ylabel('Density of snoRNA interaction', size=20)
    # plt.show()
    plt.savefig(out_file, format='svg')


def main():

    df = load_df(in_file)

    df = get_distance_from_snoRNA(df)

    graph(df)


if __name__ == '__main__':
    main()
