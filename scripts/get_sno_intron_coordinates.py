import numpy as np
import pandas as pd

intra_file = snakemake.input.intra_sno_hosts
others_file = snakemake.input.sno_host_loc

out_file = snakemake.output.sno_intron_coordinates
out_samples = snakemake.output.samples
out_units = snakemake.output.units
svg = snakemake.output.svg


def load_df(file):
    return pd.read_csv(file, sep='\t')


def create_intra_loc(intra_df):

    df = intra_df.copy(deep=True)

    cols = [
        'single_id1', 'chr2', 'start2', 'end2', 'sno_start',
        'sno_end', 'strand2', 'intron_start', 'intron_end'
    ]
    good = []
    bad = []
    for single_id1, chr, s2, e2, sno_start, sno_end, strand, intron_start, intron_end in df[cols].values:
        good_start = sno_start
        good_end = intron_end
        bad_start = intron_start
        bad_end = sno_end
        if s2 > sno_start:
            if strand == '+':
                side = 'downstream'
            else:
                side = 'upstream'
        else:
            good_start, good_end, bad_start, bad_end = bad_start, bad_end, good_start, good_end
            if strand == '+':
                side = 'upstream'
            else:
                side = 'downstream'

        good.append([single_id1, chr, good_start, good_end, strand, 'intra', 'good'])
        bad.append([single_id1, chr, bad_start, bad_end, strand, 'intra', 'bad'])



    header = ['sno_id', 'chr', 'si_start', 'si_end', 'strand', 'group', 'side']
    good_df = pd.DataFrame(good, columns=header)
    bad_df = pd.DataFrame(bad, columns=header)

    master_df = pd.concat([good_df, bad_df])
    return master_df

def create_others_loc(others_df):

    df = others_df.copy(deep=True)

    cols = [
        'gene_id', 'chr', 'start', 'end', 'strand',
        'intron_start', 'intron_end'
    ]
    upstream = []
    downstream = []
    for id, chr, start, end, strand, intron_start, intron_end in df[cols].values:
        if strand == '+':
            up_start = intron_start
            up_end = end
            down_start = start
            down_end = intron_end
        else:
            up_start = start
            up_end = intron_end
            down_start = intron_start
            down_end = end

        upstream.append([id, chr, up_start, up_end, strand, 'others', 'upstream'])
        downstream.append([id, chr, down_start, down_end, strand, 'others', 'downstream'])

    header = ['sno_id', 'chr', 'si_start', 'si_end', 'strand', 'group', 'side']
    upstream_df = pd.DataFrame(upstream, columns=header)
    downstream_df = pd.DataFrame(downstream, columns=header)

    master_df = pd.concat([upstream_df, downstream_df])
    return master_df

def check_len_stats(master_df_):

    import matplotlib.pyplot as plt
    import seaborn as sns

    THRESHOLD = 100_000

    df = master_df_.copy(deep=True)

    intra_good = df.loc[(df.group == 'intra')
                        & (df.side == 'good')]
    intra_bad = df.loc[(df.group == 'intra')
                        & (df.side == 'bad')]
    others_up = df.loc[(df.group == 'others')
                        & (df.side == 'upstream')]
    others_down = df.loc[(df.group == 'others')
                        & (df.side == 'downstream')]

    groups = [
        'Intra observed', 'Intra_other_side',
        'Others upstream', 'Others Downstream'
    ]
    colors = ['#4daf4a', '#e41a1c', '#377eb8', '#b15928']

    data = [intra_good, intra_bad, others_up, others_down]

    fig, ax = plt.subplots(figsize=(12,10))

    print(f'Number of snoRNA-intron kept with less or equal to {THRESHOLD} bp in length')
    for label, data, color in zip(groups, data, colors):
        tmp = data.loc[data.length <= THRESHOLD]
        data = data.copy(deep=True)
        data['length'] =  np.where(data['length'] > THRESHOLD, THRESHOLD, data['length'])
        print(f'{label}: {len(tmp)}/{len(data)} ({len(tmp)/len(data)*100:.1f}%)')
        sns.kdeplot(data=data.length, shade=True, linewidth=1, alpha=.3,
                    label=label, ax=ax, #bw=30,
                    color=color)
    print('----------------------------------------\n\n')

    # plt.xlim(-25, 1100)
    tick_labels = [
        int(tick_label)
        for tick_label in ax.get_xticks().tolist()
    ]

    tick_labels[-2] = str(tick_labels[-2]) + '+'
    ax.set_xticklabels(tick_labels)

    plt.title('Length of snoRNA-intron for the different groups', fontsize=18)
    plt.xlabel('Length (bp)', fontsize=15)
    plt.ylabel('Distribution', fontsize=15)
    legend = plt.legend(fontsize=12)

    plt.savefig(svg, format='svg')
    # plt.show()


def create_sample_and_units(df_):

    df_units = df_.copy(deep=True)

    # ['sno_id', 'chr', 'si_start', 'si_end', 'strand', 'group', 'side']
    df_units['sample'] = df_units.sno_id + '_' + df_units.group
    df_units['coordinates'] = 'chr' + df_units.chr + ':' \
                              + df_units.si_start.map(str) + '-'\
                              + df_units.si_end.map(str) + ':'\
                              + df_units.strand
    df_units['gene_id'] = np.nan
    df_units['fasta'] = np.nan
    df_units = df_units.rename(columns={'side': 'unit'})
    df_units = df_units[[
        'sample', 'unit', 'gene_id', 'coordinates', 'fasta'
    ]]

    df_units.to_csv(out_units, sep='\t', index=False)

    df_sample = df_units[['sample']]
    df_sample.to_csv(out_samples, sep='\t', index=False)





def main():

    # Load and drop duplicates (multiple sno-host interactions)
    intra_df = load_df(intra_file)
    intra_df.drop_duplicates(subset='single_id1', inplace=True)

    # Load and remove snoRNA already in intra
    others_df = load_df(others_file)
    others_df = others_df.loc[~(others_df.gene_id.isin(intra_df.single_id1))]

    print('intra_df: ', intra_df.columns)
    print('others_df', others_df.columns)

    print('------------------------\n')

    intra_sno_intron_df = create_intra_loc(intra_df)

    others_sno_intron_df = create_others_loc(others_df)

    master_df = pd.concat([intra_sno_intron_df, others_sno_intron_df])

    # To get the length of the entire sno-intron
    master_df['length'] = master_df.si_end - master_df.si_start

    # To check the length of the sno_introns
    check_len_stats(master_df)

    # Filter below 1000 bp the snoRNA-intron because higher than
    # that, it doesn't make sens to do the folding...
    master_df = master_df.loc[master_df.length < 1000]
    create_sample_and_units(master_df)


    master_df.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
