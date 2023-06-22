import numpy as np
import pandas as pd

data_file = snakemake.input.intra

out_file = snakemake.output.coords


def simulate_coords(df_):

    df = df_.copy(deep=True)
    print(df.columns)
    df = df[[
        'chr1', 'start1', 'end1', 'strand1', 'chr2', 'start2', 'end2',
        'DG', 'single_id1', 'single_id2', 'name1','name2',
        'E',
        'merged_name', 'sno_start', 'sno_end',
        'intron_start', 'intron_end'
    ]]

    cols = ['sno_start', 'sno_end', 'start2', 'end2', 'intron_start', 'intron_end', 'strand1']
    simulated_start = []
    simulated_end = []
    location_type = []
    for sno_start, sno_end, start, end, intron_start, intron_end, strand in df[cols].values:
        if start > sno_start:
            offset = start - sno_end
            sim_end = sno_start - offset
            sim_start = sim_end - (end - start)
        else:
            offset = sno_start - end
            sim_start = sno_end + offset
            sim_end = sim_start + (end - start)

        simulated_start.append(sim_start)
        simulated_end.append(sim_end)

        intron_length = intron_end - intron_start
        sim_length = sim_end - sim_start
        max_length = max(intron_end, sim_end) - min(intron_start, sim_start)
        if max_length <= intron_length:
            location_type.append('intronic')
        elif max_length <= intron_length + sim_length:
            location_type.append('partially_intronic')
        else:
            location_type.append('exonic')

    df['sim_start'] = simulated_start
    df['sim_end'] = simulated_end
    df['location_type'] = location_type

    view_df = df[cols + ['sim_start', 'sim_end', 'location_type']].copy(deep=True)
    view_df.drop(columns=['strand1'], inplace=True)
    view_df['min_'] = view_df.min(axis=1)
    # print(view_df)
    for col in view_df.columns[:-2]:
        view_df[col] = view_df[col] - view_df.min_

    # print(view_df)

    return df


def generate_output_df(df):

    master_list = []
    for i in df.index:
        chr = df.at[i, 'chr1']
        strand = df.at[i, 'strand1']
        DG = df.at[i, 'DG'].split('|')[0]

        # Info for snoRNA part
        start1 = df.at[i, 'start1']
        end1 = df.at[i, 'end1']

        # Info on real interaction
        start2 = df.at[i, 'start2']
        end2 = df.at[i, 'end2']

        # Info on simulated data
        sim_start = df.at[i, 'sim_start']
        sim_end = df.at[i, 'sim_end']

        # Insert info into a df
        row = [
            DG,
            f'chr{chr}:{start1}-{end1}:{strand}',
            f'chr{chr}:{start2}-{end2}:{strand}',
            f'chr{chr}:{sim_start}-{sim_end}:{strand}'
        ]
        master_list.append(row)


    cols = ['DG', 'snoRNA', 'target', 'simulated']
    master_df = pd.DataFrame(master_list, columns=cols)

    return master_df






def main():

    df = pd.read_csv(data_file, sep='\t')

    sim_df = simulate_coords(df)

    coord_df = generate_output_df(sim_df)

    print(coord_df)
    coord_df.to_csv(out_file, sep='\t', index=False)



if __name__ == '__main__':
    main()
