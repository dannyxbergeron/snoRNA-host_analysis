import pandas as pd

files = snakemake.input.raw_files
corresp_dict = snakemake.config['raw_files']

out_file = snakemake.output.initial_file

dfs = []
for file in files:
    exp_name = file.split('/')[-1].split('_')[0]
    exp_tag = corresp_dict[exp_name]

    df = pd.read_csv(file, dtype={'DG': str})
    df['DG'] = df['DG'] + '_' + exp_tag

    dfs.append(df)

master_df = pd.concat(dfs, ignore_index=True, sort=False)
master_df = master_df.loc[master_df.geometric_score >=0.01]

master_df.to_csv(out_file, index=False)
