import numpy as np
import pandas as pd

bp_file = snakemake.input.bp_distance
intra_file = snakemake.input.intra

out_best_bp = snakemake.output.best_branch_points
out_intra_wit_bp = snakemake.output.intra_with_bp


rm_quote = lambda x: x.replace('"', '')

def load_df(file):
    df = pd.read_csv(file, engine='python',
                     converters={'\"j\"': rm_quote,
                     '\"x\"': rm_quote})
    df = df.rename(columns=rm_quote)
    return df

def keep_best_bp(df_):

    df = df_.copy(deep=True)
    print(df.columns)
    df.sort_values(['id', 'branchpoint_prob'], ascending=False, inplace=True)
    df.drop_duplicates(subset='id', inplace=True)
    return df


def main():

    df = load_df(bp_file)

    df = keep_best_bp(df)

    df['bp_start'] = df.test_site - 5
    df['bp_end'] = df.test_site + 6

    df.to_csv(out_best_bp, sep='\t', index=False)

    intra_df = pd.read_csv(intra_file, sep='\t')
    intra_df['bp_start'] = intra_df.single_id1.map(dict(zip(df.id, df.bp_start)))
    intra_df['bp_end'] = intra_df.single_id1.map(dict(zip(df.id, df.bp_end)))


    intra_df.to_csv(out_intra_wit_bp, sep='\t', index=False)


if __name__ == '__main__':
    main()
