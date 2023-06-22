import numpy as np
import pandas as pd

tableS1 = snakemake.input.tableS1
snodb_file = snakemake.input.snodb
bio_function_file = snakemake.input.bio_function
gab_sno_host = snakemake.input.gab_sno_host


def main():

    df = pd.read_csv(tableS1)
    df = df[[
        'DG', 'snoRNA_id', 'host_id', 'snoRNA_name', 'host_name',
        'snoRNA_mean_tpm', 'host_mean_tpm', 'support', 'E',
        'target_mean_cons', 'splice_dist', 'extension_ratio'
    ]]

    snodb_df = pd.read_csv(snodb_file, sep='\t')
    snodb_df = snodb_df[[
        'gene_id_annot2020', 'gene_name_annot2020',
        'box type', 'rrna', 'snrna',
    ]]

    bio_function_df = pd.read_csv(bio_function_file, sep='\t')

    df['box_type'] = df.snoRNA_id.map(dict(zip(snodb_df.gene_id_annot2020, snodb_df['box type'])))
    df['target_rrna'] = df.snoRNA_id.map(dict(zip(snodb_df.gene_id_annot2020, snodb_df.rrna)))
    df['target_snrna'] = df.snoRNA_id.map(dict(zip(snodb_df.gene_id_annot2020, snodb_df.snrna)))
    df['can_target'] = [
        True if not pd.isnull(rrna) or not pd.isnull(snrna)
        else False
        for rrna, snrna in df[['target_rrna', 'target_snrna']].values
    ]
    df.drop(columns=['target_rrna', 'target_snrna'], inplace=True)

    df['host_function'] = df.host_id.map(dict(zip(bio_function_df.host_id, bio_function_df.host_function)))

    # Get the number of snoRNA per host
    sno_host_df = pd.read_csv(gab_sno_host)
    sno_per_host = sno_host_df[['gene_id', 'sno_id']].groupby('gene_id').count().reset_index()

    df['sno_per_host'] = df.host_id.map(dict(zip(sno_per_host.gene_id, sno_per_host.sno_id)))

    df.to_csv('data/tmp/classification_test.csv', index=False)


if __name__ == '__main__':
    main()
