import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

parsed_file = snakemake.input.parsed
sno_host_file = snakemake.input.prot_coding_sno_host
tpm_file = snakemake.input.tpm

out_file = snakemake.output.sno_host_loc
sno_all_transcripts = snakemake.output.sno_all_transcripts


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    if 'Ensembl' in file:
        df.rename(columns={'seqname':'chr'}, inplace=True)
    return df


def get_sno_and_host(gtf_df, sno_host_df):

    sno_dict = dict(zip(sno_host_df.sno_id, sno_host_df.host_id))

    sno_df = gtf_df.loc[(gtf_df.gene_id.isin(sno_host_df.sno_id))
                         & (gtf_df.feature == 'gene')]
    prot_cod_df = gtf_df.loc[gtf_df.gene_biotype == 'protein_coding']

    colnames = ['chr', 'start', 'end', 'gene_id', 'gene_name', 'strand']
    sno_df = sno_df[colnames]
    prot_cod_df = prot_cod_df.loc[prot_cod_df.gene_id.isin(sno_dict.values())]

    return sno_dict, prot_cod_df, sno_df


def create_introns(df):
    master_list = []
    df_values = df.values
    for i, row in enumerate(df_values):
        chr, start, end, exon_number, exon_id, strand = row
        master_list.append(row)
        if i != len(df) - 1:
            if exon_number < df_values[i+1][3]:
                intron_id = f'intron_{exon_id}'
                if strand == '+':
                    int_start = end
                    int_end = df_values[i + 1][1]
                else:
                    int_start = df_values[i + 1][2]
                    int_end = start
                intron_row = [chr, int_start, int_end,
                              exon_number, intron_id, strand]
                master_list.append(intron_row)
    return pd.DataFrame(master_list, columns=df.columns)


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id', 'gene_name', 'strand1',
                'chr2', 'start2', 'end2', 'exon_number', 'exon_id', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def get_sno_intron(snodb_host_dict, prot_cod_df, sno_df_):

    sno_df = sno_df_.copy(deep=True)
    ref_df = prot_cod_df.copy(deep=True)

    intron_start = []
    intron_end = []
    host_transcript_ids = []
    host_transcript_names = []
    exon_numbers  = []
    exon_ids = []
    to_remove = []
    for idx in sno_df.index:
        sno_id = sno_df.at[idx, 'gene_id']
        sno_name = sno_df.at[idx, 'gene_name']
        host_id = snodb_host_dict[sno_id]
        host_df = ref_df.loc[ref_df.gene_id == host_id]
        sno_data_start = sno_df.at[idx, 'start']
        sno_data_end = sno_df.at[idx, 'end']

        tmp = host_df.loc[(host_df.feature == 'transcript') &
                          (host_df.transcript_biotype == 'protein_coding')]

        if len(tmp) == 0:
            # no protein_coding transcript for this gene...
            to_remove.append(sno_id)
            continue

        tmp = tmp.sort_values('transcript_name')

        tmp_values = tmp[['start', 'end', 'transcript_id', 'transcript_name']].values
        for start, end, transcript_id, transcript_name in tmp_values:
            if sno_data_start > start and sno_data_end < end:
                host_transcript_name = transcript_name
                host_transcript_id = transcript_id
                break
        else:
            # snoRNA not in a protein_coding transcript...
            to_remove.append(sno_id)
            continue

        exon_df = host_df.loc[(host_df.transcript_id == host_transcript_id)
                              & (host_df.feature == 'exon')]
        exon_df = exon_df[['chr', 'start', 'end', 'exon_number', 'exon_id', 'strand']]
        exon_intron_df = create_introns(exon_df)

        bt1 = sno_df.loc[sno_df.index == idx]
        bt2 = exon_intron_df

        intersect_df = bedtools(bt1, bt2)

        # snoRNA not the same strand as the host...
        if len(intersect_df) < 1:
            to_remove.append(sno_id)
            continue

        # snoRNA in an exon
        if 'intron' not in intersect_df['exon_id'].values[0]:
            to_remove.append(sno_id)
            continue

        # snoRNA partly in an exon
        if len(intersect_df) > 1:
            to_remove.append(sno_id)
            continue

        int_start = intersect_df.start2.values[0]
        int_end = intersect_df.end2.values[0]
        exon_number = intersect_df.exon_number.values[0]
        exon_id = intersect_df.exon_id.values[0]

        intron_start.append(int_start)
        intron_end.append(int_end)
        host_transcript_ids.append(host_transcript_id)
        host_transcript_names.append(host_transcript_name)
        exon_numbers.append(exon_number)
        exon_ids.append(exon_id)

    to_remove_set = set(to_remove)
    filt_sno_df = sno_df.loc[~sno_df.gene_id.isin(to_remove_set)].copy(deep=True)
    filt_sno_df['host_id'] = filt_sno_df.gene_id.map(snodb_host_dict)
    filt_sno_df['host_name'] = filt_sno_df.host_id.map(dict(zip(ref_df.gene_id, ref_df.gene_name)))
    filt_sno_df['intron_start'] = intron_start
    filt_sno_df['intron_end'] = intron_end
    filt_sno_df['host_transcript_id'] = host_transcript_ids
    filt_sno_df['host_transcript_name'] = host_transcript_names
    filt_sno_df['exon_number'] = exon_numbers
    filt_sno_df['exon_id'] = exon_ids

    print('==================================')
    print('len original: {}, len filtered: {}'.format(len(sno_df), len(filt_sno_df)))
    original = sno_df.copy(deep=True)
    original['host_id'] = original.gene_id.map(snodb_host_dict)
    print('Original nb of hosts:{}, filtered nb: {}'.format(len(set(original.host_id)),
                                                            len(set(filt_sno_df.host_id))))
    print('==================================')

    return filt_sno_df


def create_tpm_dict(df_):

    df = df_.copy(deep=True)
    df.drop(columns=['gene_name', 'SKOV_nf_1', 'SKOV_nf_2'], inplace=True)
    df['avg'] = df.mean(axis=1)

    return dict(zip(df.gene_id, df['avg']))


def write_transcripts(sno_df, prot_cod_df):

    transcripts_df = prot_cod_df.loc[(prot_cod_df.transcript_id.isin(sno_df.host_transcript_id))
                                     & (prot_cod_df.feature == 'transcript')]

    transcripts_df['chr'] = transcripts_df['chr'].map(str)
    transcripts_df['chr'] = 'chr' + transcripts_df['chr']

    transcripts_df = transcripts_df.sort_values(['chr', 'start', 'end']).drop_duplicates()
    print(transcripts_df.columns)
    transcripts_df = transcripts_df[['chr', 'start', 'end', 'gene_id']]


    transcripts_df.to_csv(sno_all_transcripts, sep='\t', index=False, header=False)

    print(len(transcripts_df))


def main():

    gtf_df = load_df(parsed_file)
    sno_host_df = load_df(sno_host_file)


    snodb_host_dict, prot_cod_df, sno_df = get_sno_and_host(gtf_df, sno_host_df)

    sno_df = get_sno_intron(snodb_host_dict, prot_cod_df, sno_df)

    tpm_df = load_df(tpm_file)
    tpm_dict = create_tpm_dict(tpm_df)
    sno_df['sno_tpm'] = sno_df['gene_id'].map(tpm_dict)
    sno_df['target_tpm'] = sno_df['host_id'].map(tpm_dict)

    sno_df.to_csv(out_file, sep='\t', index=False)

    write_transcripts(sno_df, prot_cod_df)



if __name__ == '__main__':
    main()
