import os
import pandas as pd


tpm_file = snakemake.input.NMD_tpm
transId_geneId_geneName_file = snakemake.input.transId_geneId_geneName

deseq_dir = snakemake.params.deseq_dir

out_file = snakemake.output.sign_nmd_trans


def get_transcripts_info(file):
    return pd.read_csv(file, sep='\t',
                       names=['transcript_id', 'gene_id',
                              'gene_name', 'transcript_name'])

def get_tpms(file):

    df = pd.read_csv(file, sep='\t')
    df.set_index('transcript', inplace=True)

    new_columns = []
    for col in df.columns:
        if 'ctrl' in col:
            new_columns.append('ctrl')
        else:
            new_columns.append(col[:-1].replace('_', '').replace('KD',
                                                                 '_KD').replace('res', '_res'))

    df.columns = new_columns
    df = df.transpose()
    df = df.reset_index()
    df = df.groupby("index").mean().transpose()

    return df


def load_files(dir):

    files = [x for x in os.listdir(dir) if x.endswith('.csv')]

    # filters !!!
    p_value = 0.05
    fc_treshold = 1

    dfs = {}
    for f in files:
        name = f.replace('.csv', '')
        cond1, cond2 = name.split('-')[0], name.split('-')[1]

        df = pd.read_csv(os.path.join(dir, f))
        df.rename(columns={'Unnamed: 0': 'transcript'}, inplace=True)
        df.dropna(axis=0, subset=['padj'], inplace=True)

        df = df.loc[df.padj < p_value]
        if 'res' in name:
            df = df.loc[df.log2FoldChange >= fc_treshold]
        else:
            df = df.loc[df.log2FoldChange <= -fc_treshold]

        df = df.drop(columns=['lfcSE', 'stat', 'pvalue', 'baseMean'])
        df.sort_values('padj', axis=0, inplace=True)

        dfs[name] = df

    return dfs


def add_tpms(merge_df, KD, res, tpm_df):

    merge_df['ctrl_tpm'] = merge_df.transcript.map(dict(zip(tpm_df.index,
                                                            tpm_df['ctrl'])))
    merge_df['KD_tpm'.format(KD)] = merge_df.transcript.map(dict(zip(tpm_df.index,
                                                                     tpm_df[KD])))
    merge_df['res_tpm'.format(res)] = merge_df.transcript.map(dict(zip(tpm_df.index,
                                                                      tpm_df[res])))
    return merge_df


def get_significants_trans(dfs, KD, res, tpm_df):

    df_kd = dfs['ctrl-{}'.format(KD)].copy(deep=True)
    df_res = dfs['{}-{}'.format(KD, res)].copy(deep=True)

    set_kd = set(df_kd.transcript)
    set_res = set(df_res.transcript)

    intersect = set_kd.intersection(set_res)

    # ------------------STATS-------------------
    print('######################## STATS ########################')
    print('ctrl-{}-{}'.format(KD, res))
    print('KD:', len(set_kd))
    print('Res:', len(set_res))
    print('Intersection:', len(intersect))
    print('######################## STATS ########################')
    # ------------------STATS-------------------

    merge_df = df_kd.loc[df_kd.transcript.isin(intersect)].copy(deep=True)
    merge_df.columns = ['transcript', 'log2FoldChange_KD', 'padj_KD']
    merge_df['log2FoldChange_res'] = merge_df.transcript.map(dict(zip(df_res.transcript,
                                                                      df_res.log2FoldChange)))
    merge_df['padj_res'] = merge_df.transcript.map(dict(zip(df_res.transcript,
                                                            df_res.padj)))

    merge_df = add_tpms(merge_df, KD, res, tpm_df)
    # filter non-abundant --- doesn't change anything...
    # merge_df = merge_df.loc[(merge_df['ctrl_tpm'] + merge_df['KD_tpm'] + merge_df['res_tpm']) > 3]
    merge_df['kd_rescue_name'] = '{}-{}'.format(KD, res)

    return merge_df


def get_kd_res(dfs, tpm_df, info_df):

    SMG6 = get_significants_trans(dfs, 'SMG6_KD', 'SMG6_res', tpm_df)
    SMG7 = get_significants_trans(dfs, 'SMG7_KD', 'SMG7_res', tpm_df)
    UPF1 = get_significants_trans(dfs, 'UPF1_KD', 'UPF1_res', tpm_df)
    double_SMG6 = get_significants_trans(dfs, 'double_KD', 'double_resSMG6', tpm_df)
    double_SMG7 = get_significants_trans(dfs, 'double_KD', 'double_resSMG7', tpm_df)

    master_df = pd.concat([SMG6, SMG7, UPF1, double_SMG6, double_SMG7])
    master_df['gene_id'] = master_df.transcript.map(dict(zip(info_df.transcript_id,
                                                    info_df.gene_id)))
    master_df['gene_name'] = master_df.transcript.map(dict(zip(info_df.transcript_id,
                                                    info_df.gene_name)))

    master_df['transcript_name'] = master_df.transcript.map(dict(zip(info_df.transcript_id,
                                                                     info_df.transcript_name)))

    return master_df


def main():

    info_df = get_transcripts_info(transId_geneId_geneName_file)

    tpm_df = get_tpms(tpm_file)

    dfs = load_files(deseq_dir)

    kd_res_df = get_kd_res(dfs, tpm_df, info_df)

    kd_res_df.to_csv(out_file, sep='\t', index=False)



if __name__ == '__main__':
    main()
