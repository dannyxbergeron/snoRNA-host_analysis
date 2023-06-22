import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

data_file = snakemake.input.sno_host
ref_file = snakemake.input.host_ref
bedgraph_file = snakemake.input.bedgraph

out_file = snakemake.output.cons

TRESH = 8


def load_df(file):
    df = pd.read_csv(file, sep='\t')
    return df

def load_bedgraph():
    df = pd.read_csv(bedgraph_file, sep='\t')
    df['crap1'] = 'crap1'
    df['crap2'] = 'crap2'
    return df


def get_intron(df_, ref_df):

    def get_loc(num):
        loc = df.at[i, f'loc{num}']
        id = df.at[i, f'ex_int_id{num}']

        if loc == 'exon':
            return np.nan, np.nan
        if loc == 'intron' or loc == 'intron_exon':
            row = ref_df.loc[ref_df.ex_id == id].values[0]
            return row[1], row[2]
        elif loc == 'exon_intron':
            ex_int_num = df.at[i, f'ex_int_num{num}']
            transcript_id = ref_df.loc[ref_df.ex_id == id].values[0][6]
            row = ref_df.loc[(ref_df.ex_num == ex_int_num) &
                             (ref_df.trans_id == transcript_id) &
                             (ref_df.ex_id.str.contains('intron'))].values[0]
            return row[1], row[2]
        else:
            return np.nan, np.nan

    df = df_.copy(deep=True)

    res = []
    for i in df.index:
        intron_start1, intron_end1 = get_loc(1)
        intron_start2, intron_end2 = get_loc(2)
        res.append([intron_start1, intron_end1, intron_start2, intron_end2])

    res_df = pd.DataFrame(res, columns=['intron_start1', 'intron_end1',
                                        'intron_start2', 'intron_end2'])

    master_df = pd.concat([df, res_df], axis=1)
    return master_df


def compute_intron_portion(df_):

    CONSERVATION_OFFSET = 2

    df = df_.copy(deep=True)
    df['int_portion_start1'] = ''
    df['int_portion_start2'] = ''

    def get_portion(idx, num):

        if pd.isna(df.at[i, f'intron_start{num}']):
            return np.nan

        tmp_start = max(df.at[i, f'start{num}'],  df.at[i, f'intron_start{num}'])
        tmp_end = min(df.at[i, f'end{num}'],  df.at[i, f'intron_end{num}'])

        sno_start = df.at[i, 'sno_start'] - CONSERVATION_OFFSET
        sno_end = df.at[i, 'sno_end'] + CONSERVATION_OFFSET

        if tmp_start > sno_start and tmp_start < sno_end:
            tmp_start = sno_end
        elif tmp_end > sno_start and tmp_end < sno_end:
            tmp_end = sno_start
        elif tmp_start < sno_start and tmp_end > sno_end:
            if (sno_start - tmp_start) >= (tmp_end - sno_end):
                tmp_end = sno_start
            else:
                tmp_start = sno_end

        if tmp_end - tmp_start > TRESH:
            tmp_start = int(tmp_start)
            tmp_end = int(tmp_end)
            return f'{tmp_start}-{tmp_end}'
        else:
            return np.nan

    for i in df.index:
        df.at[i, 'int_portion_start1'] = get_portion(i, 1)
        df.at[i, 'int_portion_start2'] = get_portion(i, 2)

    return df


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'DG', 'side', 'strand',
                'chr2', 'start2', 'end2', 'score', 'crap1', 'crap2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def get_cons(df_, bedgraph_df):

    df = df_.copy(deep=True)

    def prep_df(num):

        tmp = df.copy(deep=True)
        tmp = tmp.loc[~tmp[f'int_portion_start{num}'].isna()]
        tmp['new_start'] = [x.split('-')[0] for x in tmp[f'int_portion_start{num}']]
        tmp['new_end'] = [x.split('-')[1] for x in tmp[f'int_portion_start{num}']]
        tmp['side'] = num
        tmp = tmp[[f'chr{num}', 'new_start', 'new_end', 'DG', 'side', f'strand{num}']]
        tmp.columns = ['chr', 'start', 'end', 'DG', 'side', 'strand']
        tmp['chr'] = 'chr' + tmp['chr']

        intersect_df = bedtools(tmp, bedgraph_df)
        intersect_df['product'] = intersect_df.overlap * intersect_df.score

        wanted_cols = ['DG', 'start1', 'end1', 'product']
        tmp_intersect = intersect_df[wanted_cols].groupby('DG').agg({'start1': 'min',
                                                                    'end1':'max',
                                                                    'product':'sum'}).reset_index()
        tmp_intersect['length'] = tmp_intersect['end1'] - tmp_intersect['start1']
        tmp_intersect['mean'] = tmp_intersect['product'] / tmp_intersect['length']

        length_dict = dict(zip(tmp_intersect['DG'], tmp_intersect['length']))
        mean_dict = dict(zip(tmp_intersect['DG'], tmp_intersect['mean']))
        return length_dict, mean_dict

    sno_len_dict, sno_mean_dict = prep_df(1)
    other_len_dict, other_mean_dict = prep_df(2)

    df['sno_len_cons'] = df.DG.map(sno_len_dict)
    df['sno_mean_cons'] = df.DG.map(sno_mean_dict)

    df['other_len_cons'] = df.DG.map(other_len_dict)
    df['other_mean_cons'] = df.DG.map(other_mean_dict)

    df.sort_values(['other_mean_cons', 'sno_mean_cons'], inplace=True, ascending=False)

    return df


def create_tpm_dict(df_):

    df = df_.copy(deep=True)
    df.drop(columns=['gene_name', 'SKOV_nf_1', 'SKOV_nf_2'], inplace=True)
    df['avg'] = df.mean(axis=1)

    return dict(zip(df.gene_id, df['avg']))


def write(df_):

    df = df_.copy(deep=True)
    # df = df.loc[~((df.sno_mean_cons.isna()) & (df.other_mean_cons.isna()))]
    print(df.columns)

    df.drop(columns=['ex_int_num1', 'ex_int_id1', 'ext_pb1',
                     'ex_int_num2', 'ex_int_id2', 'ext_pb2',
                     'target_trans_id', 'target_trans_name',
                     #'merged_name', 'sno_start', 'sno_end',
                     'sno_length',
                     'intron_start1', 'intron_end1',
                     'intron_start2', 'intron_end2'], inplace=True)

    df = df.round({'sno_len_cons': 0, 'sno_mean_cons': 2, 'other_len_cons': 0,
                   'other_mean_cons': 2, 'sno_tpm': 2, 'target_tpm': 2})
    print(df)
    print(df.columns)

    df.to_csv(out_file, sep='\t', index=False)

def main():

    df = load_df(data_file)
    ref_df = load_df(ref_file)
    bedgraph_df = load_bedgraph()

    df = get_intron(df, ref_df)

    df = compute_intron_portion(df)

    df = get_cons(df, bedgraph_df)

    write(df)


if __name__ == '__main__':
    main()
