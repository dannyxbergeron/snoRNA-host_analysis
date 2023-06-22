import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

file = snakemake.input.data_file
parsed_gtf_file = snakemake.input.parsed_gtf
sno_host_loc_file = snakemake.input.sno_host_loc

out_file = snakemake.output.sno_host
out_ref = snakemake.output.host_ref

def load_df():
    df = pd.read_csv(file, sep='\t')
    df = df.loc[df.host_interaction]
    df = df[['chr1', 'start1', 'end1', 'strand1', 'chr2', 'start2', 'end2',
       'strand2', 'support', 'DG', 'single_id1', 'single_id2', 'name1',
       'name2', 'biotype2', 'offset1', 'offset2', 'E']]
    df = df.loc[df.biotype2 == 'protein_coding']
    host_ids = list(df.single_id2)
    sno_ids = list(df.single_id1)
    return df, host_ids, sno_ids

def load_sno_host_loc(file):
    df = pd.read_csv(sno_host_loc_file, sep='\t')
    return df

def load_ref(host_ids, sno_ids):
    df = pd.read_csv(parsed_gtf_file, sep='\t', dtype={'seqname':'str'})
    df = df[['seqname', 'feature', 'gene_id', 'start', 'end', 'strand', 'transcript_name',
            'transcript_id', 'transcript_biotype', 'exon_number', 'exon_id']]
    df.rename(columns={'seqname': 'chr'}, inplace=True)
    df = df.loc[df.gene_id.isin(host_ids+sno_ids)]

    return df


def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=True, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'exon_number', 'exon_id', 'strand',
                'chr2', 'start2', 'end2', 'type', 'DG', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df

def get_pos(df_):

    def process(tmp_df):
        if len(tmp_df) == 1:
            type, ex_num, ex_id, overlap = tmp_df.values[0]
            loc = 'intron' if 'intron' in ex_id else 'exon'
            overlap = str(overlap)
        elif len(tmp_df) == 2:
            if len(set(tmp_df.exon_number)) == 1:
                loc = 'exon_intron'
            else:
                loc = 'intron_exon'
            ex_num = min(tmp_df.exon_number)
            ex_id = tmp_df.values[0][2]
            overlap = '|'.join([str(x) for x in tmp_df.overlap])

        info = (loc, ex_num, ex_id, overlap)
        return info

    df = df_[['type', 'exon_number', 'exon_id', 'overlap']].copy(deep=True)
    sno_df = df.loc[df.type == 'sno']
    target_df = df.loc[df.type == 'target']

    sno_info = process(sno_df)
    tar_info = process(target_df)

    return sno_info, tar_info


def get_intron(df, ref_df, snoHostLoc_df):

    def create_introns(df):
        master_list = []
        df_values = df.values
        for i,row in enumerate(df_values):
            chr,start,end,exon_number,exon_id,strand = row
            master_list.append(row)
            if i != len(df) -1:
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


    sno_host_dict = dict(zip(snoHostLoc_df.gene_id, snoHostLoc_df.host_id))
    sno_host_trans_id_dict = dict(zip(snoHostLoc_df.gene_id,
                                      snoHostLoc_df.host_transcript_id))
    sno_host_trans_name_dict = dict(zip(snoHostLoc_df.gene_id,
                                      snoHostLoc_df.host_name))

    df = df.loc[df.single_id1.isin(snoHostLoc_df.gene_id)]

    bt_df_ready = df[['chr2', 'start2', 'end2', 'name2', 'DG', 'strand2']].copy(deep=True)
    bt_df_ready['name2'] = 'target'
    bt_df_sno = df[['chr1', 'start1', 'end1', 'name1', 'DG', 'strand1']].copy(deep=True)
    bt_df_sno['name1'] = 'sno'
    bt_df_sno.columns = bt_df_ready.columns

    all_target_pos = []
    all_sno_pos = []
    host_transcript_info = []
    to_remove = []
    MASTER_HOST_REF = []
    for idx in df.index:
        DG = df.at[idx, 'DG']
        sno_id = df.at[idx, 'single_id1']
        host_id = df.at[idx, 'single_id2']

        test_host_id = sno_host_dict[sno_id]

        if host_id != test_host_id:
            to_remove.append(DG)
            continue

        host_transcript_id = sno_host_trans_id_dict[sno_id]
        host_transcript_name = sno_host_trans_name_dict[sno_id]
        host_df = ref_df.loc[ref_df.transcript_id == host_transcript_id]

        row = host_df.values[0]
        bed = [row[0], row[3], row[4], 0, 'transcript', row[5], row[7]]
        MASTER_HOST_REF.append(bed)

        exon_df = host_df.loc[host_df.feature == 'exon']
        exon_df = exon_df[['chr', 'start', 'end', 'exon_number', 'exon_id', 'strand']]
        exon_intron_df = create_introns(exon_df)

         # for the host_ref file
        host_ref_df = exon_intron_df.copy(deep=True)
        host_ref_df['transcript_id'] = host_transcript_id
        MASTER_HOST_REF += list(host_ref_df.values)

        data_target_bt = bt_df_ready.loc[bt_df_ready.index == idx]
        data_sno_bt = bt_df_sno.loc[bt_df_sno.index == idx]
        data_bt = pd.concat([data_target_bt, data_sno_bt])
        intersect_df = bedtools(exon_intron_df, data_bt)

        if len(intersect_df) < 2:
            # print(host_transcript_name)
            to_remove.append(DG)
            continue

        sno_pos, target_pos = get_pos(intersect_df)
        all_sno_pos.append(sno_pos)
        all_target_pos.append(target_pos)

        host_transcript_info.append([host_transcript_id, host_transcript_name])

    all_sno_df = pd.DataFrame(all_sno_pos, columns=['loc1', 'ex_int_num1','ex_int_id1', 'ext_pb1'])
    all_target_df = pd.DataFrame(all_target_pos, columns=['loc2', 'ex_int_num2','ex_int_id2', 'ext_pb2'])

    filtered_df = df.loc[~df.DG.isin(to_remove)].copy(deep=True).reset_index(drop=True)
    filtered_df = pd.concat([filtered_df, all_sno_df, all_target_df], axis=1)

    cols = ['chr', 'start', 'end', 'ex_num', 'ex_id', 'strand', 'trans_id']
    MASTER_HOST_DF = pd.DataFrame(MASTER_HOST_REF, columns=cols)
    MASTER_HOST_DF['chr'] = 'chr' + MASTER_HOST_DF['chr']
    MASTER_HOST_DF.sort_values(['chr', 'start', 'end'], inplace=True)
    MASTER_HOST_DF.drop_duplicates(inplace=True)
    MASTER_HOST_DF.to_csv(out_ref, sep='\t', index=False)

    host_transcript_info = np.array(host_transcript_info)
    filtered_df['target_trans_id'] = host_transcript_info[:,0]
    filtered_df['target_trans_name'] = host_transcript_info[:,1]

    return filtered_df


def get_types(df_):

    df = df_.copy(deep=True)
    vals = df[['loc1', 'ex_int_num1', 'ex_int_id1',
               'loc2', 'ex_int_num2', 'ex_int_id2', 'DG']].values
    type_dict = {}

    for l1,n1,id1,l2,n2,id2,dg in vals:

        if (n1 == n2 and l1 == l2) or \
           (n1+1 == n2 and l1 == 'exon_intron' and l2 == 'intron') or \
           (n1 == n2 and l1 == 'intron_exon' and l2 == 'intron'):
            type_dict[dg] = 'intra'
        elif (n1 == n2-1 and l2 == 'exon') or \
             (n1 == n2 and l2 == 'intron_exon'):
            type_dict[dg] = 'next_exon'
        elif (n1 == n2) or \
             (n1 == n2+1 and l2 == 'intron_exon'):
            type_dict[dg] = 'previous_exon'
        else:
            type_dict[dg] = 'distant'

    df['interaction_type'] = df.DG.map(type_dict)

    # =============== STATS ====================
    from collections import Counter
    counter = Counter(type_dict.values())
    print(counter.most_common())

    test = df.copy(deep=True)
    test = test.loc[test.interaction_type == 'intra']
    test['merged_name'] = test.name1 + '_' + test.name2
    print(len(set(test.merged_name)))
    # =============== STATS ====================

    return df

def analysis_intra(df_, ref_df):

    df = df_.copy(deep=True)
    # df = df.loc[df.interaction_type == 'intra']
    df['merged_name'] = df.name1 + '_' + df.name2

    ref_df = ref_df.loc[(ref_df.gene_id.isin(df.single_id1)) &
                        (ref_df.feature == 'gene')]

    df['sno_start'] = df.single_id1.map(dict(zip(ref_df.gene_id,
                                                 ref_df.start)))
    df['sno_end'] = df.single_id1.map(dict(zip(ref_df.gene_id,
                                                 ref_df.end)))
    df['sno_length'] = df.sno_end - df.sno_start

    df.to_csv(out_file, sep='\t', index=False)


def main():

    data_df, host_ids, sno_ids = load_df()
    snoHostLoc_df = load_sno_host_loc(sno_host_loc_file)

    ref_df = load_ref(host_ids, sno_ids)

    data_df = get_intron(data_df, ref_df, snoHostLoc_df)

    data_df = get_types(data_df)

    analysis_intra(data_df, ref_df)
    print(data_df)
    print(data_df.columns)


if __name__ == '__main__':
    main()
