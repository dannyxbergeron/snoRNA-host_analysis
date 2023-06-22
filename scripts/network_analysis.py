from collections import defaultdict

import pandas as pd

from pybedtools import BedTool as bt

data_file = snakemake.input.merged_double_inta
clip_file = snakemake.input.merged_clip

out_file = snakemake.output.full_merge


def load_and_filter():
    """ Load the Clip data and filter for score >= 1000 or pValue == -1 (Zavolan) """
    df = pd.read_csv(clip_file, sep='\t')
    df = df.loc[(df.score == 1000) | (df.pValue == -1)]  # CHANGED
    # df = df.loc[(df.score == 1000) & (df.pValue == 400)]
    return df


def load_net_data():
    """ Load the network data and add chr to the chr column """
    df = pd.read_csv(data_file, sep='\t')
    df['chr1'] = 'chr' + df['chr1']
    df['chr2'] = 'chr' + df['chr2']

    return df


def get_sno_dict(df_, sno_ids):
    """ Create a dict({sno_id: {prot: (best_score, best_pValue)}}) and filter
        the clip df to keep only the entry that are in experiments where proteins
        are binding snoRNAs """
    def get_cor_pval(df):
        master_dict = {}
        for i in df.index:
            exp = df.at[i, 'exp'].split('_')[0]
            id = df.at[i, 'gene_id']
            score = df.at[i, 'score']
            pval = df.at[i, 'pValue']

            if id not in master_dict:
                master_dict[id] = {}

            if exp not in master_dict[id]:
                master_dict[id][exp] = (score, pval)
            else:
                master_dict[id][exp] = (max(master_dict[id][exp][0], score),
                                        max(master_dict[id][exp][1], pval))
        return master_dict

    df = df_.copy(deep=True)

    sno_df = df.loc[df.gene_id.isin(sno_ids)]
    sno_dict = get_cor_pval(sno_df)

    df = df.loc[df.exp.isin(set(sno_df.exp))]

    return df, sno_dict


def merge(df_, net_data, sno_prot):
    """ Map the clip binding sites to the network """

    def create_dict_list(df, col):
        """ Create a dict of all the proteins for a DG """
        dict_list = defaultdict(list)
        for i in df.index:
            DG = df.at[i, 'DG']
            val = df.at[i, col]
            dict_list[DG].append(val)
        return dict_list

    def keep_best(exp_dict, score_dict, pValue_dict):
        """ Create a list of tuples of the best score and pValue (prot, score, pval)
            for an interaction for each dg"""
        master_dict = {}
        for dg, prots, scores, pvals in zip(exp_dict.keys(),
                                            exp_dict.values(),
                                            score_dict.values(),
                                            pValue_dict.values()):
            p = ''
            s = 0
            pv = -1
            tmp = []
            for prot, score, pval in zip(prots, scores, pvals):
                if prot != p:
                    tmp.append((prot, score, pval))
                    p, s, pv = prot, score, pval
                else:
                    s = max(score, tmp[-1][1])
                    pv = max(pval, tmp[-1][2])
                    tmp[-1] = (prot, s, pv)
            master_dict[dg] = tmp

        return master_dict

    def remove_sno_not_in_dataset(master_dict, dg_snoid_dict, sno_prot):
        """ Remove the data for snoRNAs that are not bount by proteins,
            which means removing prot mapping to target but not snoRNA """
        new_dict = {}
        for dg, prot_info in master_dict.items():
            snoid = dg_snoid_dict[dg]
            tmp = []
            for prot, score, pv in prot_info:
                if snoid in sno_prot and prot in sno_prot[snoid]:
                    sno_score = sno_prot[snoid][prot][0]
                    sno_pv = sno_prot[snoid][prot][1]
                    tmp.append((prot, sno_score, sno_pv, score, pv))
            if tmp != []:
                new_dict[f'{dg}|{snoid}'] = tmp
        return new_dict

    df = df_.copy(deep=True)
    df['idx'] = [x for x in range(len(df))]
    df['exp'] = df['exp'].str.replace('adrenal_gland', 'adrenalGland')
    df[['clip_prot', 'clip_cell_line']] = df['exp'].str.split('_', expand=True)

    df1 = df[['chr', 'start', 'end', 'idx', 'clip_prot', 'strand']]
    df2 = net_data[['chr2', 'start2', 'end2', 'single_id1', 'DG', 'strand2']]

    df1_bt = bt.from_dataframe(df1)
    df2_bt = bt.from_dataframe(df2)
    intersect = df1_bt.intersect(df2_bt, wo=True, s=True, sorted=False)
    new_cols = ['chr', 'start1', 'end1', 'idx', 'clip_prot', 'strand',
                'chr2', 'start2', 'end2', 'single_id1', 'DG', 'strand2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})

    dg_snoid_dict = dict(zip(intersect_df.DG, intersect_df.single_id1))

    intersect_df['score'] = intersect_df.idx.map(dict(zip(df.idx, df['score'])))
    intersect_df['pValue'] = intersect_df.idx.map(dict(zip(df.idx, df['pValue'])))

    intersect_df = intersect_df[['DG', 'clip_prot', 'score', 'pValue']]
    intersect_df.drop_duplicates(inplace=True)

    exp_dict = create_dict_list(intersect_df, 'clip_prot')
    score_dict = create_dict_list(intersect_df, 'score')
    pValue_dict = create_dict_list(intersect_df, 'pValue')

    master_dict = keep_best(exp_dict, score_dict, pValue_dict)
    master_dict = remove_sno_not_in_dataset(master_dict, dg_snoid_dict, sno_prot)

    return master_dict


def main():

    # Load the clip and the full network data
    df = load_and_filter()
    net_data = load_net_data()

    # Filter the df based on the exp that contains snoRNA
    # and create a dict with {snoid: {prot: (score, pval)}}
    df, sno_prot = get_sno_dict(df, list(net_data['single_id1']))

    # Map the interactions for single_id2, and get all the protein, score
    # and pValue that map to a DG in the network data. The info was added
    # to each row in a list of tupple
    # (prot, sno_score, sno_pValue, target_score, target_pValue)
    dg_prot_mapped_dict = merge(df, net_data, sno_prot)

    net_data['DG_snoid'] = net_data['DG'] + '|' + net_data['single_id1']
    net_data['prot_info'] = net_data.DG_snoid.map(dg_prot_mapped_dict)
    net_data.drop(columns=['DG_snoid'], inplace=True)

    net_data.to_csv(out_file, sep='\t', index=False)


if __name__ == '__main__':
    main()
