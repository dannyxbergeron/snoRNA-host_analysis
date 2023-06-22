from collections import Counter

import numpy as np
import pandas as pd

from pybedtools import BedTool as bt

gene_bed = snakemake.input.gene_bed
merged_double = snakemake.input.merged_double
sno_db_file = snakemake.input.snodb

out_file = snakemake.output.prot_coding_sno_host

def bedtools(df1, df2):

    first = bt.from_dataframe(df1)
    second = bt.from_dataframe(df2)
    intersect = first.intersect(second, wo=True, s=False, sorted=False)
    new_cols = ['chr1', 'start1', 'end1', 'gene_id1', 'gene_name1', 'strand1', 'biotype1',
                'chr2', 'start2', 'end2', 'gene_id2', 'gene_name2', 'strand2', 'biotype2', 'overlap']
    intersect_df = intersect.to_dataframe(names=new_cols, index_col=False,
                                          dtype={'chr': str, 'chr2': str})
    return intersect_df


def remove_dups(df, network_ids, snodb_df_, id_name_dict, id_biotype_dict):

    snodb_df = snodb_df_.copy(deep=True)
    snodb_df_pc = snodb_df.loc[~(pd.isnull(snodb_df['host gene id']))
                               & (snodb_df.biotype == 'protein_coding')]
    snodb_dict_pc = dict(zip(snodb_df_pc.gene_id_annot2020, snodb_df_pc['host gene id']))

    snodb_host_dict = {
        k:v for (k,v) in snodb_dict_pc.items()
        if v in id_biotype_dict
    }
    print(snodb_host_dict['ENSG00000201823'])

    snodb_df_other = snodb_df.loc[~(pd.isnull(snodb_df['host gene id']))
                                    & ~(snodb_df.biotype == 'protein_coding')]
    sno_db_other_list = list(snodb_df_other.gene_id_annot2020)

    counter = Counter(list(df.gene_id1))
    single = [
        id for (id, val) in counter.most_common()
        if val == 1
        and id not in snodb_host_dict
        and id not in sno_db_other_list
    ]
    single_df = df.loc[df.gene_id1.isin(single)]
    single_dict = dict(zip(single_df.gene_id1, single_df.gene_id2))

    multiple = [
        id for (id, val) in counter.most_common()
        if val > 1
        and id not in snodb_host_dict
        and id not in sno_db_other_list
    ]

    multiple_df = df.loc[df.gene_id1.isin(multiple)]

    good_host = {}
    for m_id in multiple:
        tmp = multiple_df.loc[multiple_df.gene_id1 == m_id]

        intersection = network_ids.intersection(set(tmp.gene_id2))
        if len(intersection) == 1 and m_id not in snodb_host_dict:
            # print(m_id, id_name_dict[m_id])
            good_host[m_id] = list(intersection)[0]

        else:
            pass

    man_curated = {
            'snoDB1210': 'ENSG00000042429',
            'ENSG00000207523': 'ENSG00000122406',
            'ENSG00000252050': 'ENSG00000156313',
            'ENSG00000238961': 'ENSG00000132846',
            'snoDB621': 'ENSG00000128245',
            'ENSG00000252213': None,
            'snoDB286': 'ENSG00000154710',
            'snoDB1003': 'ENSG00000172738',
            'snoDB1929': 'ENSG00000260456',
            'snoDB1753': 'ENSG00000138032',
            'ENSG00000288603': 'ENSG00000204152',
            'ENSG00000207109': 'ENSG00000272886',
            'ENSG00000201384': 'ENSG00000165548',
            'snoDB598': 'ENSG00000092871',
            'snoDB1065': 'ENSG00000070759',
            'ENSG00000200959': None,
            'snoDB257': 'ENSG00000196975',
            'ENSG00000206913': 'ENSG00000214517',
            'snoDB1662': 'ENSG00000171488',
            'ENSG00000206680': 'ENSG00000122406',
            'ENSG00000238832': 'ENSG00000075790',
            'ENSG00000206755': 'ENSG00000080603',
            'ENSG00000207496': 'ENSG00000144713',
            'NR_145781': 'ENSG00000108107',
            'ENSG00000252129': 'ENSG00000241322',
            'snoDB1689': 'ENSG00000111490',
            'ENSG00000221125': 'ENSG00000161640',
            'NR_145797': 'ENSG00000002822',
            'snoDB765': 'ENSG00000165156',
    }

    good_host.update(snodb_host_dict)
    good_host.update(single_dict)
    good_host.update(man_curated)
    good_host = {k:v for (k,v) in good_host.items() if v is not None}

    res_df = pd.DataFrame.from_dict(good_host, orient='index', columns=['host_id'])
    res_df = res_df.reset_index().rename(columns={'index':'sno_id'})
    res_df['sno_name'] = res_df.sno_id.map(id_name_dict)
    res_df['host_name'] = res_df.host_id.map(id_name_dict)

    res_df = res_df[['sno_name', 'sno_id', 'host_name', 'host_id']]

    return res_df


def main():

    ref_df = pd.read_csv(gene_bed, sep='\t')
    data_df = pd.read_csv(merged_double, sep='\t')
    snodb_df = pd.read_csv(sno_db_file, sep='\t')

    network_ids = set(data_df.single_id2)
    id_biotype_dict = dict(zip(ref_df.gene_id, ref_df.gene_biotype))
    id_name_dict = dict(zip(ref_df.gene_id, ref_df.gene_name))

    sno_df = ref_df.loc[ref_df.gene_biotype == 'snoRNA']
    prot_coding = ref_df.loc[ref_df.gene_biotype == 'protein_coding']

    intersect_df = bedtools(sno_df, prot_coding)

    result = remove_dups(intersect_df, network_ids, snodb_df,
                         id_name_dict,id_biotype_dict)

    result['biotype'] = result.host_id.map(id_biotype_dict)
    result = result.loc[result.biotype == 'protein_coding']
    result.drop(columns=['biotype'], inplace=True)

    result.to_csv(out_file, sep='\t', index=False)




if __name__ == '__main__':
    main()
