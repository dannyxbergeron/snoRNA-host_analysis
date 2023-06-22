import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey
import seaborn as sns

snodb_file = '/home/danx/Documents/projects/network_with_clip/data/references/snoDB1.1.1_for_hg38_Ensembl_V101_Scottlab_2020.tsv'
sno_host_file = '/home/danx/Documents/projects/snoDB/data/sno_host/final_host_gene_101.csv'
network_file = '/home/danx/Documents/projects/network_with_clip/data/processed/merged_P-L-S_double_sno.tsv'
host_interaction_file = '/home/danx/Documents/projects/network_with_clip/data/sno_host/full_sno_host.tsv'

snodb_df = pd.read_csv(snodb_file, delimiter='\t')
sno_host_df = pd.read_csv(sno_host_file)
net_df = pd.read_csv(network_file, sep='\t')
host_inteaction_df = pd.read_csv(host_interaction_file, sep='\t')


snodb_df = snodb_df[[
    'gene_id_annot2020', 'gene_name_annot2020',
]]

snodb_df['host_id'] = snodb_df.gene_id_annot2020.map(dict(zip(sno_host_df.sno_id, sno_host_df.gene_id)))
snodb_df['host_biotype'] = snodb_df.gene_id_annot2020.map(dict(zip(sno_host_df.sno_id, sno_host_df.gene_biotype)))

snodb_df = snodb_df.dropna(subset=['gene_id_annot2020'])
snodb_df_host = snodb_df.dropna(subset=['host_id'])
snodb_df_host_prot_cod = snodb_df_host.loc[snodb_df_host.host_biotype == 'protein_coding']
snodb_df_host_prot_cod_detected = snodb_df_host_prot_cod.loc[snodb_df_host_prot_cod.gene_id_annot2020.isin(net_df.single_id1)]


num_sno = len(snodb_df)
intronic_sno = len(snodb_df_host)
prot_cod_sno = len(snodb_df_host_prot_cod)
detected_prot_cod_sno = len(snodb_df_host_prot_cod_detected)
detected_interacting = len(set(host_inteaction_df.single_id1))
same_intron = len(set(host_inteaction_df.loc[host_inteaction_df.interaction_type == "intra"].single_id1))


print(f'Number of unique snoRNAs: {num_sno}')
print(f'Number of unique intronic snoRNAs: {intronic_sno}')
print(f'Number of unique intronic snoRNAs in protein_coding: {prot_cod_sno}')
print(f'Number of unique intronic snoRNAs in protein_coding and detected in P-L-S: {}')
print(f'Number of unique intronic snoRNAs in protein_coding and detected in P-L-S and interacting with their host: {detected_interacting}')
print(f'Number of unique intronic snoRNAs in protein_coding and detected in P-L-S and interacting with their host in the same intron: {same_intron}')
print('----------------------------------------')


## Graph sankey
values = [
   num_sno,
   intronic_sno,
   prot_cod_sno,
   detected_prot_cod_sno,
   detected_interacting,
   same_intron
]
 
groups = [
    'Total',
    'Intronic',
    'In protein_coding',
    'In non-coding',
    'In protein coding',
    'Filtered out',
    'Filtered',
    'Not in network',
    'In network',
    'Interaction elsewhere',
    'Same intron interaction'
]
                                        


                                     
                    