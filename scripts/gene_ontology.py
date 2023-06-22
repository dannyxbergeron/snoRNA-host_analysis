from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu as mwu

import matplotlib.pyplot as plt
import seaborn as sns

data_file = snakemake.input.cons
sno_host_loc_file = snakemake.input.sno_host_loc
bio_function_file = snakemake.input.bio_function

def load_df(file):
    return pd.read_csv(file, sep='\t')

def main():

    df = load_df(data_file)
    sno_host_df = load_df(sno_host_loc_file)
    bio_function_df = load_df(bio_function_file)

    bio_function_dict = dict(zip(bio_function_df.host_name,
                                 bio_function_df.host_function))

    net_host = list(set(df.name2))
    other_host = list(set(sno_host_df.host_name) - set(df.name2))

    # Get the missing ones
    missing = [x for x in net_host+other_host if x not in list(bio_function_df.host_name)]

    print(f'Number of missing: {len(missing)}')
    print(f'Number of host in the network: {len(set(df.single_id2))}')
    print(f'Number of host not in the network: {len(set(sno_host_df.host_id) - set(df.single_id2))}')

    filtered_net_host = [x for x in net_host if x in bio_function_dict]
    filtered_other_host = [x for x in other_host if x in bio_function_dict]

    print('------------------------------------')
    print(len(filtered_net_host))
    print(len(filtered_other_host))
    print('------------------------------------')

    print('----------------- for network info --------------------')
    counter = Counter([bio_function_dict[x] for x in filtered_net_host])
    print(counter.most_common())
    counter = Counter([bio_function_dict[x] for x in filtered_other_host])
    print(counter.most_common())





if __name__ == '__main__':
    main()
