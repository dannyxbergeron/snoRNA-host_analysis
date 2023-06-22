rule get_sno_intron_coordinates:
    """
    Get the sno_intron coordinates for snoRNA interacting with their
    host intron or get the upstream and downstream from the others
    Cuttoff = 1000 pb for snoRNA-intron
    """
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        intra_sno_hosts = join(config['path']['sno_host_data'],
                               config['file']['alt_splice']),
    output:
        samples = join(config['path']['sno_intron'], 'samples.tsv'),
        units = join(config['path']['sno_intron'], 'units.tsv'),
        sno_intron_coordinates = join(config['path']['sno_intron'],
                                      config['file']['sno_intron_coord']),
        svg = '/data/articles/SNORD2_article/svgs/snoRNA_intron_length_all.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/get_sno_intron_coordinates.py"


'''
The mfe was calculated using the samples and units
with the folowing pipline:
    https://github.com/dannyxbergeron/fold_and_vizualize_RNA
'''

FILES, = glob_wildcards("data/sno_intron/mfe/{id}.txt")

rule combine_mfe_and_plot:
    input:
        sno_intron_coordinates = join(config['path']['sno_intron'],
                                      config['file']['sno_intron_coord']),
        mfe_files = expand(
            'data/sno_intron/mfe/{id_group_side}.txt',
            id_group_side=FILES
        ),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        merged_mfe = join(config['path']['sno_intron'], config['file']['merged_mfe']),
        svg = '/data/articles/SNORD2_article/svgs/mfe_plot.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/combine_mfe.py"
