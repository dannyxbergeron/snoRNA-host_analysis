rule merge_raw_and_filter:
    input:
        raw_files = expand(join(config['path']['gab_raw'],
                                '{SRR}_sno_interaction_snobothside.coco.csv'),
                           SRR=config['raw_files'].keys())
    output:
        initial_file = join(config['path']['raw'],
                            config['file']['initial'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/merge_raw_files.py"


rule get_single_id:
    """ Get single_id in case of more that one """
    input:
        initial_file = join(config['path']['raw'],
                            config['file']['initial']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        man_curated = join(config['path']['raw'],
                           config['file']['man_curated']),
    output:
        single_id = join(config['path']['processed'],
                         config['file']['single_id'])
    params:
        multi = join(config['path']['tmp'],
                     config['file']['single_id'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_single_ID.py"


rule get_single_dg:
    """ Resolve the multiple match per dg """
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id'])
    output:
        single_id_and_dg = join(config['path']['processed'],
                                config['file']['single_id_and_dg'])
    params:
        min_length = 8,
        host_max_offset = 10,
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_single_dg_coco.py"


rule check_overlaps:
    """ Merge the windows together """
    input:
        single_id = join(config['path']['processed'],
                         config['file']['single_id_and_dg']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype'])
    output:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/check_overlaps.py"


rule add_offsets:
    """ Add the offsets (in the host gene or the target host/sno) for the
        sno and the targets """
    input:
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_windows'])
    output:
        merged_with_offset = join(config['path']['processed'],
                                  config['file']['merged_offset'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_offsets.py"


rule double_and_filter:
    """ Double the entries for snoRNA-snoRNA for further analysis and filter
        for interactions < 8 nt and remove intergenic """
    input:
        merged_windows = join(config['path']['processed'],
                              config['file']['merged_offset']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        snoDB = join(config['path']['ref'],
                     config['file']['snoDB'])
    output:
        merged_double = join(config['path']['processed'],
                             config['file']['merged_double'])
    params:
        min_length = 9
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/analyse_merged.py"


subworkflow RNAplex_analyser:
    workdir:
        "../RNAplex_analyser"
    configfile:
        "../RNAplex_analyser/config.json"


rule intaRNA:
    input:
        merged_double = join(config['path']['processed'],
                             config['file']['merged_double']),
        intaRNA = RNAplex_analyser(
            "../network_with_clip/data/processed/merged_P-L-S_double_sno_intaRNA.tsv")
    output:
        intaRNA_tok = 'data/tmp/intaRNA_tok',
        merged_double_inta = join(config['path']['processed'],
                                  config['file']['merged_double_inta']),
    shell:
        "echo 'intaRNA done !' && touch {output.intaRNA_tok}"
