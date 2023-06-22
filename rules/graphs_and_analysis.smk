rule viz:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf']),
        multiple_list = join(config['path']['ref'],
                             config['file']['multiple_gene_name']),
        gene_id = join(config['path']['ref'],
                       config['file']['name_id']),
        sno_host = join(config['path']['ref'],
                        config['file']['sno_host']),
        trans_tpm = join(config['path']['count_matrices'],
                         config['file']['kallisto_trans_tpm']),
    output:
        'viz.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/viz_transcripts.py"


rule clip_candidates:
    input:
        full_merged = join(config['path']['processed'],
                           config['file']['full_network_clip'])
    output:
        tmp_candidates = join(config['path']['tmp'],
                              config['file']['clip_candidates']),
        bed_file = join(config['path']['tmp'],
                        config['file']['clip_candidates_bed'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/clip_candidates.py"


rule get_stats:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
    output:
        'stats.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/stats_and_graphs.py"


rule pearson_corr:
    input:
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf']),
    output:
        pearson_corr = join(config['path']['count_matrices'],
                            config['file']['pearson'])
    threads:
        6
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/multiprocessing_pearson.py"


rule correlation_sno_host:
    input:
        pearson_corr = join(config['path']['count_matrices'],
                            config['file']['pearson']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        'correlation.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/correlation_sno_host.py"


rule gene_ontology:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons']),
        bio_function = join(config['path']['ref'],
                            config['file']['bio_function']),
    output:
        'gene_ontology.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/gene_ontology.py"


rule graph_merged:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
    output:
        'graph_merged.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/graph_merged.py"

rule conservation_simulation_graph:
    output:
        svg = '/data/articles/SNORD2_article/svgs/other_conservation.svg',
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/conservation_simulation_graph.py"


rule sankey:
    input:
        sno_host = join(config['path']['ref'],
                        config['file']['prot_coding_sno_host']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        gene_bed = join(config['path']['ref'],
                        config['file']['gene_bed_biotype']),
        full_network_interactions =  join(config['path']['processed'],
                                          config['file']['merged_double']),
        net_sno_host = join(config['path']['sno_host_data'],
                            config['file']['sno_host_data']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/sankey.svg'
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/sankey.py"


rule get_significant_NMD:
    input:
        transId_geneId_geneName = join(config['path']['ref'],
                                       config['file']['transId_geneId_geneName']),
        NMD_tpm = join(config['path']['count_matrices'],
                       config['file']['NMD_tpm_transcript']),
    output:
        sign_nmd_trans = join(config['path']['NMD'],
                              config['file']['sign_trans_nmd']),
    params:
        deseq_dir = join(config['path']['NMD'], 'DESeq2'),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_all_significants_transcripts.py"


rule NMD_graph:
    input:
        NMD_tpm = join(config['path']['count_matrices'],
                       config['file']['NMD_tpm_transcript']),
        transId_geneId = join(config['path']['ref'],
                              config['file']['transId_geneId_geneName']),
        sign_nmd_trans = join(config['path']['NMD'],
                              config['file']['sign_trans_nmd']),
    output:
        'NMD.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/NMD.py"


rule EIF4A3_analysis:
    input:
        EIF4A3_CLIP = join(config['path']['tmp'],
                           config['file']['EIF4A3_CLIP']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
    output:
        'EIF4A3_analysis.tok'
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/EIF4A3_analysis.py"


rule tableS1:
    input:
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
    output:
        tableS1 = join(config['path']['tmp'],
                       config['file']['tableS1']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/table_S1.py"


rule upset_plot:
    input:
        tableS1 = join(config['path']['tmp'],
                       config['file']['tableS1']),
    output:
        svg = '/data/labmeetings/host_interactions/upset_plot.svg'
    conda:
        "../envs/upset_plot.yaml"
    script:
        "../scripts/upset_plot.py"


# New figures asked by Michelle -------------------------------------------------------------
rule dist_num_interactions:
    input:
        tableS1 = join(config['path']['tmp'],
                       config['file']['tableS1']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/dist_num_interactions.svg'
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/dist_num_interactions.py"

rule distance_target_from_snoRNA:
    input:
        intra_sno_hosts = join(config['path']['sno_host_data'],
                               config['file']['alt_splice']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/distance_target_from_snoRNA.svg'
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/distance_target_from_snoRNA.py"


rule upstream_vs_downstream:
    input:
        intra_sno_hosts = join(config['path']['sno_host_data'],
                               config['file']['sno_host_data']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/upstream_vs_downstream.svg'
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/upstream_vs_downstream.py"

rule cannonical_targets:
    """ For now it's same intron vs different locations vs other in network
        **Only intronic snoRNA considered"""
    input:
        merged_double = join(config['path']['processed'],
                             config['file']['merged_double']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
        bio_function = join(config['path']['ref'],
                            config['file']['bio_function_completed']),
    output:
        svg_can_targets = '/data/articles/SNORD2_article/svgs/cannonical_targets.svg',
        svg_bio_functions = '/data/articles/SNORD2_article/svgs/bio_functions.svg',
        svg_box_type = '/data/articles/SNORD2_article/svgs/box_type.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/cannonical_targets.py"


rule SNORA12_targets:
    input:
        initial = join(config['path']['raw'], config['file']['initial']),
        parsed_gtf = join(config['path']['ref'], config['file']['light_parsed_gtf']),
        multiple_list = join(config['path']['ref'],
                             config['file']['multiple_gene_name']),
        gene_id = join(config['path']['ref'],
                       config['file']['name_id']),
        sno_host = join(config['path']['ref'],
                        config['file']['sno_host']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/SNORA12_interactions.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/SNORA12_targets.py"

rule piechart_percentage_sno_host_intron:
    input:
        full_sno_host = join(config['path']['sno_host_data'],
                             config['file']['sno_host_data'])
    output:
        svg = '/data/articles/SNORD2_article/svgs/piechart_percentage_sno_host_intron.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/piechart_percentage_sno_host_intron.py"

rule pearson_correlation:
    input:
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm']),
    output:
        svg = '/data/articles/SNORD2_article/svgs/SNORD2_EIF4A2_correlation.svg',
    conda:
        "../envs/python_scikit.yaml"
    script:
        "../scripts/pearson_correlation.py"

rule SNORD2_intron_psi_cov:
    input:
        psi_ext = '/home/danx/Documents/projects/bam_reads_viz/results/merged_psi_ext/EIF4A2/merged_psi_ext.tsv'
    output:
        svg_SNORD2 = '/data/articles/SNORD2_article/svgs/SNORD2_psi_ext.svg',
        svg_SNORD2_intron = '/data/articles/SNORD2_article/svgs/SNORD2_intron_psi_ext.svg',
    conda:
        "../envs/python_scikit.yaml"
    script:
        "../scripts/SNORD2_intron_psi_cov.py"

rule snoRNA_intron_length:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        full_sno_host = join(config['path']['sno_host_data'],
                             config['file']['sno_host_data'])
    output:
        svg = '/data/articles/SNORD2_article/svgs/snoRNA_intron_length.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/snoRNA_intron_length.py"
