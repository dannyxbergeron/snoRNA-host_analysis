rule create_sno_host_pair:
    input:
        gene_bed = join(config['path']['ref'],
                        config['file']['gene_bed_biotype']),
        merged_double = join(config['path']['processed'],
                             config['file']['merged_double']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
    output:
        prot_coding_sno_host = join(config['path']['ref'],
                                    config['file']['prot_coding_sno_host']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_sno_host_pair.py"


rule get_sno_intron_loc:
    input:
        prot_coding_sno_host = join(config['path']['ref'],
                                    config['file']['prot_coding_sno_host']),
        parsed = join(config['path']['ref'],
                      config['file']['light_parsed_gtf']),
        tpm = join(config['path']['count_matrices'],
                   config['file']['coco_tpm'])
    output:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        sno_all_transcripts = join(config['path']['ref'],
                                   config['file']['sno_all_transcripts']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/get_sno_intron_loc.py"

# ==========================================================================

rule host_interacting:
    input:
        data_file = join(config['path']['processed'],
                         config['file']['merged_double_inta']),
        parsed_gtf = join(config['path']['ref'],
                          config['file']['light_parsed_gtf']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/host_interaction_analysis.py"


rule create_bedgraph_cons:
    input:
        phastconst = '/data/heavy_stuff/phastconst/phastCons100way.BedGraph',
        sno_all_transcripts = join(config['path']['ref'],
                                   config['file']['sno_all_transcripts']),
    output:
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['host_bedgraph']),
    shell:
        "bedtools intersect -wb -sorted "
        "-a {input.sno_all_transcripts} "
        "-b {input.phastconst} | "
        "awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}{{print$1, $2, $3, $8}}' "
        "| sort -k1,1 -k2,2n -k3,3n "
        "| uniq "
        "> {output.bedgraph}"


rule create_simplified_bedgraph_cons:
    input:
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['host_bedgraph']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        simple_bedgraph = join(config['path']['sno_host_data'],
                               config['file']['simplified_host_bedgraph']),
    shell:
        "awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}"
        "{{if(NR > 1)print \"chr\"$1,$9-100,$10+100,$7}}' {input.sno_host_loc} "
        "| sort -k1,1 -k2,2n -k3,3n > tmp && "
        "bedtools intersect -wb -sorted "
        "-a tmp "
        "-b {input.bedgraph} | "
        "awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}{{print$1, $2, $3, $8}}' "
        "| sort -k1,1 -k2,2n -k3,3n "
        "| uniq "
        "> {output.simple_bedgraph} && "
        "rm tmp"


rule conservation:
    input:
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
        host_ref = join(config['path']['sno_host_data'],
                        config['file']['host_ref']),
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['simplified_host_bedgraph']),
    output:
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/conservation.py"


rule other_conservation:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['simplified_host_bedgraph']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        'tok_'
    params:
        threshold = 0.5
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/other_conservation.py"


rule other_conservation_bootstrap:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        bedgraph = join(config['path']['sno_host_data'],
                        config['file']['simplified_host_bedgraph']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        'tok'
    params:
        threshold = 0.5
    threads:
        6
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/other_conservation_bootstrap.py"


rule transcript_per_gene:
    input:
        parsed = join(config['path']['ref'],
                      config['file']['parsed_gtf']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        sno_host = join(config['path']['sno_host_data'],
                        config['file']['sno_host_data']),
    output:
        # svg = '/data/labmeetings/host_interactions/transcript_per_gene.svg'
        svg = '/data/articles/SNORD2_article/svgs/transcript_per_gene.svg',
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/transcript_per_gene.py"

rule alternative_splicing_intron:
    input:
        parsed = join(config['path']['ref'],
                      config['file']['light_parsed_gtf']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        cons = join(config['path']['sno_host_data'],
                    config['file']['cons'])
    output:
        alt_splice = join(config['path']['sno_host_data'],
                          config['file']['alt_splice']),
        # graph1 = '/data/labmeetings/host_interactions/alt_splicing_bar.svg',
        graph2 = '/data/articles/SNORD2_article/svgs/alternative_splicing_intron.svg',
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/alternative_splicing_intron.py"


rule prepare_bedgraph_search:
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        introns = join(config['path']['tmp'],
                       config['file']['introns_sno_in_host']),
    shell:
        "colTab -f {input.sno_host_loc} -c chr,intron_start,intron_end "
        "| awk '{{if(NR > 1){{print \"chr\" $0}}}}' "
        "| sort -k1,1 -k2,2n -k3,3 "
        "| uniq "
        "| awk 'BEGIN{{FS=\"\t\";OFS=\"\t\"}}{{print$1, $2-100, $3+100}}' "
        "> {output.introns}"

rule get_intron_bed:
    input:
        bed = join(config['path']['beds'], 'sorted_clean_{bed}.bedgraph'),
        introns = join(config['path']['tmp'],
                       config['file']['introns_sno_in_host']),
    output:
        bg = join(config['path']['tissues_cell_bg'],
                  'intron_{bed}.bedgraph'),
    conda:
        "../envs/python.yaml"
    shell:
        "set +o pipefail; cat {input.bed} | scripts/getIntronBG {input.introns} "
        "> {output.bg}"


rule reads_in_extensions:
    input:
        alt_splice = join(config['path']['sno_host_data'],
                          config['file']['alt_splice']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
        bg = expand(join(config['path']['tissues_cell_bg'], 'intron_{bed}.bedgraph'),
                    bed=config['bedgraphs'])
    output:
        bed_viz = join(config['path']['tmp'],
                       config['file']['bedgraph_viz']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/reads_in_extensions.py"

# ---------------- TEST FOR NOW --------------------------------
rule host_int_multiple_sno_in_host:
    input:
        merged_file = join(config['path']['processed'],
                           config['file']['merged_double_inta'])
    output:
        'tok'
    conda:
        "../envs/python_seabord.yaml"
    script:
        "../scripts/host_int_multiple_sno_in_host.py"
