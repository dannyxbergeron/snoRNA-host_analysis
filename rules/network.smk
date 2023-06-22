subworkflow clip_analysis:
    workdir:
        "../Tommy_stuff/dan_new_analysis"

rule process_clip:
    input:
        merge_clip = clip_analysis("data_processed/all_merged.bed")
    output:
        merge_clip = join(config['path']['clip'],
                          config['file']['clip_data'])
    shell:
        "cp {input.merge_clip} {output.merge_clip}"


rule merge_network_and_clip:
    input:
        merged_double_inta = join(config['path']['processed'],
                                  config['file']['merged_double_inta']),
        merged_clip = join(config['path']['clip'],
                           config['file']['clip_data'])
    output:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip'])
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/network_analysis.py"


rule build_network:
    input:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
        snodb = join(config['path']['ref'],
                     config['file']['snoDB']),
    output:
        interactions = join(config['path']['cytoscape'],
                            config['network_file']['interactions']),
        nodes = join(config['path']['cytoscape'],
                     config['network_file']['nodes']),
        edges = join(config['path']['cytoscape'],
                     config['network_file']['edges']),
        mapped_clip_only = join(config['path']['processed'],
                                config['file']['mapped_clip_only']),
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_network_V2_with_ENCODE_data.py"
