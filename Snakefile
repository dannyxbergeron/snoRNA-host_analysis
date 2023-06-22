from os.path import join

configfile: "config.json"

rule all:
    input:
        full_merge = join(config['path']['processed'],
                          config['file']['full_network_clip']),
        interactions = join(config['path']['cytoscape'],
                            config['network_file']['interactions']),
        edges = join(config['path']['cytoscape'],
                     config['network_file']['edges']),
        beds = expand(join(config['path']['beds'], 'sorted_clean_{bed}.bedgraph'),
                      bed=config['bedgraphs']),
        intron_tok = expand(join(config['path']['tissues_cell_bg'],
                                        'intron_{bed}.bedgraph'),
                            bed=config['bedgraphs']),
        ext_ratio = join(config['path']['sno_host_data'],
                         config['file']['ext_ratio']),
        dist_num_interactions = '/data/articles/SNORD2_article/svgs/dist_num_interactions.svg',
        intra_sno_hosts = '/data/articles/SNORD2_article/svgs/distance_target_from_snoRNA.svg',
        upstream_vs_downstream = '/data/articles/SNORD2_article/svgs/upstream_vs_downstream.svg',
        cannonical_targets = '/data/articles/SNORD2_article/svgs/cannonical_targets.svg',
        bp_distance = join(config['path']['branchpoint'], config['file']['branchpoint']),
        best_branch_points = join(config['path']['branchpoint'],
                                  config['file']['best_branchpoint']),
        svg_branch_point = '/data/articles/SNORD2_article/svgs/branchpoint_binding.svg',
        sno_intron_coordinates = join(config['path']['sno_intron'],
                                      config['file']['sno_intron_coord']),
        merged_mfe = join(config['path']['sno_intron'], config['file']['merged_mfe']),
        mfe_svg = '/data/articles/SNORD2_article/svgs/mfe_plot.svg',
        svg = '/data/articles/SNORD2_article/svgs/SNORA12_interactions.svg',
        alternative_splicing_intron = '/data/articles/SNORD2_article/svgs/alternative_splicing_intron.svg',
        other_conservation = '/data/articles/SNORD2_article/svgs/other_conservation.svg',
        piechart_percentage_sno_host_intron = '/data/articles/SNORD2_article/svgs/piechart_percentage_sno_host_intron.svg',
        coords = join(config['path']['intaRNA_simulation'], config['file']['intaRNA_simul_coord']),
        intaRNA_prediction = '/data/articles/SNORD2_article/svgs/intaRNA_simulation_distribution.svg',
        snord2_eif4a2_corr_svg = '/data/articles/SNORD2_article/svgs/SNORD2_EIF4A2_correlation.svg',
        svg_SNORD2_intron = '/data/articles/SNORD2_article/svgs/SNORD2_intron_psi_ext.svg',
        intron_length = '/data/articles/SNORD2_article/svgs/snoRNA_intron_length.svg',




# Include raw processing
include: "rules/process_raw.smk"

# Include building of the network
include: "rules/network.smk"

# Include files
include: "rules/process_ref.smk"
include: "rules/host_interactions.smk"
include: "rules/graphs_and_analysis.smk"
include: "rules/branch_point_prediction.smk"
include: "rules/mfe_intron_folding.smk"
include: "rules/intaRNA_simulation.smk"

# Include rule for sno_host classification ?
include: "rules/classification.smk"
