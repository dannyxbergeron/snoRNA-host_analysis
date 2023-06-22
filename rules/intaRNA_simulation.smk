rule create_intaRNA_simulated_windows:
    input:
        intra = join(config['path']['sno_host_data'],
                     config['file']['alt_splice'])
    output:
        coords = join(config['path']['intaRNA_simulation'],
                      config['file']['intaRNA_simul_coord'])
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/create_intaRNA_simulated_windows.py"


rule generate_fasta_from_coords:
    input:
        coords = join(config['path']['intaRNA_simulation'],
                      config['file']['intaRNA_simul_coord'])
    output:
        fasta_dir = directory(config['path']['intaRNA_simulation_fasta'])
    params:
        chr_dir = '/home/danx/Documents/projects/fold_and_viz_RNA/workflow/resources/chrom/2bit',
        rev_complement_script = '/home/danx/Documents/projects/fold_and_viz_RNA/workflow/scripts/reverse_complement.py'
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/generate_fasta_from_coords.py"

rule intaRNA_prediction:
    input:
        coords = join(config['path']['intaRNA_simulation'],
                      config['file']['intaRNA_simul_coord'])
    output:
        svg = '/data/articles/SNORD2_article/svgs/intaRNA_simulation_distribution.svg',
        cumsum_svg = '/data/articles/SNORD2_article/svgs/intaRNA_simulation_cumsum.svg',
    params:
        fasta_dir = config['path']['intaRNA_simulation_fasta']
    conda:
        "../envs/intaRNA.yaml"
    script:
        "../scripts/intaRNA_prediction.py"
