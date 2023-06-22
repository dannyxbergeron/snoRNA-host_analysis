#https://bioconductor.org/packages/release/bioc/manuals/branchpointer/man/branchpointer.pdf
#https://www.bioconductor.org/packages/release/bioc/vignettes/branchpointer/inst/doc/branchpointer.pdf

rule generate_intron_file:
    """
    query regions must contain a branchpoint window - that is the region located at -18 to -44
    from the 3’ splice site.
    """
    input:
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        intron_file = join(config['path']['ref'],
                           config['file']['intron_file'])
    shell:
        """
        cat {input.sno_host_loc} \
        | colTab -f IN -c gene_id,chr,intron_start,intron_end,strand \
        | sed -e 's:gene_id:id:g' -e 's:chr:chromosome:g' \
        | sed -e 's:intron_start:start:g' -e 's:intron_end:end:g' \
        | awk 'BEGIN{{OFS="\t";FS="\t"}}\
        {{if(NR > 1){{if($5=="+"){{$3=$4-44;$4=$4-18}}\
        else{{$4=$3+44;$3=$3+18}}}}\
        print$0}}' \
        > {output.intron_file}
        """

rule generate_smaller_gtf:
    """
    Generate a smaller gtf since it does not currently work with the full gtf...
    """
    input:
        gtf = join(config['path']['ref'],
                   config['file']['gtf']),
        sno_host_loc = join(config['path']['ref'],
                            config['file']['prot_cod_sno_host_loc']),
    output:
        mini_gtf = join(config['path']['branchpoint'],
                        config['file']['mini_gtf']),
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/create_mini_gtf.py"

rule branchpointer:
    """ Find the branchpoint in all introns of snoRNAs """
    input:
        mini_gtf = join(config['path']['branchpoint'],
                        config['file']['mini_gtf']),
        intron_file = join(config['path']['ref'],
                           config['file']['intron_file']),
        genome = "/data/heavy_stuff/genome_fasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        bp_distance = join(config['path']['branchpoint'],
                           config['file']['branchpoint']),
    params:
        bedtools = "/usr/bin/bedtools",
        tmp_dir = config['path']['tmp']
    conda:
        "../envs/r_branchpointer.yaml"
    script:
        "../scripts/branchpoint_prediction.R"

rule keep_best_branchpoints:
    """
    Extract and keep the best branchpoint position
    for the branchpointer prediction
    """
    input:
        bp_distance = join(config['path']['branchpoint'],
                           config['file']['branchpoint']),
        intra = join(config['path']['sno_host_data'],
                     config['file']['alt_splice']),
    output:
        best_branch_points = join(config['path']['branchpoint'],
                                  config['file']['best_branchpoint']),
        intra_with_bp = join(config['path']['branchpoint'],
                             config['file']['alt_splice_with_bp']),
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/keep_best_branchpoints.py"

rule branchpoint_binding:
    """ Plot a piechart with proportion of snoRNA binding the branchpoint """
    input:
        intra_with_bp = join(config['path']['branchpoint'],
                             config['file']['alt_splice_with_bp']),
    output:
        svg_branch_point = '/data/articles/SNORD2_article/svgs/branchpoint_binding.svg',
    conda:
        "../envs/python_seaborn.yaml"
    script:
        "../scripts/branchpoint_binding.py"
