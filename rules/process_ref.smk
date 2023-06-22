rule parse_and_process_ref:
    input:
        gtf = join(config['path']['ref'],
                   config['file']['gtf']),
    output:
        parsed = join(config['path']['ref'],
                      config['file']['parsed_gtf']),
        gene_bed_biotype = join(config['path']['ref'],
                                config['file']['gene_bed_biotype']),
    shell:
        "gtfParser find_parse {input.gtf} > {output.parsed} && "
        "sed 's/seqname/chr/g' {output.parsed} "
        "| awk '$3 == \"gene\" || $1 == \"chr\"' "
        "| colTab -f IN -c chr,start,end,gene_id,gene_name,strand,gene_biotype "
        "> {output.gene_bed_biotype}"
