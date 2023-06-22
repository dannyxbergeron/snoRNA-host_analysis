import sys
import os
import pandas as pd

from snakemake.shell import shell

input_file = snakemake.input.coords
chr_dir = snakemake.params.chr_dir
rev_complement = snakemake.params.rev_complement_script

output_dir = snakemake.output.fasta_dir


coords = pd.read_csv(input_file, sep='\t')
print(coords.columns)
print(coords)

def extract(coord):
    chr, start_end, strand = coord.split(':')
    start, end = start_end.split('-')
    return chr.replace('chr', ''), start, end, strand

# Create the output folder if it doesn't exists
shell(f'mkdir -p {output_dir}')

for DG, snoRNA, target, simulated in coords.values:
    for type, coord in zip(coords.columns[1:], [snoRNA, target, simulated]):
        chr, start, end, strand = extract(coord)

        # gtf to bed conversion
        start = int(start) - 1

        # To reverse the fasta if needed
        rev = '_reverse' if strand == '-' else ''

        chrom_file = os.path.join(chr_dir, f'{chr}.2bit')
        name = f'{DG}-{type}{rev}.fa'
        cmd = 'twoBitToFa {}:{}:{}-{} {}/{}'.format(chrom_file,
                                                    chr,
                                                    start,
                                                    end,
                                                    output_dir,
                                                    name)
        shell(cmd)

        if rev == '_reverse':
            new_name = f'{DG}-{type}.fa'
            new_cmd = f'python {rev_complement} {output_dir}/{name} > {output_dir}/{new_name}'
            shell(new_cmd)
            shell('rm {output_dir}/{name}')
