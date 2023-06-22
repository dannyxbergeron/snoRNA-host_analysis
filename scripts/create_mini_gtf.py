import numpy as np
import pandas as pd

gtf = snakemake.input.gtf
sno_host_loc = snakemake.input.sno_host_loc

out_file = snakemake.output.mini_gtf

def main():

    df = pd.read_csv(sno_host_loc, sep='\t')
    host_ids = list(df.host_id)
    transcript_ids = list(df.host_transcript_id)

    to_keep = []
    with open(gtf, 'r') as f:
        with open(out_file, 'w') as out:
            for line in f.read().splitlines():
                for host_id, transcript_id in zip(host_ids, transcript_ids):
                    if host_id in line:
                        if '\tgene\t' in line \
                                or transcript_id in line:
                            if line not in to_keep:
                                to_keep.append(line)
                                out.write(line)
                                out.write('\n')
                            continue



if __name__ == '__main__':
    main()
