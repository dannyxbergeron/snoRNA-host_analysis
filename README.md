# Workflow for the analysis of snoRNA-host interactions from RNA-RNA interaction datasets

__Author__ : Danny Bergeron

__Email__ :  _<danny.bergeron@usherbrooke.ca>_


## Description
This workflow was developped to investigate and characterize the interactions,
between snoRNAs and their host genes, found in RNA-RNA interaction datasets from
PARIS (PMID: 27180905), LIGR-Seq (PMID: 27184080) and SPLASH (PMID: 27184079).
The first part of the analysis requires to process the raw data with this pipeline:\
<a href="https://github.com/Gabrielle-DF/paris">`https://github.com/Gabrielle-DF/paris`</a>

## Requirements
Conda (Miniconda3) (https://docs.conda.io/en/latest/miniconda.html)
Snakemake (https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Set up
```bash
conda create -c conda-forge -c bioconda -n snakemake snakemake
```

Before running Snakemake, you have to initialize the environment
```bash
conda activate snakemake
```

## Run
To run the workflow locally simply run the following command in the Snakemake conda environment, where `$CORES` is the number of available cores.
```bash
snakemake --use-conda --cores=$CORES
```