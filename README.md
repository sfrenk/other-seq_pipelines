# Pipelines for performing NGS analysis in C elegans

See rna-seq_pipelines repository for RNA-seq specific pipelines

These pipelines have been written for use with the SLURM scheduling system, but they can be used on any system with a few slight modifications.

## Script types

### Snakemake

I would recommend using the Snakemake scripts for running these pipelines. To run a pipeline, use setup_dir.sh to copy in the appropriate snakefile and create a bash script for running the pipeline. NOTE: setup_dir.sh can be found in the rna-seq_pipelines repository

setup_dir.sh -d (directory containing samples) -p (pipeline name)

### Bash scripts

Use these as a simple alternative to Snakemake (warning - these pipelines may not have been updated recently).
