#!/bin/sh
source activate rsem-pipeline

snakemake --resources mem_mb=8000 --until run_hbadeals --cores 32  --config samples_table=$1 fastq_dir=$2
