#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --job-name="pps"
#SBATCH -o ppr-%j.out
#SBATCH -e ppr-%j.err 
#SBATCH --mail-user=youremail@email.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem-per-cpu=8G
#SBATCH --array=3,7,9,38,50,51,52,53,54,55,56,58,72,74,76,77

cd $SLURM_SUBMIT_DIR/fastq_$SLURM_ARRAY_TASK_ID

module load singularity

cp ../Snakefile .

singularity exec ../sing.sif bash ../run_snakemake.sh ../case_control_c${SLURM_ARRAY_TASK_ID}.tsv /Snakemake/fastq_$SLURM_ARRAY_TASK_ID
















	







