#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --job-name="pps"
#SBATCH -o cs-%j.out
#SBATCH -e cs-%j.err 
#SBATCH --mail-user=youremail@email.org
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem-per-cpu=4G


cd $SLURM_SUBMIT_DIR

module load singularity

singularity run http://s3-far.jax.org/builder/builder sing_img.def sing.sif















	







