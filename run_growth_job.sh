#!/bin/bash
#SBATCH -J growth
#SBATCH -p batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -t 5-00:00:00
#SBATCH --array=1-43
#SBATCH -o logs/growth_%A_%a.out
#SBATCH -e logs/growth_%A_%a.err

module --force purge
module load StdEnv/2023
module load r/4.5.0

cd $SLURM_SUBMIT_DIR

mkdir -p logs

Rscript run_growth_job.R ${SLURM_ARRAY_TASK_ID}