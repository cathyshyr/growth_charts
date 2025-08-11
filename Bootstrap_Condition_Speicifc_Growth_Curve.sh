#!/bin/bash
#SBATCH -J Bootstrap
#SBATCH -p batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 02:00:00
#SBATCH --array=1-5000
#SBATCH -o /dev/null
#SBATCH -e /dev/null

module --force purge
module load StdEnv/2023
module load r/4.5.0

cd $SLURM_SUBMIT_DIR
Rscript Bootstrap_Condition_Specific_Growth_Curve.R $SLURM_ARRAY_TASK_ID
