#!/bin/bash
#SBATCH --job-name=perm3
#SBATCH --partition=cpu_dev
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=perm3_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
Rscript perm3vs0.R
