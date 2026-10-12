#!/bin/bash
#SBATCH --job-name=rocrf
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=rocrf_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
Rscript make_roc_rfxgb_nonint.R
