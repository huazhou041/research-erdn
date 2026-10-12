#!/bin/bash
#SBATCH --job-name=rocni
#SBATCH --partition=cpu_dev
#SBATCH --time=03:50:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --output=rocni_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
Rscript make_roc_glm_nonint.R
Rscript make_roc_rfxgb_nonint.R
echo ALL_ROC_DONE
