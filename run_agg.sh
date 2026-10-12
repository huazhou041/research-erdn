#!/bin/bash
#SBATCH --job-name=agg
#SBATCH --partition=cpu_short
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --output=agg_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
Rscript aggregate_fsmethods.R
echo AGG_GLM_DONE
Rscript aggregate_rfxgb.R
echo AGG_RFXGB_DONE
