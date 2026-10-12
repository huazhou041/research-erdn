#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -J smk_post
#SBATCH -o logs_smoking/post_smoking.o
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
mkdir -p logs_smoking
echo "=== [1/5] aggregate GLM ==="; Rscript aggregate_fsmethods_smoking.R
echo "=== [2/5] aggregate RF/XGB ==="; Rscript aggregate_rfxgb_smoking.R
echo "=== [3/5] GLM ROC ==="; Rscript make_roc_curves_fsmethods_smoking.R
echo "=== [4/5] RF/XGB ROC ==="; Rscript make_roc_rfxgb_smoking.R
echo "=== [5/5] build dashboard ==="; Rscript codes_20260713/make_dashboard_fsmethods_smoking.R
echo "=== POST_SMOKING_DONE ==="
