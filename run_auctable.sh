#!/bin/bash
#SBATCH --job-name=auctable
#SBATCH --partition=cpu_short
#SBATCH --time=03:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --output=auctable_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
Rscript -e "cat(\"pkgs:\", all(c(\"glmnet\",\"xgboost\",\"randomForest\",\"pROC\",\"caret\") %in% rownames(installed.packages())),\"\n\")"
Rscript auc_table_nonint.R
