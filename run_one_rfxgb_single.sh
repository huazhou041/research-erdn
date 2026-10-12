#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --exclude=cn-0008,cn-0043,cn-0023
set -euo pipefail
f1="$1"; pre="$2"; model="$3"; vc="$4"
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
if [ "$model" = "rf" ]; then
  slurm/deep-drug-response/tune_drug_rf.R recurrence "$f1" slurm/deep-drug-response/hyperparams/grid_rf.json "$pre" -v "$vc" -k 10 --nboot 50
else
  slurm/deep-drug-response/tune_drug_xgb.R recurrence "$f1" slurm/deep-drug-response/hyperparams/grid_xgb.json "$pre" -v "$vc" -k 10 --nboot 50 --max_models 25 --inner_kfold 10 --nthread 4
fi
