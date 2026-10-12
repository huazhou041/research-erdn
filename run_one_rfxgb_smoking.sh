#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=24G
#SBATCH --time=48:00:00
#SBATCH --exclude=cn-0008,cn-0043
set -euo pipefail
datapath="${1}"; label="${2}"; model="${3}"; subset="${4}"
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
mkdir -p res_rfxgb_smoking
for varcut in 50 100 200 500 1000; do
  for dataname in clinics module20 module20_clinics fpkm fpkm_clinics; do
    f1=./${datapath}/${dataname}_${subset}.RDS
    [ ! -f "$f1" ] && continue
    out=./res_rfxgb_smoking/${label}_${dataname}_${subset}_${model}_v${varcut}_cvSummary.csv
    [ -f "$out" ] && continue
    pre=./res_rfxgb_smoking/${label}_${dataname}_${subset}_${model}_v${varcut}
    if [ "$model" = "rf" ]; then
      slurm/deep-drug-response/tune_drug_rf.R recurrence "$f1" slurm/deep-drug-response/hyperparams/grid_rf.json "$pre" -v "$varcut" -k 10 --nboot 50
    else
      slurm/deep-drug-response/tune_drug_xgb.R recurrence "$f1" slurm/deep-drug-response/hyperparams/grid_xgb.json "$pre" -v "$varcut" -k 10 --nboot 50 --max_models 25 --inner_kfold 10 --nthread 4
    fi
  done
done
touch ./res_rfxgb_smoking/${label}_${subset}_${model}_group.done
