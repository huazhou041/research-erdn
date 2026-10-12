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
mkdir -p res_rfxgb
# RF/XGB use MAD feature selection only (no --fs_method flag), matching the GLM `mad` arm.
# XGB does nested inner-CV hyperparameter selection (--inner_kfold), the fix behind
# archive_xgb_biased_20260706. Both already carry the direction="<" AUC fix.
# NOTE: a single degenerate cell (e.g. novel/3vs0 fpkm_clinics, where too few events
# make XGB's select_best_params cut() fail with "invalid number of intervals") must NOT
# abort the whole group -- otherwise every later varcut/dataname is lost and group.done
# is never written. Such cells are logged to *_failed.log and skipped; the group still
# completes so its computable cells and group.done are produced.
faillog=./res_rfxgb/${label}_${subset}_${model}_failed.log
for varcut in 50 100 200 500 1000; do
  for dataname in clinics module20 module20_clinics fpkm fpkm_clinics; do
    f1=./${datapath}/${dataname}_${subset}.RDS
    [ ! -f "$f1" ] && continue
    out=./res_rfxgb/${label}_${dataname}_${subset}_${model}_v${varcut}_cvSummary.csv
    [ -f "$out" ] && continue
    pre=./res_rfxgb/${label}_${dataname}_${subset}_${model}_v${varcut}
    if [ "$model" = "rf" ]; then
      slurm/deep-drug-response/tune_drug_rf.R recurrence "$f1" \
        slurm/deep-drug-response/hyperparams/grid_rf.json "$pre" \
        -v "$varcut" -k 10 --nboot 50 \
        || { echo "$(date '+%F %T') CELL_FAILED $out" >> "$faillog"; }
    else
      slurm/deep-drug-response/tune_drug_xgb.R recurrence "$f1" \
        slurm/deep-drug-response/hyperparams/grid_xgb.json "$pre" \
        -v "$varcut" -k 10 --nboot 50 --max_models 25 --inner_kfold 10 --nthread 4 \
        || { echo "$(date '+%F %T') CELL_FAILED $out" >> "$faillog"; }
    fi
  done
done
touch ./res_rfxgb/${label}_${subset}_${model}_group.done
