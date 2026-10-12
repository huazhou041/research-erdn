#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --exclude=cn-0008,cn-0043
set -euo pipefail
datapath="${1}"; label="${2}"; endpoint="${3}"; fsm="${4}"
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
mkdir -p res_cox
# Cox uses the FULL cohort (123vs0 RDS = all patients); the endpoint decides which
# events count. endpoint=recur censors second primaries instead of dropping them.
for varcut in 50 100 200 500 1000; do
  for dataname in clinics fpkm fpkm_clinics module20 module20_clinics; do
    f1=./${datapath}/${dataname}_123vs0.RDS
    [ ! -f "$f1" ] && continue
    out=./res_cox/${label}_${dataname}_cox_${endpoint}_${fsm}_v${varcut}_summary.csv
    [ -f "$out" ] && continue
    slurm/deep-drug-response/tune_drug_cox.R "$f1" surv_table_202606.csv \
      ./res_cox/${label}_${dataname}_cox_${endpoint}_${fsm}_v${varcut} \
      --endpoint "$endpoint" --varcut "$varcut" --kfold 10 --nboot 50 --fs_method "$fsm"
  done
done
touch ./res_cox/${label}_cox_${endpoint}_${fsm}_group.done
