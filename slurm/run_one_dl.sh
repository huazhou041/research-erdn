#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --exclude=cn-0008,cn-0043
set -euo pipefail
datapath="${1}"; label="${2}"; subset="${3}"
set +u
module load miniconda3/4.6.14
source activate r_env
unset JAVA_HOME          # profile sets it to the literal string ":-", which breaks h2o.init
set -u
mkdir -p res_dl
for varcut in 50 100 200 500 1000; do
  for dataname in clinics module20 module20_clinics fpkm fpkm_clinics; do
    f1=./${datapath}/${dataname}_${subset}.RDS
    [ ! -f "$f1" ] && continue
    out=./res_dl/${label}_${dataname}_${subset}_dl_v${varcut}_summary.csv
    [ -f "$out" ] && continue
    slurm/deep-drug-response/tune_drug_dl.R recurrence "$f1" \
      slurm/deep-drug-response/hyperparams/grid_dl.json \
      ./res_dl/${label}_${dataname}_${subset}_dl_v${varcut} \
      -v "$varcut" -k 10 --nboot 10 --max_models 8 --inner_nfolds 3 --max_mem 12G --cores 2
  done
done
touch ./res_dl/${label}_${subset}_dl_group.done
