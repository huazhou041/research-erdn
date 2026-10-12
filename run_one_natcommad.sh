#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --exclude=cn-0008,cn-0043,cn-0023
set -euo pipefail
datapath="${1}"; label="${2}"; subset="${3}"
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
[ ! -d res ] && mkdir -p res
RANK="${ERDN_RANK_GENES:-raw_data/natcom_mad_ranked.csv}"
# natcommad = FIXED external nature_com MAD gene list; only gene-based feature sets apply.
for varcut in 50 100 200 500 1000; do
  for dataname in fpkm fpkm_clinics; do
    data=${dataname}_${subset}
    f1=./${datapath}/${data}.RDS
    [ ! -f "$f1" ] && continue
    out=./res/${label}_${data}_glm_natcommad_v${varcut}_cvSummary.csv
    [ ! -f "$out" ] && slurm/deep-drug-response/tune_drug.R recurrence "$f1" 'slurm/deep-drug-response/hyperparams/grid_glm.json' ./res/${label}_${data}_glm_natcommad_v${varcut} -v "$varcut" --learner 'GLM' -k 10 --nboot 100 --fs_method natcommad --rank_genes_file "$RANK"
  done
done
touch ./res/${label}_${subset}_glm_natcommad_group.done
