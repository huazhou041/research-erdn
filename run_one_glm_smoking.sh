#!/bin/bash
#SBATCH --partition=cpu_medium,cpu_long
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --exclude=cn-0008,cn-0043
set -euo pipefail
datapath="${1}"; label="${2}"
fs_method="${3:?fs_method (mad|cor|mad_nosexmt)}"; subset="${4:?subset (12vs0|3vs0|123vs0)}"
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
[ ! -d "res_smoking" ] && mkdir -p res_smoking
extra_args=()
if [ "$fs_method" = "mad_nosexmt" ]; then extra_args+=(--exclude_sexmt_annot raw_data/gene_type.csv); fi
for varcut in 50 100 200 500 1000; do
  for dataname in clinics fpkm fpkm_clinics module20 module20_clinics; do
    data=${dataname}_${subset}
    f1=./${datapath}/${data}.RDS
    [ ! -f "$f1" ] && continue
    out=./res_smoking/${label}_${data}_glm_${fs_method}_v${varcut}_cvSummary.csv
    [ ! -f "$out" ] && slurm/deep-drug-response/tune_drug.R recurrence "$f1" 'slurm/deep-drug-response/hyperparams/grid_glm.json' ./res_smoking/${label}_${data}_glm_${fs_method}_v${varcut} -v "$varcut" --learner 'GLM' -k 10 --nboot 100 --fs_method "$fs_method" "${extra_args[@]}"
  done
done
touch ./res_smoking/${label}_${subset}_glm_${fs_method}_group.done
