#!/bin/bash
set -euo pipefail
sets="${1:-prev,current}"
cohorts="${2:-wholeset,wholeset_coding,stage1,stage1_coding,intersection,intersection_coding,nonintersection,nonintersection_coding}"
fs_methods="${3:-mad,cor,mad_nosexmt}"
mkdir -p logs_smoking res_smoking
IFS=',' read -r -a SL <<< "$sets"; IFS=',' read -r -a CL <<< "$cohorts"; IFS=',' read -r -a FL <<< "$fs_methods"
n=0
for s in "${SL[@]}"; do for co in "${CL[@]}"; do for dt in normalized fpkm; do
  dp=data/model_data_smoking/${s}/${co}/${dt}; lab=${s}_${co}_${dt}
  [ ! -d "$dp" ] && { echo "skip missing $dp" >&2; continue; }
  for sub in 12vs0 3vs0 123vs0; do for fsm in "${FL[@]}"; do
    [ -f "./res_smoking/${lab}_${sub}_glm_${fsm}_group.done" ] && continue
    sbatch --nice=1000 -J "smk_glm_${fsm}" -o "logs_smoking/${lab}_${sub}_${fsm}.o_glm" \
      run_one_glm_smoking.sh "$dp" "$lab" "$fsm" "$sub"
    n=$((n+1))
  done; done
done; done; done
echo "submitted GLM jobs: $n"
