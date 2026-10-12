#!/bin/bash
set -euo pipefail
sets="${1:-prev,current}"
cohorts="${2:-wholeset,wholeset_coding,stage1,stage1_coding,intersection,intersection_coding,nonintersection,nonintersection_coding,novel,novel_coding}"
models="${3:-rf,xgb}"
mkdir -p logs_rfxgb res_rfxgb
IFS=',' read -r -a SL <<< "$sets"; IFS=',' read -r -a CL <<< "$cohorts"; IFS=',' read -r -a ML <<< "$models"
for s in "${SL[@]}"; do for co in "${CL[@]}"; do for dt in normalized fpkm; do
  dp=data/model_data/${s}/${co}/${dt}; lab=${s}_${co}_${dt}
  [ ! -d "$dp" ] && continue
  for sub in 12vs0 3vs0 123vs0; do for m in "${ML[@]}"; do
    [ -f "./res_rfxgb/${lab}_${sub}_${m}_group.done" ] && continue
    sbatch --nice=1000 -J "erdn_${m}" -o "logs_rfxgb/${lab}_${sub}.o_${m}" \
      slurm/run_one_rfxgb.sh "$dp" "$lab" "$m" "$sub"
  done; done
done; done; done
