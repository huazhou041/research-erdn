#!/bin/bash
set -euo pipefail
sets="${1:-prev,current}"
cohorts="${2:-wholeset,wholeset_coding,stage1,stage1_coding,intersection,intersection_coding,nonintersection,nonintersection_coding}"
models="${3:-rf,xgb}"
mkdir -p logs_rfxgb_smoking res_rfxgb_smoking
IFS=',' read -r -a SL <<< "$sets"; IFS=',' read -r -a CL <<< "$cohorts"; IFS=',' read -r -a ML <<< "$models"
n=0
for s in "${SL[@]}"; do for co in "${CL[@]}"; do for dt in normalized fpkm; do
  dp=data/model_data_smoking/${s}/${co}/${dt}; lab=${s}_${co}_${dt}
  [ ! -d "$dp" ] && { echo "skip missing $dp" >&2; continue; }
  for sub in 12vs0 3vs0 123vs0; do for m in "${ML[@]}"; do
    [ -f "./res_rfxgb_smoking/${lab}_${sub}_${m}_group.done" ] && continue
    sbatch --nice=1000 -J "smk_${m}" -o "logs_rfxgb_smoking/${lab}_${sub}.o_${m}" \
      run_one_rfxgb_smoking.sh "$dp" "$lab" "$m" "$sub"
    n=$((n+1))
  done; done
done; done; done
echo "submitted RF/XGB jobs: $n"
