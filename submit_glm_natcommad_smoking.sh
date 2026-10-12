#!/bin/bash
# GLM with fs_method=natcommad on the EVER-SMOKER subset: FIXED nature_com top-MAD
# gene list transferred to the smoking-subset prev/current data. Gene feature sets
# (fpkm, fpkm_clinics) only. 8 cohorts (no novel, matching the smoking build).
# Requires raw_data/natcom_mad_ranked.csv (build_natcom_madgenes.R).
set -euo pipefail
sets="${1:-prev,current}"
cohorts="${2:-wholeset,wholeset_coding,stage1,stage1_coding,intersection,intersection_coding,nonintersection,nonintersection_coding}"
mkdir -p logs_smoking res_smoking
[ -f raw_data/natcom_mad_ranked.csv ] || { echo "ERROR: raw_data/natcom_mad_ranked.csv missing -- run build_natcom_madgenes.R first" >&2; exit 1; }
IFS=',' read -r -a SL <<< "$sets"; IFS=',' read -r -a CL <<< "$cohorts"
n=0
for s in "${SL[@]}"; do for co in "${CL[@]}"; do for dt in normalized fpkm; do
  dp=data/model_data_smoking/${s}/${co}/${dt}; lab=${s}_${co}_${dt}
  [ ! -d "$dp" ] && continue
  for sub in 12vs0 3vs0 123vs0; do
    [ -f "./res_smoking/${lab}_${sub}_glm_natcommad_group.done" ] && continue
    sbatch --nice=1000 -J smk_natcommad -o "logs_smoking/${lab}_${sub}_natcommad.o_glm" run_one_natcommad_smoking.sh "$dp" "$lab" "$sub"
    n=$((n+1))
  done
done; done; done
echo "submitted smoking natcommad GLM jobs: $n"
