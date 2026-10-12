#!/bin/bash
# GLM-only comparison of two gene-prefilter methods (mad = original unsupervised
# variance, cor = supervised |correlation-with-outcome|, both fold-local) across
# the erdn_202607 prev/current model_data layout.
#
# Usage:
#   ./run_submit_20260709_glm_fsmethods.sh [sets] [cohorts] [fs_methods]
#   sets       : comma list (default: prev,current)
#   cohorts    : comma list (default: all 6 cohorts)
#   fs_methods : comma list (default: mad,cor)

set -euo pipefail

sets="${1:-prev,current}"
cohorts="${2:-wholeset,wholeset_coding,stage1,stage1_coding,intersection,intersection_coding,nonintersection,nonintersection_coding,novel,novel_coding}"
fs_methods="${3:-mad,cor,mad_nosexmt}"

curname='erdn20260709_glmfs'
[ ! -d "logs" ] && mkdir -p logs
[ ! -d "res" ] && mkdir -p res

IFS=',' read -r -a set_list <<< "$sets"
IFS=',' read -r -a cohort_list <<< "$cohorts"
IFS=',' read -r -a fsm_list <<< "$fs_methods"

for set_name in "${set_list[@]}"; do
  for cohort in "${cohort_list[@]}"; do
    for datatype in normalized fpkm; do
      datapath=data/model_data/${set_name}/${cohort}/${datatype}
      label=${set_name}_${cohort}_${datatype}

      if [ ! -d "$datapath" ]; then
        echo "skip (missing): $datapath" >&2
        continue
      fi

      for subset in 12vs0 3vs0 123vs0; do
        for fsm in "${fsm_list[@]}"; do
          [ -f "./res/${label}_${subset}_glm_${fsm}_group.done" ] && continue
          sbatch --nice=1000 -J "${curname}_${fsm}" \
            -o "logs/${label}_${subset}_${fsm}.o_glm" \
            slurm/run_one_dataset_type_model.sh "$datapath" "$label" "$fsm" "$subset"
        done
      done
    done
  done
done
