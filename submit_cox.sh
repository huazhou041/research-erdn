#!/bin/bash
set -euo pipefail
mkdir -p logs_cox res_cox
for set_name in prev current; do
  for cohort in wholeset wholeset_coding stage1 stage1_coding intersection intersection_coding; do
    for datatype in normalized fpkm; do
      datapath=data/model_data/${set_name}/${cohort}/${datatype}
      label=${set_name}_${cohort}_${datatype}
      [ ! -d "$datapath" ] && continue
      for endpoint in recur all; do
        [ -f "./res_cox/${label}_cox_${endpoint}_mad_group.done" ] && continue
        sbatch --nice=1000 -J "erdn_cox_${endpoint}" \
          -o "logs_cox/${label}_${endpoint}.o_cox" \
          slurm/run_one_cox.sh "$datapath" "$label" "$endpoint" mad
      done
    done
  done
done
