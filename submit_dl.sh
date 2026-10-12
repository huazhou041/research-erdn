#!/bin/bash
set -euo pipefail
mkdir -p logs_dl res_dl
for set_name in prev current; do
  for cohort in wholeset stage1 intersection; do
    for datatype in normalized fpkm; do
      datapath=data/model_data/${set_name}/${cohort}/${datatype}
      label=${set_name}_${cohort}_${datatype}
      [ ! -d "$datapath" ] && continue
      for subset in 12vs0 123vs0; do
        [ -f "./res_dl/${label}_${subset}_dl_group.done" ] && continue
        sbatch --nice=1000 -J "erdn_dl" \
          -o "logs_dl/${label}_${subset}.o_dl" \
          slurm/run_one_dl.sh "$datapath" "$label" "$subset"
      done
    done
  done
done
