#!/bin/bash
#SBATCH --job-name=agg_all
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH -o agg_all_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
echo "=== aggregate_fsmethods (full GLM+natcommad) ==="; Rscript aggregate_fsmethods.R && echo OK_FSMETHODS
echo "=== aggregate_rfxgb (full) ==="; Rscript aggregate_rfxgb.R && echo OK_RFXGB
echo "=== aggregate_fsmethods_smoking ==="; Rscript aggregate_fsmethods_smoking.R && echo OK_FSMETHODS_SMK
echo "=== aggregate_rfxgb_smoking ==="; Rscript aggregate_rfxgb_smoking.R && echo OK_RFXGB_SMK
echo "=== ALL_AGG_DONE ==="
