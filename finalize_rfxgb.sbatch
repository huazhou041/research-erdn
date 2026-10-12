#!/bin/bash
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -J finalize_rfxgb
#SBATCH -o logs_smoking/finalize_rfxgb.o
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
mkdir -p logs_smoking
# All 10 cohorts are submitted by submit_rfxgb.sh, so aggregate_rfxgb.R produces the
# COMPLETE rf/xgb summary (rf_xgb_auc_202607.csv) directly -- no merge with an old CSV.
echo "=== aggregate RF/XGB (all cohorts incl novel) -> rf_xgb_auc_202607.csv ==="
Rscript aggregate_rfxgb.R
echo "=== FINALIZE_DONE ==="
