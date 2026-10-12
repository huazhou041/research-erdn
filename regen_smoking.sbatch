#!/bin/bash
#SBATCH -p cpu_short,cpu_medium
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH -c 4
#SBATCH -J regen_smoking
#SBATCH -o /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609/logs_smoking/regen_smoking.o
set -euo pipefail
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
mkdir -p logs_smoking
# preserve the Aug-10 (row-subset) data for a value-equality comparison
if [ -d data/model_data_smoking ] && [ ! -d data/model_data_smoking.aug10bak ]; then
  mv data/model_data_smoking data/model_data_smoking.aug10bak
fi
echo "=== regenerating data/model_data_smoking from raw_data with build_modeldata_smoking_202607.R ==="
Rscript codes_local/build_modeldata_smoking_202607.R
echo "=== REGEN_SMOKING_DONE ==="
