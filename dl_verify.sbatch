#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=2 --mem=12G --time=00:40:00
#SBATCH -o dl_verify.log
module load miniconda3/4.6.14
source activate r_env
unset JAVA_HOME
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
# Exact cell that failed 10/10 before the fix (intersection = fastest folds = tightest race window)
echo "===== VERIFY current/intersection/normalized fpkm 123vs0 ====="
Rscript slurm/deep-drug-response/tune_drug_dl.R recurrence \
  data/model_data/current/intersection/normalized/fpkm_123vs0.RDS \
  slurm/deep-drug-response/hyperparams/grid_dl.json \
  res_dl/VERIFY_isect -v 200 -k 10 --nboot 2 --max_models 4 --inner_nfolds 3 --max_mem 8G --cores 2
echo "EXIT=$?"
echo ALLDONE
