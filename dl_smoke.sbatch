#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=4 --mem=24G --time=01:00:00
#SBATCH -o dl_smoke.log
module load miniconda3/4.6.14
source activate r_env
unset JAVA_HOME
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
echo "===== A) clinics (5 feats), nboot=2, max_models=4 ====="
time Rscript slurm/deep-drug-response/tune_drug_dl.R recurrence \
  data/model_data/prev/wholeset/normalized/clinics_123vs0.RDS \
  slurm/deep-drug-response/hyperparams/grid_dl.json \
  res_dl/SMOKE_clinics -v 100 -k 10 --nboot 2 --max_models 4 --max_mem 16G --cores 4
echo "===== B) gene panel varcut=200, nboot=1, max_models=4 (timing probe) ====="
time Rscript slurm/deep-drug-response/tune_drug_dl.R recurrence \
  data/model_data/prev/wholeset/normalized/fpkm_123vs0.RDS \
  slurm/deep-drug-response/hyperparams/grid_dl.json \
  res_dl/SMOKE_fpkm -v 200 -k 10 --nboot 1 --max_models 4 --max_mem 16G --cores 4
echo "ALLDONE"
