#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=2 --mem=16G --time=00:40:00
#SBATCH -o glm3_verify.log
module load miniconda3/4.6.14
source activate r_env
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
echo "===== VERIFY novel/3vs0 (4 events) now runs under the nfolds fallback ====="
Rscript slurm/deep-drug-response/tune_drug.R recurrence \
  data/model_data/prev/novel/normalized/fpkm_3vs0.RDS \
  slurm/deep-drug-response/hyperparams/grid_glm.json \
  res/VERIFY3_novel -v 100 --learner GLM -k 10 --nboot 3 --fs_method mad
echo "EXIT=$?"
echo ALLDONE
