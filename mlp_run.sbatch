#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=4 --mem=16G --time=02:00:00
#SBATCH -o mlp_run.log
module load miniconda3/4.6.14
source activate r_env
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
mkdir -p mlp_tmp
for s in prev current; do for co in nonintersection novel; do
  f=data/model_data/$s/$co/normalized/fpkm_123vs0.RDS
  csv=mlp_tmp/${s}_${co}.csv
  [ -f "$csv" ] || Rscript -e "m<-readRDS('$f'); m<-m[!is.na(m\$recurrence),]; write.csv(m, '$csv', row.names=FALSE)" 2>/dev/null
done; done
echo "=== MLP (sklearn) pooled OOF AUC + permutation null ==="
for s in prev current; do for co in nonintersection novel; do
  python3 mlp_pipeline.py mlp_tmp/${s}_${co}.csv 25
done; done
echo ALLDONE
