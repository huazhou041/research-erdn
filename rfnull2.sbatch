#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=2 --mem=8G --time=00:30:00
#SBATCH -o rf_null2.log
module load miniconda3/4.6.14
source activate r_env
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript rf_null2.R
echo ALLDONE
