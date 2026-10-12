#!/bin/bash
#SBATCH --partition=cpu_short
#SBATCH --nodes=1 --ntasks=4 --mem=16G --time=00:45:00
#SBATCH -o rf_null.log
module load miniconda3/4.6.14
source activate r_env
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript rf_null.R
echo ALLDONE
