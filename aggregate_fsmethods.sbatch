#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH -J agg_fsmethods
#SBATCH -o logs/aggregate_fsmethods.o

set +u
module load miniconda3/4.6.14
source activate r_env
set -u

cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
./aggregate_fsmethods.R
