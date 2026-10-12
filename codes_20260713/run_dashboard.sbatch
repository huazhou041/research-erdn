#!/bin/bash
#SBATCH --partition=cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -J dash_fsmethods
#SBATCH -o logs/dashboard_fsmethods.o

set +u
module load miniconda3/4.6.14
source activate r_env
set -u

cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods.R
