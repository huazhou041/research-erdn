#!/bin/bash
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=24G
#SBATCH --time=00:40:00
#SBATCH -J smk_cmp
#SBATCH -o logs_smoking/build_cmp.o
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods_cmp.R
echo "=== BUILD_CMP_DONE ==="
