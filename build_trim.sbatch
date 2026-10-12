#!/bin/bash
#SBATCH --job-name=build_trim
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH -o build_trim_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods_cmp_trim.R && echo OK_BUILD_TRIM || echo FAIL_BUILD_TRIM
echo "=== BUILD_TRIM_DONE ==="
