#!/bin/bash
#SBATCH --job-name=calibplatt
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH -o calib_platt_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript calib_cv.R && python3 calib_to_compact.py && echo PLATT_PIPELINE_DONE
