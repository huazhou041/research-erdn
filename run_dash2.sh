#!/bin/bash
#SBATCH --job-name=dash
#SBATCH --partition=cpu_dev
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --output=dash_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate r_env
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods.R
echo DASH_DONE
