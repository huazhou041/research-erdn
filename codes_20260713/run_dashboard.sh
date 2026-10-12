#!/bin/bash
set +u
module load miniconda3/4.6.14
source activate r_env
set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods.R
