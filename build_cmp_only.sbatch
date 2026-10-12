#!/bin/bash
#SBATCH --job-name=build_cmp
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH -o build_cmp_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods_cmp.R && echo OK_BUILD_CMP || echo FAIL_BUILD_CMP
H=htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking.html
echo "natcommad section present: $(grep -c 'id=\"natcommad\"' "$H")  renderNatcommad: $(grep -c 'function renderNatcommad' "$H")  ncTable: $(grep -c 'id=\"ncTable\"' "$H")"
echo "=== BUILD_CMP_ONLY_DONE ==="
