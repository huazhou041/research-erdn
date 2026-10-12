#!/bin/bash
#SBATCH --job-name=build_only
#SBATCH --partition=cpu_short,cpu_medium
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH -o build_only_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
Rscript codes_20260713/make_dashboard_fsmethods.R         && echo OK_BUILD_FULL   || echo FAIL_BUILD_FULL
Rscript codes_20260713/make_dashboard_fsmethods_smoking.R && echo OK_BUILD_SMK    || echo FAIL_BUILD_SMK
Rscript codes_20260713/make_dashboard_fsmethods_cmp.R     && echo OK_BUILD_CMP    || echo FAIL_BUILD_CMP
echo "=== natcommad label present in outputs? ==="
for h in htmls_20260713/auc_dashboard_20260713_fsmethods.html htmls_20260713/auc_dashboard_20260713_fsmethods_smoking.html htmls_20260713/auc_dashboard_20260713_fsmethods_fullVsSmoking.html; do
  echo "$h : GLM+ncMAD=$(grep -c 'GLM+ncMAD' "$h" 2>/dev/null) glm_natcommad=$(grep -c 'glm_natcommad' "$h" 2>/dev/null)"
done
echo "=== BUILD_ONLY_DONE ==="
