#!/bin/bash
#SBATCH --job-name=roc_build
#SBATCH --partition=cpu_medium,cpu_short
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH -o roc_build_%j.log
set +u; module load miniconda3/4.6.14; source activate r_env; set -u
cd /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
step(){ echo "### $1"; }
step "ROC full GLM";     Rscript codes_20260713/make_roc_curves_fsmethods.R        && echo OK_ROC_GLM      || echo FAIL_ROC_GLM
step "ROC smoking GLM";  Rscript make_roc_curves_fsmethods_smoking.R               && echo OK_ROC_GLM_SMK  || echo FAIL_ROC_GLM_SMK
step "ROC full rfxgb";   Rscript make_roc_rfxgb_full.R                             && echo OK_ROC_RFXGB    || echo FAIL_ROC_RFXGB
step "ROC smoking rfxgb";Rscript make_roc_rfxgb_smoking.R                          && echo OK_ROC_RFXGB_SMK|| echo FAIL_ROC_RFXGB_SMK
echo "===== ROC CSVs present ====="
for f in results_roc_curves_fsmethods.csv results_roc_curves_fsmethods_smoking.csv rf_xgb_roc_202607.csv rf_xgb_roc_smoking.csv; do [ -f "$f" ] && echo "HAVE $f ($(wc -l <$f))" || echo "MISS $f"; done
step "BUILD full";       Rscript codes_20260713/make_dashboard_fsmethods.R         && echo OK_BUILD_FULL   || echo FAIL_BUILD_FULL
step "BUILD smoking";    Rscript codes_20260713/make_dashboard_fsmethods_smoking.R && echo OK_BUILD_SMK    || echo FAIL_BUILD_SMK
step "BUILD cmp toggle"; Rscript codes_20260713/make_dashboard_fsmethods_cmp.R     && echo OK_BUILD_CMP    || echo FAIL_BUILD_CMP
echo "===== HTML outputs ====="
ls -la htmls_20260713/*.html 2>/dev/null || echo "(no html yet)"
echo "===== ROC_BUILD_ALL_DONE ====="
