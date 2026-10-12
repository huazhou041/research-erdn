#!/bin/bash
#SBATCH -p cpu_short,cpu_medium
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -c 4
#SBATCH -J clone_202609
#SBATCH -o /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609/logs_smoking/clone_202609.o
set -euo pipefail
SRC=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
DST=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
echo "=== phase1 rsync (code + stable data/results; exclude logs/tmp + changing rfxgb) ==="
rsync -a \
  --exclude="logs" --exclude="logs_*" --exclude="mlp_tmp" --exclude="__pycache__" \
  --exclude="res_rfxgb/" \
  --exclude="rf_xgb_auc_202607.csv" --exclude="rf_xgb_auc_nonint_202607.csv" \
  --exclude="rf_xgb_roc_202607.csv" --exclude="rf_xgb_roc_nonint_202607.csv" \
  "$SRC/" "$DST/"
echo "=== repoint self-referential absolute paths 20260709 -> 202609 in code ==="
grep -rlZ "rnaseq_202609" "$DST" --include="*.R" --include="*.sh" --include="*.py" --include="*.sbatch" 2>/dev/null | xargs -0 -r sed -i "s#rnaseq_202609#rnaseq_202609#g"
echo "=== phase1 size ==="; du -sh "$DST" 2>/dev/null
echo "CLONE_PHASE1_DONE"
