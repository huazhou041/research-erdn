#!/bin/bash
#SBATCH -p cpu_short,cpu_medium
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH -c 4
#SBATCH -J clone_cd_202609
#SBATCH -o /gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609/logs_smoking/clone_code_data.o
set -euo pipefail
SRC=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
DST=/gpfs/data/tsirigoslab/home/zhouh05/segal/rnaseq_202609
mkdir -p "$DST"
echo "=== step 1: CODE files only ==="
rsync -am \
  --exclude="res/" --exclude="res_cox/" --exclude="res_dl/" --exclude="res_rfxgb/" \
  --exclude="res_rfxgb_smoking/" --exclude="res_smoking/" --exclude="results/" \
  --exclude="htmls_20260713/" --exclude="logs/" --exclude="logs_*/" \
  --exclude="mlp_tmp/" --exclude="__pycache__/" --exclude="data/" \
  --include="*/" \
  --include="*.R" --include="*.sh" --include="*.py" --include="*.sbatch" \
  --include="*.json" --include="*.html" --include="*.js" --include="ever_smoker_ids.txt" \
  --include="*.md" \
  --exclude="*" \
  "$SRC/" "$DST/"
echo "=== step 2: DATA dir in full ==="
rsync -a "$SRC/data/" "$DST/data/"
echo "=== step 3: repoint absolute paths 20260709 -> 202609 in code ==="
grep -rlZ "rnaseq_202609" "$DST" --include="*.R" --include="*.sh" --include="*.py" --include="*.sbatch" 2>/dev/null | xargs -0 -r sed -i "s#rnaseq_202609#rnaseq_202609#g"
echo "=== sizes / counts ==="
du -sh "$DST" "$DST/data" 2>/dev/null
echo -n "code files: "; find "$DST" -type f \( -name "*.R" -o -name "*.sh" -o -name "*.py" -o -name "*.sbatch" \) -not -path "*/data/*" | wc -l
echo -n "residual result dirs (want 0): "; ls -d "$DST"/res* "$DST"/results "$DST"/htmls_20260713 2>/dev/null | wc -l
echo "CLONE_CD_DONE"
