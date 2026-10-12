#!/bin/bash
#SBATCH --job-name=mlp_nonint
#SBATCH --partition=cpu_short
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs_mlp_%j.log
source /gpfs/share/apps/anaconda3/gpu/2022.10/etc/profile.d/conda.sh
conda activate bioinfo
python -c "import sklearn; print(\"sklearn\", sklearn.__version__)"
for s in prev current; do
  echo "=== $s ==="
  python mlp_pipeline.py scratch_mlp_${s}_nonint_123vs0.csv 20
done
echo ALLDONE
