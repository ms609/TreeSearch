#!/bin/bash
# Recipe-level whole-search test (route B'): drift scorer choice in wall-limited
# thorough on TRAINING-split MorphoBank (EW-recoded). One catalogue key per task.
# Reuses lib-driftexact (no rebuild). TRAINING split only (sequestration).
#SBATCH --job-name=ts-driftexact-recipe
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-7
#SBATCH --output=/nobackup/%u/TreeSearch/logs/driftexact_recipe_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/driftexact_recipe_%A_%a.err
set -e
module load r/4.5.1; module load gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
OUTDIR=/nobackup/$USER/TreeSearch/drift-exactness/out
mkdir -p "$OUTDIR"; cd "$REPO"
# Deterministic mid-size (65-135t) training-split keys, spread across sizes.
KEYS=(project3617 project589 project4624 project2124 project5201 project3807 project3602 project4614)
KEY=${KEYS[$SLURM_ARRAY_TASK_ID]}
GATE_LIB="$LIB" GATE_KEY="$KEY" GATE_CAT=/nobackup/$USER/floor/mbank_catalogue.csv \
  GATE_SEEDS=6 GATE_REPS=2,4,8 \
  GATE_OUT="$OUTDIR/recipe_${KEY}.csv" \
  Rscript dev/profiling/drift-exactness-recipe-bench.R
echo "DONE recipe $KEY"
