#!/bin/bash
#SBATCH --job-name=kick-anytime
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=1:30:00
#SBATCH --array=1-125
#SBATCH --output=/nobackup/%u/TreeSearch/logs/kick_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/kick_%A_%a.err

# Mission A precursor -- ratchetPerturbMaxMoves 5 -> auto/scaled corpus ANYTIME test.
# One array task = one (matrix, seed) cell (25 training matrices x 5 seeds = 125),
# each runs 3 arms (fixed5 / auto / scaled) on the same node.
# BUILDLESS test: the engine LIB is built once BEFORE this array (build_and_smoke.sh);
# this script only runs the harness against that pre-built lib.

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

REPO=/nobackup/$USER/TreeSearch-kick          # current source (has harness + Jul-14 catalogue)
LIB=/nobackup/$USER/TreeSearch/kick-lib       # freshly-built current engine
DEPS=/nobackup/$USER/TreeSearch/lib           # complete dependency library (ape/TreeTools/...)
OUTDIR=/nobackup/$USER/TreeSearch/kick_anytime_out
NEOTRANS=/nobackup/$USER/neotrans/inst/matrices
mkdir -p "$OUTDIR" /nobackup/$USER/TreeSearch/logs

[ -f "$LIB/TreeSearch/DESCRIPTION" ] || { echo "FATAL: engine lib not pre-built at $LIB"; exit 1; }

cd "$REPO"
export R_LIBS_USER="$LIB:$DEPS"               # current TreeSearch first, then deps
TS_LIB="$LIB" NEOTRANS_DIR="$NEOTRANS" \
  CAT_CSV="$REPO/dev/benchmarks/mbank_catalogue.csv" \
  OUT_DIR="$OUTDIR" TASK_ID="${SLURM_ARRAY_TASK_ID:-1}" N_SEEDS=5 \
  Rscript dev/benchmarks/kick_anytime.R 2>&1

echo "Task ${SLURM_ARRAY_TASK_ID} done: $(date)"
