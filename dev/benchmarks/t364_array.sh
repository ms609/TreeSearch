#!/bin/bash
# T-364/T-370 three-arm constrained-search battery: one array task = one
# (matrix, shape, seed) cell, all THREE arms run back-to-back on the SAME node
# inside one wall-clock window.
#
# WHY INTERLEAVED, NOT THREE SEPARATE ARRAYS: a sequential A/B across separately
# scheduled jobs compares different nodes at different times under different
# neighbour load.  Same node, same window, adjacent processes is the only fair
# wall comparison available on a shared partition.
#
# ARM ORDER IS ROTATED BY TASK ID so that no arm systematically runs first (first
# run pays page-cache warm-up on the matrix file and the .so).
#
# Submit with:  sbatch --array=1-N --export=ALL,TIER_SET=main|xlarge t364_array.sh
#SBATCH --job-name=t364-batt
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t364batt_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t364batt_%A_%a.err

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BASE=/nobackup/$USER/TreeSearch
export NEOTRANS_DIR=/nobackup/$USER/neotrans/inst/matrices
export CAT_CSV=$BASE-t364arm3/dev/benchmarks/mbank_catalogue.csv
export COMMON=$BASE/t364/t364_common.R
export CONS_DIR=$BASE/t364/cons
export OUT_DIR=$BASE/t364/batt
export TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
mkdir -p "$OUT_DIR"

DEPLIB=$BASE/lib
cd "$BASE/t364" || exit 1

# Rotate: task 1 -> 1,2,3 ; task 2 -> 2,3,1 ; task 3 -> 3,1,2 ; ...
off=$(( (TASK_ID - 1) % 3 ))
for k in 0 1 2; do
  a=$(( (off + k) % 3 + 1 ))
  echo "--- task $TASK_ID, run $((k+1))/3: ARM $a ---"
  TS_LIB=$BASE/t364-lib$a ARM=$a \
    R_LIBS="$BASE/t364-lib$a:$DEPLIB" \
    Rscript "$BASE/t364/t364_cell.R" || echo "ARM $a FAILED (rc=$?)"
done
echo "task $TASK_ID done"
