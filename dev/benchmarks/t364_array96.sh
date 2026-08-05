#!/bin/bash
# T-364 battery, SUPPLEMENTARY sweep at the PRODUCTION replicate budget (96).
#
# Sharp mechanistic test, not a repeat of the main battery.  The main battery runs
# maxReplicates = 24 and shows arm 2 exhausting it while arms 1 and 3 converge
# early.  If arm 2's cost really is "a blocked replicate scores no hit, so the
# convergence rule never trips", then its cost must GROW with the budget -- the
# ratio is bounded by maxrep/reps_converged.  If instead arm 2 were uniformly
# slower per unit work, the ratio would be budget-INVARIANT.  The two hypotheses
# make opposite predictions here, which is the point of running it.
#
# Restricted to small+medium so 96 replicates x 3 arms stays affordable.
#SBATCH --job-name=t364-b96
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t364b96_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t364b96_%A_%a.err

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BASE=/nobackup/$USER/TreeSearch
export NEOTRANS_DIR=/nobackup/$USER/neotrans/inst/matrices
export CAT_CSV=$BASE-t364arm3/dev/benchmarks/mbank_catalogue.csv
export COMMON=$BASE/t364/t364_common.R
export CONS_DIR=$BASE/t364/cons
export OUT_DIR=$BASE/t364/batt96
export TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
export MAXREP=96
export CAP_S=1800
export TIERS=small,medium
mkdir -p "$OUT_DIR"

DEPLIB=$BASE/lib
cd "$BASE/t364" || exit 1

off=$(( (TASK_ID - 1) % 3 ))
for k in 0 1 2; do
  a=$(( (off + k) % 3 + 1 ))
  echo "--- task $TASK_ID, run $((k+1))/3: ARM $a ---"
  TS_LIB=$BASE/t364-lib$a ARM=$a \
    R_LIBS="$BASE/t364-lib$a:$DEPLIB" \
    Rscript "$BASE/t364/t364_cell.R" || echo "ARM $a FAILED (rc=$?)"
done
echo "task $TASK_ID done"
