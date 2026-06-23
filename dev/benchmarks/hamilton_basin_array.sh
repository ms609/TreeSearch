#!/bin/bash
#SBATCH --job-name=ts_basin
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-61%32
# Grid = BD_INDEP_K(30) + 2*BD_PROD_SEEDS(16) = 62 cells => 0-61.
# NB: time budget is per task. INDEP cells (maxReplicates=1) finish in seconds;
# PROD cells run for BD_PROD_SECONDS wall-clock (+ build/teardown), so they set
# the ceiling. Raise --time if BD_PROD_SECONDS is increased.
#SBATCH --output=/nobackup/%u/TreeSearch/logs/basin_%A_%a.out
#
# Basin-diversity array (audit item B1).  One task per grid cell:
#   cells 0 .. BD_INDEP_K-1                 : INDEP restarts (maxReplicates=1)
#   next  2 * BD_PROD_SEEDS cells           : PROD prod_on / prod_off
# Grid size = BD_INDEP_K + 2*BD_PROD_SEEDS.  Set --array=0-(N-1) to match
#   (default below: 50 + 2*8 = 66 => 0-65).
#
# Pipeline (run from a Hamilton login node):
#   bid=$(sbatch --parsable dev/benchmarks/hamilton_build_once.sh)
#   aid=$(sbatch --parsable --dependency=afterok:$bid dev/benchmarks/hamilton_basin_array.sh)
#   # then collect: pull $PARTIAL_DIR back and run BD_MODE=analyze locally
#   # (analyze also runs the TNT optimal-island arm; keep it off the cluster).
set -euo pipefail

module load r/4.5.1 gcc/14.2

# REPO (checkout with the harness scripts) and LIB (installed package) are
# env-overridable so the same script serves the shared layout and an isolated
# per-measurement layout. Defaults match hamilton_basin_build.sh.
: "${REPO:=/nobackup/$USER/TreeSearch-basin}"
: "${LIB:=/nobackup/$USER/TreeSearch-basin-lib}"
: "${DEPS_LIB:=/nobackup/$USER/TreeSearch/lib}"
export R_LIBS_USER="$LIB:$DEPS_LIB"  # TreeSearch from $LIB, deps from $DEPS_LIB
export TS_LIB=$LIB
cd "$REPO"

export BD_MODE=cell
export BD_DATASET=${BD_DATASET:-Zanol2014}
export BD_INDEP_K=${BD_INDEP_K:-30}
export BD_PROD_REPS=${BD_PROD_REPS:-9999}     # high cap; the real budget is wall-clock
export BD_PROD_SEEDS=${BD_PROD_SEEDS:-16}
# B1 pairwise-fuse A/B (2026-06-22): prod_on enables TS_FUSE_PAIRWISE recombination
# into suboptimal pool recipients; prod_off never fuses. Both arms hold
# poolSuboptimal identical, and are TIME-MATCHED on WALL-CLOCK (not replicates),
# so fuse-ON wins only if recombination beats the restarts it crowds out.
export BD_PROD_SECONDS=${BD_PROD_SECONDS:-30}  # per-run wall-clock budget
export BD_POOL_SUBOPT=${BD_POOL_SUBOPT:-5}     # suboptimal-recipient headroom (both arms)
export BD_PARTIAL_DIR=${BD_PARTIAL_DIR:-$REPO/basin_partials_${BD_DATASET}}
export BD_TNT=0   # TNT executable lives on the workstation, not Hamilton

mkdir -p "$BD_PARTIAL_DIR" "/nobackup/$USER/TreeSearch/logs"

Rscript dev/benchmarks/basin_cell.R "$SLURM_ARRAY_TASK_ID"
