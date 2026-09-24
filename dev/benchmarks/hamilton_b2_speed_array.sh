#!/bin/bash
#SBATCH --job-name=b2_speed
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --array=0-895%32
# Grid = 2 arms x 14 keys x 32 seeds = 896 cells (0-895). One cell = one timed
# (arm, key, seed) MaximizeParsimony recording the anytime (elapsed, best)
# trajectory.  10 CHAR-POOR training datasets (high collapse density at the
# optimum + expensive search = where the lever can plausibly win) + 4 char-RICH
# controls (expect ON==OFF).  Does TS_COLLAPSE_AGGRESSIVE reach the optimum
# SOONER?  Pull b2_speed_partials and run b2_speed_analyze.R.
#SBATCH --output=/nobackup/%u/TreeSearch/logs/b2_%A_%a.out
set -euo pipefail

module load r/4.5.1 gcc/14.2

: "${REPO:=/nobackup/$USER/TreeSearch-basin}"
: "${LIB:=/nobackup/$USER/TreeSearch-basin-lib}"
: "${DEPS_LIB:=/nobackup/$USER/TreeSearch/lib}"
export R_LIBS_USER="$LIB:$DEPS_LIB"
export TS_LIB=$LIB
cd "$REPO"

export MBANK_DIR=${MBANK_DIR:-/nobackup/$USER/neotrans/inst/matrices}
export MBANK_CAT=${MBANK_CAT:-/nobackup/$USER/floor/mbank_catalogue.csv}
export B2_SEEDS=${B2_SEEDS:-32}
export B2_BUDGET=${B2_BUDGET:-90}
export B2_STRATEGY=${B2_STRATEGY:-thorough}
export PARTIAL_DIR=${PARTIAL_DIR:-$REPO/b2_speed_partials}
# 10 char-poor (ratio<1.5, 54-188 tips) + 4 char-rich controls.  Order MUST match
# the analysis; keys verbatim from mbank_catalogue.csv (training split).
export B2_KEYS="project784 project3354_(1) project4626 project3512_(2) project954 project3906 project2144 project2771 project4598 project3287_Cassidulidae_without_partial_uncertainties project2084_(1) project3887_(1) syab07203 project2289"

mkdir -p "$PARTIAL_DIR" "/nobackup/$USER/TreeSearch/logs"
Rscript dev/benchmarks/b2_speed_cell.R "$SLURM_ARRAY_TASK_ID"
