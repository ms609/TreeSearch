#!/bin/bash
#SBATCH --job-name=ts_fe
#SBATCH --partition=shared
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:40:00
#SBATCH --array=0-95%32
# Grid = 2 arms x FE_SEEDS(48) = 96 cells => 0-95.  Each cell = one (seed, arm)
# timed run of FE_CAP_SECONDS, recording the anytime (elapsed, best) trajectory.
#SBATCH --output=/nobackup/%u/TreeSearch/logs/fe_%A_%a.out
#
# Fuse EFFICIENCY array: does enabling pairwise fuse+pool reach the known optimum
# SOONER (wall-clock) than plain restarts?  Reuses the SAME isolated lib the
# basin array built (it already carries the TS_FUSE_PAIRWISE code) -- NO rebuild.
#   aid=$(sbatch --parsable --export=ALL,FE_DATASET=Dikow2009 dev/benchmarks/hamilton_fe_array.sh)
# Pull fe_partials_<ds> back and run FeAnalyze() locally with the CSV target.
set -euo pipefail

module load r/4.5.1 gcc/14.2

: "${REPO:=/nobackup/$USER/TreeSearch-basin}"
: "${LIB:=/nobackup/$USER/TreeSearch-basin-lib}"
: "${DEPS_LIB:=/nobackup/$USER/TreeSearch/lib}"
export R_LIBS_USER="$LIB:$DEPS_LIB"
export TS_LIB=$LIB
cd "$REPO"

export FE_MODE=cell
export FE_DATASET=${FE_DATASET:-Dikow2009}
export FE_SEEDS=${FE_SEEDS:-48}
export FE_CAP_SECONDS=${FE_CAP_SECONDS:-60}
export FE_POOL_SUBOPT=${FE_POOL_SUBOPT:-5}
export FE_PARTIAL_DIR=${FE_PARTIAL_DIR:-$REPO/fe_partials_${FE_DATASET}}

mkdir -p "$FE_PARTIAL_DIR" "/nobackup/$USER/TreeSearch/logs"

Rscript dev/benchmarks/fuse_efficiency.R "$SLURM_ARRAY_TASK_ID"
