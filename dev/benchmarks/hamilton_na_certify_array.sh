#!/bin/bash
# NA-certification gate panel: one array task per (dataset x seed x tabuSize),
# each running all five arms on the same node so quality and wall comparisons
# are not contaminated by node variance.
# Consumes the $LIB built by hamilton_na_certify_build.sh.
# Submit: sbatch --dependency=afterok:<buildjob> hamilton_na_certify_array.sh
#
# Grid = expand.grid(dataset, seed, tabu) in bench_na_certify_cell.R's order:
#   30 datasets x 5 seeds x 2 tabu = 300 tasks (0-299).
# All 30 bundled inapplicable matrices, NATIVE (never recoded "-"->"?"): the NA
# path is the whole point.  Four matrices would not support a paired test at
# matrix level (pair-on-matrices-not-seeds), which is why the panel is all 30.
#
# tabuSize is an arm dimension: do_reroot -- the gate on exact_verify_sweep --
# requires tabu_size == 0, and the shipped presets set 100/200.  Running both
# separates "the lever is worthless" from "the lever is unreachable under the
# default preset".
#
# No %N concurrency cap: let SLURM and fairshare pace the array (no-array-throttle).
#SBATCH --job-name=ts-nacert
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
# Headroom for the worst cell: 88-tip Dikow2009 at tabuSize = 0 with
# certification ON runs the O(n^3) sweep at every whole-tree convergence, and the
# matched-wall arms then re-spend that same budget.  Well under `shared`'s 72 h,
# so no `long` partition (hamilton-partition-choice).
#SBATCH --time=12:00:00
#SBATCH --array=0-299
#SBATCH --output=/nobackup/%u/TreeSearch/logs/nacert_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/nacert_%A_%a.err

module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DEPLIB=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-nacert
REPO=/nobackup/$USER/TreeSearch-nacert
export R_LIBS_USER="$LIB:$DEPLIB"
export TS_LIB=$LIB
export PARTIAL_DIR=/nobackup/$USER/TreeSearch/na_certify_partials
export TS_REPS=8
export TS_SEEDS="1 2 3 4 5"
export TS_TABU="100 0"
# EW (concavity = Inf) is the MaximizeParsimony default and shows the effect
# (Zanol2014: 1321 certified vs 1326 gated).  Set TS_CONCAVITY=10 for an IW pass,
# which is the regime dev/profiling/na-exact-verify-dominates.md profiled.
export TS_CONCAVITY=Inf
# TS_DATASETS unset => all 30 bundled inapplicable matrices.

mkdir -p "$PARTIAL_DIR" /nobackup/$USER/TreeSearch/logs
cd "$REPO" || exit 1
Rscript dev/benchmarks/bench_na_certify_cell.R "$SLURM_ARRAY_TASK_ID"
