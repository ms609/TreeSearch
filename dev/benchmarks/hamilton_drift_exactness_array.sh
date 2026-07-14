#!/bin/bash
# Authoritative matched-wall A/B gate for the exact-directional EW drift scorer.
# One array task per dataset; each runs the committed harness
# dev/profiling/drift-exactness-gate-bench.R with production seeds/budgets,
# union (default) vs TS_DRIFT_EXACT, from identical per-seed local-opt starts.
# Submit chained: sbatch --dependency=afterok:<build-jobid> this.sh
#SBATCH --job-name=ts-driftexact-gate
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-2
#SBATCH --output=/nobackup/%u/TreeSearch/logs/driftexact_gate_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/driftexact_gate_%A_%a.err

set -e
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
OUTDIR=/nobackup/$USER/TreeSearch/drift-exactness/out
mkdir -p "$OUTDIR"
cd "$REPO"                       # data-raw/*.nex + the harness live here

DATASETS=(Zhu2013 Zanol2014 Dikow2009)
DS=${DATASETS[$SLURM_ARRAY_TASK_ID]}

# 30 seeds; nCycles sweep extended into the saturated regime to resolve the
# plateau-vs-climb crossover the 8-seed local pilot showed.
GATE_LIB="$LIB" \
GATE_SEEDS=30 \
GATE_NCYC=8,16,32,64,128,256 \
GATE_DATA="$DS" \
GATE_OUT="$OUTDIR/${DS}.csv" \
  Rscript dev/profiling/drift-exactness-gate-bench.R

echo "DONE $DS -> $OUTDIR/${DS}.csv"
