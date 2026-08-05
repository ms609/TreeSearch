#!/bin/bash
# Sector-godrift isolation gate (route B). One array task per dataset; runs the
# committed harness dev/profiling/drift-exactness-sector-bench.R (3 arms:
# none/union/exact) at production seeds/replicates. Reuses lib-driftexact (code
# unchanged since f7f17e16 -- docs-only commits since), so NO rebuild: just
# `git pull` the repo for the new harness, then submit this.
#SBATCH --job-name=ts-driftexact-sector
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH --time=8:00:00
#SBATCH --array=0-3
#SBATCH --output=/nobackup/%u/TreeSearch/logs/driftexact_sector_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/driftexact_sector_%A_%a.err
set -e
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
OUTDIR=/nobackup/$USER/TreeSearch/drift-exactness/out
mkdir -p "$OUTDIR"
cd "$REPO"
DATASETS=(Zhu2013 Zanol2014 Dikow2009 mbank_X30754)
DS=${DATASETS[$SLURM_ARRAY_TASK_ID]}
# Fewer seeds for the 180-tip mbank (each thorough search is far slower).
SEEDS=15; if [ "$DS" = "mbank_X30754" ]; then SEEDS=8; fi
GATE_LIB="$LIB" GATE_SEEDS=$SEEDS GATE_REPS=1,2,4 GATE_DATA="$DS" \
  GATE_OUT="$OUTDIR/sector_${DS}.csv" \
  Rscript dev/profiling/drift-exactness-sector-bench.R
echo "DONE sector $DS"
