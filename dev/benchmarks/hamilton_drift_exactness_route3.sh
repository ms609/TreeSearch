#!/bin/bash
# Route-3 marginal-value experiment. One array task per dataset; runs the
# committed harness dev/profiling/drift-exactness-route3-bench.R (arms:
# ratchet / drift_exact / drift_union) from a shared per-seed local optimum,
# cycles knob, matched-wall frontier. mbank_X30754 (180t) is the large/hard
# anchor where a TBR local-opt has real room; the 74-88t sets are easy controls.
# Reuses lib-driftexact (no rebuild).
#SBATCH --job-name=ts-driftexact-route3
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH --time=8:00:00
#SBATCH --array=0-3
#SBATCH --output=/nobackup/%u/TreeSearch/logs/driftexact_route3_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/driftexact_route3_%A_%a.err
set -e
module load r/4.5.1
module load gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
OUTDIR=/nobackup/$USER/TreeSearch/drift-exactness/out
mkdir -p "$OUTDIR"; cd "$REPO"
DATASETS=(Zhu2013 Zanol2014 Dikow2009 mbank_X30754)
DS=${DATASETS[$SLURM_ARRAY_TASK_ID]}
SEEDS=20; if [ "$DS" = "mbank_X30754" ]; then SEEDS=10; fi
GATE_LIB="$LIB" GATE_SEEDS=$SEEDS GATE_CYC=1,2,4,8,16 GATE_DATA="$DS" \
  GATE_OUT="$OUTDIR/route3_${DS}.csv" \
  Rscript dev/profiling/drift-exactness-route3-bench.R
echo "DONE route3 $DS"
