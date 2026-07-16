#!/bin/bash
# sprFirst gate -- JOB 1: warmup -> first-TBR handoff on TRAINING-split MorphoBank
# (EW-recoded), size-stratified & weighted to the LARGE end (the discriminator).
# One catalogue key per array task.  Reuses lib-driftexact (exact-SPR default +
# TS_SPR_UNION kill-switch) -- NO rebuild.  TRAINING split only (sequestration).
#
# Submit (from a Hamilton login node):
#   sbatch dev/benchmarks/hamilton_sprfirst_gate.sh
# Then aggregate locally: Rscript dev/profiling/sprfirst-analyze.R
#SBATCH --job-name=ts-sprfirst
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --array=0-5
#SBATCH --output=/nobackup/%u/TreeSearch/logs/sprfirst_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/sprfirst_%A_%a.err
set -e
module load r/4.5.1; module load gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/$USER/TreeSearch/lib
LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
OUTDIR=/nobackup/$USER/TreeSearch/sprfirst/out
mkdir -p "$OUTDIR"; cd "$REPO"
# Size-stratified, weighted large; degenerate low-char keys excluded (project3212
# 10c, project861 32c -> trivial convergence / zero signal).
KEYS=(project2183 project6403 project3763 project1024 project2124 project3617)
KEY=${KEYS[$SLURM_ARRAY_TASK_ID]}
GATE_LIB="$LIB" GATE_KEY="$KEY" GATE_CAT=/nobackup/$USER/floor/mbank_catalogue.csv \
  GATE_SEEDS=5 \
  GATE_OUT="$OUTDIR/sprfirst_${KEY}.csv" \
  Rscript dev/profiling/sprfirst-gate-bench.R
echo "DONE sprfirst $KEY"
