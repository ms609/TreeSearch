#!/bin/bash
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/timing_%x_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/timing_%x_%j.err

# Per-dataset TreeSearch-vs-TNT wall-clock timing (run-only: reuses the
# pre-built lib + staged 64-bit TNT). One dataset per job (TS_DATASET via
# --export) so results land independently. TNT is static -> cache the CSV.
module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export LD_LIBRARY_PATH=/nobackup/$USER/TreeSearch/tnt/TNT-bin:$LD_LIBRARY_PATH
export TERM=xterm
export TNT_EXE=/nobackup/$USER/TreeSearch/tnt/TNT-bin/tnt

LIB=/nobackup/$USER/TreeSearch/lib
OUTDIR=/nobackup/$USER/TreeSearch/timing_results
HARNESS=/nobackup/$USER/TreeSearch/scripts/hamilton_timing.R
mkdir -p "$OUTDIR" /nobackup/$USER/TreeSearch/logs

echo "=== Timing: ${TS_DATASET} | $(date) | node $(hostname) ==="
echo "TreeSearch: $(Rscript -e ".libPaths(c(\"$LIB\",.libPaths())); cat(as.character(packageVersion(\"TreeSearch\")))" 2>/dev/null)"
echo "TNT: $TNT_EXE"

TS_LIB="$LIB" TS_DATASET="$TS_DATASET" OUTDIR="$OUTDIR" NSEED="${NSEED:-3}" \
  Rscript "$HARNESS"

echo "Completed: $(date)"
ls -lh "$OUTDIR/timing_${TS_DATASET}.csv" 2>/dev/null
