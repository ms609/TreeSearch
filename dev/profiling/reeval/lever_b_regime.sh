#!/bin/bash
# lever_b_regime.sh -- Mode (ii) / thorough GO-NO-GO probe (zero engine change).
#
# THE GATING QUESTION (advisor, 2026-07-19): does ANY real dataset show restart
# breadth-discovery FLATTENING mid-run while replicates keep coming?  The mode-(ii)
# restart-stop lever ("stop restarts once new island discovery stalls") has VALUE
# ONLY in an INTERMEDIATE landscape -- multi-island BUT saturating.  The two ends
# are already settled: 510 saturates instantly (auto/enumerator covers it); 3733
# climbs indefinitely (nothing to reclaim).  We do NOT know an intermediate dataset
# is in the corpus.  This probe reads the uncapped restart-breadth curve for a
# spread of candidates + both controls and looks for a mid-run flatten.
#
# METHOD = the uncapped no-code trick (poolCap huge + never-stop): pool_size-vs-rep
# during the main loop IS the uncapped restart-breadth curve.  analyze_disc.R segments
# pool growth by best-so-far.  Read each curve:
#   * instant flatten (saturates by a few reps)        -> like 510  (auto already covers)
#   * indefinite climb (new distinct/rep never stalls) -> like 3733 (nothing to reclaim)
#   * MID-RUN flatten while reps keep coming            -> THE REGIME -> justifies a build
# Verdict logic: NONE flatten => mode (ii) closes by empty-regime (no instrument written);
#                ONE flattens => that dataset justifies the shadow-counter build.
#
# Cheapest matrix first (file size = per-rep-cost proxy): 4138 (37KB) runs the most
# reps -> the only cell that can expose a LARGE-k flatten (cf. its rep99->205 dead
# window).  Bigger files run fewer reps; "still climbing when budget ran out" is still
# a valid (if weaker) read = not-saturating-in-budget.
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH -t 03:15:00
#SBATCH --array=0-7
#SBATCH --job-name=leverbreg
#SBATCH --output=/nobackup/pjjg18/leverb/logs/regime_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/leverb/logs/regime_%A_%a.err
set -e
module load r/4.5.1 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/pjjg18/TreeSearch/lib
export TS_LIB=/nobackup/pjjg18/tsmarch/lib
M=/nobackup/pjjg18/neotrans/inst/matrices
OUT=/nobackup/pjjg18/leverb/regime
mkdir -p "$OUT"

# (dataset seed role) -- cheapest-first; controls flagged
CELLS=(
  "project4138.nex 1 intermediate-cheap"   # 37KB, known 327-plateau+dead-window; observe at TRUE best
  "project4138.nex 2 intermediate-cheap"   # seed replicate
  "project2183.nex 1 candidate"            # 252KB, 2nd-cheapest, structure unknown
  "project3733.nex 1 NEG-control"          # climbs indefinitely (0.85-0.95 distinct/rep)
  "project510.nex  1 POS-control"          # saturates ~rep 4 (single island)
  "project2668.nex 1 candidate"            # 398KB moderate
  "project970.nex  1 candidate"            # 623KB, multi-island (breadth swings 6<->48)
  "project4327.nex 1 candidate"            # 414KB, spread
)
CELL=(${CELLS[$SLURM_ARRAY_TASK_ID]})
DSET=${CELL[0]}; SEED=${CELL[1]}; ROLE=${CELL[2]}

echo "regime probe cell=$SLURM_ARRAY_TASK_ID dset=$DSET seed=$SEED role=$ROLE start=$(date)"
# args: matrix mode strategy maxRep seed outdir poolCap maxSec targetHits(NA=never)
Rscript /nobackup/pjjg18/leverb/lever_b_curve.R "$M/$DSET" ew thorough 500 "$SEED" "$OUT" 20000 10800 NA
echo "regime probe cell=$SLURM_ARRAY_TASK_ID done=$(date)"
