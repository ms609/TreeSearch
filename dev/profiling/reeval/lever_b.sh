#!/bin/bash
# lever_b.sh -- Mission A / Lever B Stage-1 redundancy-curve array (NO engine change).
# Bounded-wall anytime framing: run to maxSeconds under the THOROUGH breadth recipe,
# never converge-stop (targetHits huge) so we observe the full over-search regime.
# Uncapped pool (poolCap 20000) + nThreads 1 => pool_size == count_at_best, exact dedup.
# Reads breadth-vs-wall, first-hit timing, and the main-loop-vs-terminal-enumerator split.
# Lib: tsmarch/lib = current cpp-search tip (7fd61271), -march fast, trajectory-identical.
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH -t 01:45:00
#SBATCH --array=0-11
#SBATCH --job-name=leverb
#SBATCH --output=/nobackup/pjjg18/leverb/logs/leverb_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/leverb/logs/leverb_%A_%a.err
set -e
module load r/4.5.1 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/pjjg18/TreeSearch/lib   # deps (TreeTools, ...)
export TS_LIB=/nobackup/pjjg18/tsmarch/lib           # current-tip TreeSearch
BASE=/nobackup/pjjg18/leverb
MATRICES=/nobackup/pjjg18/neotrans/inst/matrices
DSETS=(project3733 project970 project2668 project510)
SEEDS=(1 2 3)
i=$SLURM_ARRAY_TASK_ID
d=${DSETS[$((i % 4))]}
s=${SEEDS[$((i / 4))]}
echo "cell $i: dset=$d seed=$s host=$(hostname) start=$(date)"
Rscript $BASE/lever_b_curve.R $MATRICES/$d.nex ew thorough 1000 $s $BASE/out 20000 3600
echo "cell $i done=$(date)"
