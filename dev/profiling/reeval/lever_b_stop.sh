#!/bin/bash
# lever_b_stop.sh -- Mode (i) "certain it's an MPT": targetHits confidence-stop sweep.
# The score IMPROVES over the first ~6-12 reps (first-hit != true optimum), so a safe
# mode-(i) stop must fire AFTER the score plateaus + is re-confirmed. targetHits = stop
# after N reps hit the current best without improvement (improvement resets the count).
# TEST: does a tight targetHits reach the SAME optimum FASTER on easy-reach (510/970)
# WITHOUT stopping prematurely (worse score) on HARD-reach (4138/5432)?
# poolCap=100 (production; breadth is not the point here). maxSeconds=2700 = safety bound.
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=6G
#SBATCH -t 01:05:00
#SBATCH --array=0-23
#SBATCH --job-name=leverbstop
#SBATCH --output=/nobackup/pjjg18/leverb/logs/stop_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/leverb/logs/stop_%A_%a.err
set -e
module load r/4.5.1 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/pjjg18/TreeSearch/lib
export TS_LIB=/nobackup/pjjg18/tsmarch/lib
DSETS=(project510 project970 project4138 project5432)
THITS=(4 12 999)
i=$SLURM_ARRAY_TASK_ID
seed=$(( i % 2 + 1 ))
th=${THITS[$(( (i / 2) % 3 ))]}
d=${DSETS[$(( i / 6 ))]}
echo "cell $i: dset=$d targetHits=$th seed=$seed start=$(date)"
Rscript /nobackup/pjjg18/leverb/lever_b_curve.R /nobackup/pjjg18/neotrans/inst/matrices/$d.nex ew thorough 1000 $seed /nobackup/pjjg18/leverb/out 100 2700 $th
echo "cell $i done=$(date)"
