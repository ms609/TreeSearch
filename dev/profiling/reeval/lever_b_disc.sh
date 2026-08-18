#!/bin/bash
# lever_b_disc.sh -- Mode (i) DECIDER: is an adaptive discovery-rate stop buildable, or
# is cheap certainty on hard instances definitively impossible?
# From the mode-(i) sweep: 4138/s1 sat at score 327 across a 62-rep plateau (rep 7->69),
# then IMPROVED to 326 at rep 363 (300 reps later). The pool there was capped at 100, so we
# could not see whether NEW distinct 327-topologies were still being discovered during the
# plateau. This reruns 4138/s1 UNCAPPED + never-stop, past rep 363, so pool[rep] == cumulative
# distinct MPTs at the current best. Read: did distinct discovery stay ELEVATED through rep 362
# (=> "keep going while discovering" would have survived to catch 326 => adaptive stop viable)
# or go DEAD by rep 69 (=> no local/trajectory signal can catch 326 => document-only)?
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH -t 02:30:00
#SBATCH --job-name=leverbdisc
#SBATCH --output=/nobackup/pjjg18/leverb/logs/disc_%A.out
#SBATCH --error=/nobackup/pjjg18/leverb/logs/disc_%A.err
set -e
module load r/4.5.1 2>/dev/null || true
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export R_LIBS_USER=/nobackup/pjjg18/TreeSearch/lib
export TS_LIB=/nobackup/pjjg18/tsmarch/lib
M=/nobackup/pjjg18/neotrans/inst/matrices
echo "disc probe: 4138/s1 uncapped never-stop start=$(date)"
# args: matrix mode strategy maxRep seed outdir poolCap maxSec targetHits(NA=never)
Rscript /nobackup/pjjg18/leverb/lever_b_curve.R $M/project4138.nex ew thorough 500 1 /nobackup/pjjg18/leverb/out 20000 7000 NA
echo "disc probe done=$(date)"
