#!/bin/bash
# T-364 battery step 1 on Hamilton: freeze the constraint set, one array task per
# matrix.  Runs against ARM 3 only -- the constraints must not be chosen by the arm
# they will be used to judge, and every arm then consumes the identical frozen set.
#SBATCH --job-name=t364-gen
#SBATCH -p shared
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/t364gen_%A_%a.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/t364gen_%A_%a.err

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BASE=/nobackup/$USER/TreeSearch
export NEOTRANS_DIR=/nobackup/$USER/neotrans/inst/matrices
export CAT_CSV=$BASE-t364arm3/dev/benchmarks/mbank_catalogue.csv
export COMMON=$BASE/t364/t364_common.R
export OUT_DIR=$BASE/t364/cons
export TS_LIB=$BASE/t364-lib3
export R_LIBS="$BASE/t364-lib3:$BASE/lib"
export TASK_ID=${SLURM_ARRAY_TASK_ID:-1}
export REF_REP=${REF_REP:-12}
mkdir -p "$OUT_DIR"
cd "$BASE/t364" || exit 1
Rscript "$BASE/t364/t364_gen.R"
