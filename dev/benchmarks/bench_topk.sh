#!/bin/bash
#SBATCH --job-name=bench-topk
#SBATCH --output=/nobackup/pjjg18/ts-bench/results/bench_topk_%j.log
#SBATCH --error=/nobackup/pjjg18/ts-bench/results/bench_topk_%j.err
#SBATCH -n 1
#SBATCH --time=2:00:00
#SBATCH --mem=8000M
#SBATCH -p shared

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BENCH=/nobackup/pjjg18/ts-bench
LIB=$BENCH/lib-cid
export R_LIBS="$LIB:$BENCH/lib-baseline"

cd $BENCH/results

echo "Starting topK benchmark at $(date)"
Rscript $BENCH/TreeSearch/dev/benchmarks/bench_topk.R
echo "Finished at $(date)"
