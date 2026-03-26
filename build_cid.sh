#!/bin/bash
#SBATCH --job-name=build-cid
#SBATCH --output=/nobackup/pjjg18/ts-bench/results/build_cid.log
#SBATCH --error=/nobackup/pjjg18/ts-bench/results/build_cid.err
#SBATCH -n 4
#SBATCH --time=0:30:00
#SBATCH --mem=4000M
#SBATCH -p shared

module load r/4.5.1
module load gcc/14.2

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

BENCH=/nobackup/pjjg18/ts-bench
REPO=$BENCH/TreeSearch
LIB=$BENCH/lib-cid

# Make deps visible to R CMD INSTALL
export R_LIBS="$LIB:$BENCH/lib-baseline"

mkdir -p "$LIB"

rm -f "$REPO"/src/*.o "$REPO"/src/*.so

TMPDIR=$(mktemp -d)
echo "Building tarball in $TMPDIR..."
(cd "$TMPDIR" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$REPO")

echo "Installing to $LIB..."
R CMD INSTALL --library="$LIB" "$TMPDIR"/TreeSearch_*.tar.gz
BUILD_STATUS=$?

rm -rf "$TMPDIR"

echo "Build exit status: $BUILD_STATUS"

if [ $BUILD_STATUS -eq 0 ]; then
  echo "Verifying build..."
  Rscript -e "
    .libPaths(c('$LIB', '$BENCH/lib-baseline', .libPaths()))
    library(TreeSearch)
    cat('TreeSearch version:', as.character(packageVersion('TreeSearch')), '\n')
    cat('InfoConsensus available:', exists('InfoConsensus', where='package:TreeSearch', mode='function'), '\n')
    library(TreeTools)
    trees <- lapply(1:10, function(i) RandomTree(20))
    class(trees) <- 'multiPhylo'
    cat('Smoke test: scoring 10 random 20-tip trees...\n')
    result <- InfoConsensus(trees, neverDrop = TRUE, maxSeconds = 10)
    cat('Result class:', class(result), '\n')
    cat('Result tips:', Ntip(result), '\n')
    cat('Build verification PASSED\n')
  "
else
  echo "BUILD FAILED"
fi
