#!/bin/bash
# Build TreeSearch from the certify_unrooted branch ONCE into a dedicated
# read-only library, so the panel array never recompiles.  Deps (TreeTools/Rcpp/
# ape/...) are reused from the existing populated $DEPLIB.
#
# The branch adds a field to TBRParams in src/ts_tbr.h, so this MUST be a clean
# build -- a stale .o against the old struct layout is a recorded trap
# (stale-object-abi-gotcha).  `rm -f src/*.o src/*.so` below is that guarantee;
# do not replace it with an incremental install.
#
# Submit first; chain the array on afterok of this job.
#SBATCH --job-name=ts-nacert-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/nacert_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/nacert_build_%j.err

module load r/4.5.1
module load gcc/14.2

BRANCH=${TS_BRANCH:-claude/angry-driscoll-25bfe3}
REPO=/nobackup/$USER/TreeSearch-nacert
DEPLIB=/nobackup/$USER/TreeSearch/lib            # has TreeTools/Rcpp/ape/...
LIB=/nobackup/$USER/TreeSearch/lib-nacert        # fresh target for this build
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs
export R_LIBS_USER="$LIB:$DEPLIB"

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO" || { echo "FATAL: clone failed"; exit 1; }
fi
cd "$REPO" || { echo "FATAL: no $REPO"; exit 1; }
git fetch origin "$BRANCH" && (git checkout "$BRANCH" && git reset --hard "origin/$BRANCH")
echo "Git HEAD: $(git log --oneline -1)"

# Header change => clean build, not incremental (stale-object-abi-gotcha).
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"

# Fail loudly if the field never made it in: every arm would silently be the
# baseline and the panel would report a null result for the wrong reason.
Rscript -e '
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  library(TreeSearch)
  data("inapplicable.phyData", package = "TreeSearch")
  d <- inapplicable.phyData[["Vinther2008"]]
  Sys.setenv(TS_NA_NOCERTIFY = "1")
  set.seed(1)
  r <- MaximizeParsimony(d, maxReplicates = 1L, targetHits = 9999L,
                         nThreads = 1L, verbosity = 0L, tabuSize = 0L)
  nd <- attr(r, "naDiag")
  if (is.null(nd)) stop("FATAL: no naDiag -- wrong build")
  cat("gate smoke: n_evs =", nd$n_evs, " n_evs_skipped =", nd$n_evs_skipped, "\n")
  if (nd$n_evs_skipped <= 0) stop("FATAL: gate never fired -- wrong build")
' || exit 1

exit $rc
