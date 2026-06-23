#!/bin/bash
# Build TreeSearch from the NA-IW-x4 branch (claude/confident-gates-0f627e) ONCE
# into a dedicated read-only library, so the A/B array never recompiles.
# Deps (TreeTools/Rcpp/ape/...) are reused from the existing populated $DEPLIB.
# Submit first; chain the array on afterok of this job.
#SBATCH --job-name=ts-naiw-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/naiw_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/naiw_build_%j.err

module load r/4.5.1
module load gcc/14.2

BRANCH=claude/confident-gates-0f627e
REPO=/nobackup/$USER/TreeSearch-t29
DEPLIB=/nobackup/$USER/TreeSearch/lib          # has TreeTools/Rcpp/ape/...
LIB=/nobackup/$USER/TreeSearch/lib-naiw        # fresh target for this build
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs
export R_LIBS_USER="$LIB:$DEPLIB"

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO" || { echo "FATAL: clone failed"; exit 1; }
fi
cd "$REPO" || { echo "FATAL: no $REPO"; exit 1; }
git fetch origin "$BRANCH" && (git checkout "$BRANCH" && git reset --hard "origin/$BRANCH")
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"
exit $rc
