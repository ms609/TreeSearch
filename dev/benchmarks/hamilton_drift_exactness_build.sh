#!/bin/bash
# Build TreeSearch from the drift-exactness branch into a DEDICATED shared lib
# (does NOT clobber the cpp-search shared lib other sessions use). Chain the
# gate array on afterok of this job.
#SBATCH --job-name=ts-driftexact-build
#SBATCH -p shared
#SBATCH -n 4
#SBATCH --mem=8G
#SBATCH --time=0:45:00
#SBATCH --output=/nobackup/%u/TreeSearch/logs/driftexact_build_%j.out
#SBATCH --error=/nobackup/%u/TreeSearch/logs/driftexact_build_%j.err

set -e
module load r/4.5.1
module load gcc/14.2

BRANCH=claude/happy-jennings-e7a165
LIB=/nobackup/$USER/TreeSearch/lib-driftexact
REPO=/nobackup/$USER/TreeSearch-driftexact
mkdir -p "$LIB" /nobackup/$USER/TreeSearch/logs

if [ ! -d "$REPO/.git" ]; then
  git clone https://github.com/ms609/TreeSearch.git "$REPO"
fi
cd "$REPO"
git fetch origin "$BRANCH"
git checkout "$BRANCH"
git reset --hard "origin/$BRANCH"
echo "Git HEAD: $(git log --oneline -1)"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" TreeSearch_*.tar.gz
rc=$?
rm -f TreeSearch_*.tar.gz
echo "INSTALL exit: $rc; version: $(R_LIBS_USER=$LIB Rscript -e 'cat(as.character(packageVersion("TreeSearch")))' 2>/dev/null)"
exit $rc
