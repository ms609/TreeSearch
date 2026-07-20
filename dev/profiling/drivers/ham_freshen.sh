#!/bin/bash
set -e
REPO=/nobackup/$USER/TreeSearch-t29
cd "$REPO"
echo "fetching origin/cpp-search ..."
git fetch origin cpp-search 2>&1 | tail -3
# Working tree is clean (verified); create/switch local cpp-search tracking origin.
git checkout cpp-search 2>&1 | tail -2 || git checkout -b cpp-search origin/cpp-search 2>&1 | tail -2
git merge --ff-only origin/cpp-search 2>&1 | tail -2
echo "HEAD: $(git log --oneline -1)"
echo "contains 25e35be7 (unrooted-default): $(git merge-base --is-ancestor 25e35be7 HEAD 2>/dev/null && echo YES || echo NO)"
echo "contains 2b299e4b (Wagner fix): $(git merge-base --is-ancestor 2b299e4b HEAD 2>/dev/null && echo YES || echo NO)"
git log --oneline -6
