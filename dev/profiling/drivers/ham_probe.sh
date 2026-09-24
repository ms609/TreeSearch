#!/bin/bash
echo "===REPOS==="
for d in /nobackup/$USER/TreeSearch /nobackup/$USER/TreeSearch-a /nobackup/$USER/TreeSearch-t29; do
  if [ -d "$d/.git" ]; then
    echo "$d: $(git -C "$d" log --oneline -1 2>/dev/null) [$(git -C "$d" rev-parse --abbrev-ref HEAD 2>/dev/null)]"
  else
    echo "$d: (no .git or absent)"
  fi
done
echo "===LIB DESCRIPTION==="
desc=/nobackup/$USER/TreeSearch/lib/TreeSearch/DESCRIPTION
if [ -f "$desc" ]; then grep -E "^Version|^Packaged|^Built" "$desc"; else echo "no installed TreeSearch at lib"; fi
echo "===TNT 64bit==="
file /nobackup/$USER/TreeSearch/tnt/TNT-bin/tnt 2>/dev/null
echo "===DATASETS in t29 repo==="
ls /nobackup/$USER/TreeSearch-t29/dev/benchmarks/bench_tnt_headtohead.R 2>/dev/null && echo "harness present"
