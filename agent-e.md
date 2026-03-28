# Agent E — Progress Log

## Current Task
- **Status:** PARKED — GHA 23690338955 (feature/tbr-batch); Hamilton down

### T-289f — PR NNI polish cost reduction (2026-03-28)

Root cause of Stage 4 failure identified: full TBR convergence on the full
tree after every PR cycle (~7s/cycle × 5 = ~35s, before outer TBR runs again).

Added two new SearchControl() params:
- `pruneReinsertNni = TRUE`: NNI instead of TBR for full-tree polish (~5x cheaper)
- `pruneReinsertFullMoves = N`: limit full-tree TBR moves (0 = converge, backward compat)

7 files changed: ts_prune_reinsert.h/.cpp, ts_driven.h/.cpp, ts_rcpp.cpp,
R/SearchControl.R, man/SearchControl.Rd. commit 09c93468.

Stage 5 benchmark script created (bench_pr_stage5_nni.R + t289f_stage5_hamilton.sh).
3 configs × 5 datasets × 2 budgets × 10 seeds = 300 runs.
Committed aa3f16ea. Hamilton unreachable; submit when back:
  sbatch /nobackup/pjjg18/TreeSearch-a/dev/benchmarks/t289f_stage5_hamilton.sh

S-RED on new NNI branch: clean. No bugs. Constraint-staleness non-issue
(tbr_search re-syncs cd at entry). Timeout handling correct.

### T-289 COMPLETE (2026-03-28)

Stage 4 (multi-dataset validation) results: PR adds ~90% per-replicate overhead
at 206 tips. syab07205/206t: 0 replicates at 60s budget. pruneReinsertCycles=0L
in large preset. commit 74698524.

### Codoc fix — SearchControl.Rd (E-003, 2026-03-28)
Rd missing pruneReinsertTbrMoves. Fixed manually. commit fdf25673.
GHA 23687210711 PASSED.

### T-291 — bench_framework.R interface (E-004)
benchmark_run() updated to three structured lists. commit f1ed5dfc.

### S-RED E-005 — ts_strategy.h + ts_temper (no bugs)
### S-RED E-002 — ts_rng + ts_parallel (no bugs)
