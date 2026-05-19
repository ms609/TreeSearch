# TreeSearch — hand-off 2026-05-19

## What this repo is
R package implementing maximum-parsimony phylogenetic tree search. Core search (TBR/SPR/NNI + ratchet/drift/fuse heuristics) implemented in C++ (`src/ts_*.cpp`); R layer provides Shiny UI, inapplicable-character scoring, and implied-weights. Branch `cpp-search` is a draft PR (`cpp-search → main`, PR #210) consolidating C++ search engine improvements.

## Where we left off

- **c504ea87** `fix(drift)`: Missing `drift_restore_topology(tree, snap)` in the second `drift_apply_tbr_move` failure handler (`src/ts_drift.cpp:680`). The EW RFD re-apply path (line 677) lacked the restore call that the first apply's handler (line 589-593) had correctly. Theoretically unreachable (same args succeed on first apply → succeed on second), but robustness fix applied.
- **b7303ee** `revert(T-300)`: T-300 (incremental SPR rescore) had a consistent `diff=-3` bug. Root cause found (analytically): overlapping ancestor paths in `fitch_incremental_downpass(nz) + fitch_incremental_downpass(nx)` when regraft edge is above nz. Chain 2 subtracts already-updated `above` local_cost (from chain 1) instead of original. Reverted to `full_rescore(tree, ds)` at the accept path in `tbr_search` (`src/ts_tbr.cpp:1138`).
- **b186e801** `fix(progress)`: Guards `R_FlushConsole()` behind `extern Rboolean R_Interactive` in the parallel progress loop (`src/ts_parallel.cpp:435`). Root cause of 6h GHA hangs.
- **44d929a8** `perf(sector)`: `copy flat_blocks` and `all_weight_one` in `build_reduced_dataset` — landed and passing.

## Pending jobs

| Type | ID / ref | Status | ETA | On completion |
|------|----------|--------|-----|---------------|
| GHA  | 26075916824 | in_progress | ~2-3h total | ubuntu-release + macos-intel + windows PASSED. arm/4.1/devel/macOS-latest still at "Check package" step (started 04:18 UTC). If all pass: R_Interactive fix confirmed; cherry-pick/merge PR #244 then proceed to T-300 redesign. If arm/devel hang again: replace `R_Interactive` with `isatty(fileno(stdout))`. |
| GHA  | 26075916857 | in_progress | ~2-3h total | gcc-ASAN for same push. Check after 26075916824 completes. |
| PR   | #244 (T-302) | open, CI pending | — | Requeue after cpp-search CI green (push no-op commit) |
| PR   | #242 (T-298) | open | — | Windows DLL crash at `library("TreeSearch")` — GCOV instrumentation failure; needs rebase + diagnosis |

## Open items / next steps

1. **Monitor 26075916824** — arm/4.1/devel/macOS-latest still in progress. If all pass: R_Interactive fix confirmed. If arm/devel hang again: replace with `isatty(fileno(stdout))` (`<unistd.h>` / `<io.h>`).
2. **PR #244 (T-302)** — once cpp-search CI is green, push a no-op commit to requeue checks.
3. **T-300 redesign** — root cause: overlapping ancestor paths. Fix: accumulate dirty_set = union(nz→root, nx→root), walk postorder once updating each node exactly once. Start at `fitch_incremental_downpass` at `src/ts_fitch.cpp:134`.
4. **PR #242 (T-298 quartet-concordance)** — rebase `feature/quartet-concordance-opt` from main; diagnose Windows DLL crash (GCOV instrumentation).
5. **Red-team area #5** (Data pipeline & simplification) — next in rotation after area #4 (Parallelism & RNG, clean).
6. **Profiling round 1** — scaffold done (`dev/profiling/`); area #1 or #2 is next target; needs VTune on Windows. VTune updated to 2026.0 at `C:\Program Files (x86)\Intel\oneAPI\vtune\2026.0`.

## Technical pointers

- PR #210: cpp-search → main (draft)
- Main C++ sources: `src/ts_tbr.cpp`, `src/ts_parallel.cpp`, `src/ts_fitch.cpp`, `src/ts_tree.cpp`
- Progress reporting + R_Interactive guard: `src/ts_parallel.cpp` lines 424–455
- TBR accept path: `src/ts_tbr.cpp:1138` (`double actual = full_rescore(tree, ds)`)
- Incremental rescore functions: `fitch_incremental_downpass` at `src/ts_fitch.cpp:134`, `fitch_incremental_uppass` at `src/ts_fitch.cpp:199`
- PreallocUndo: `src/ts_tree.h:86`; capacity = `3 * n_node`; `grow()` at overflow
- Drift EW RFD re-apply path: `src/ts_drift.cpp` lines 677-685 (second apply, failure handler)
- Profiling scaffold: `dev/profiling/focus-areas.md`, `dev/profiling/log.md`
- Red-team log: `.positai/expertise/red-team.md` (gitignored), current `last_focus: area 4`
- Shiny server: `inst/Parsimony/server.R`, search module `inst/Parsimony/server/mod_search.R`
- VTune 2026.0: `C:\Program Files (x86)\Intel\oneAPI\vtune\2026.0`

## Things ruled out

- **T-300 incremental chain1+chain2**: Root cause FOUND: When regraft edge (above, below) is ABOVE nz in the tree, chain 1 `fitch_incremental_downpass(nz)` walks nz→above→root, UPDATES prelim/local_cost at `above`. Chain 2 `fitch_incremental_downpass(nx)` then subtracts the already-UPDATED `above` local_cost rather than the original — giving a delta relative to chain-1 output. The -3 is the accumulated double-subtraction error over the shared ancestor path. FIX: dirty_set = union(nz→root, nx→root), walk postorder once.
- **Pre-hoc constraint screening via `impose_constraint`**: Move rejection is O(1) post-move; pre-screening via constraint repair is only needed for NNI perturbation, Fuse, and Drift — not TBR.
- **R_FlushConsole removal entirely**: Would lose real-time progress display in interactive sessions. Guard on R_Interactive is the right approach.
- **isatty() as alternative to R_Interactive**: Valid fallback if R_Interactive is not accessible on some runners, but R_Interactive is the correct R-native test. Use only if arm/devel hang again.

## Suggested first action

```r
# Check if the remaining GHA jobs completed:
gh run view 26075916824 --json status,conclusion,jobs | \
  python -c "import json,sys; d=json.load(sys.stdin); [print(j['name'], j['status'], j.get('conclusion','')) for j in d['jobs']]"
```

If arm/4.1/devel/macOS-latest all succeeded: cherry-pick the T-302 fix from PR #244 into cpp-search.
If any hung again: apply the `isatty(fileno(stdout))` fallback to `src/ts_parallel.cpp:14,435`.
