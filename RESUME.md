# TreeSearch — hand-off 2026-05-19

## What this repo is
R package implementing maximum-parsimony phylogenetic tree search. Core search (TBR/SPR/NNI + ratchet/drift/fuse heuristics) implemented in C++ (`src/ts_*.cpp`); R layer provides Shiny UI, inapplicable-character scoring, and implied-weights. Branch `cpp-search` is a draft PR (`cpp-search → main`, PR #210) consolidating C++ search engine improvements.

## Where we left off

- **b67db1a1** `chore(profile)`: Round 2 baseline filed.  VTune confirmed `full_rescore` is 19.2 % self-time on Zhu2013 thorough preset (`fitch_na_score` 0.585 s + `load_tip_states` 0.032 s) — all flowing through `tbr_search:1138` accept path.  Zhu2013 has NA so the initial T-300 commit does not benefit it; an NA-aware variant is the next follow-up.
- **f531bbcd** `perf(T-300)`: Dirty-set incremental rescore for SPR accept.  Replaces the reverted overlap-chain implementation (1e3fc9a7).  Each affected node visited exactly once in postorder, sidestepping the −3 bug analytically.  Gated on `is_spr && !has_na`; DEBUG_RESCORE cross-check enabled.  92 local tests pass with 0 mismatches across EW + IW.
- **c504ea87** `fix(drift)`: Missing `drift_restore_topology(tree, snap)` in the second `drift_apply_tbr_move` failure handler (`src/ts_drift.cpp:680`). The EW RFD re-apply path (line 677) lacked the restore call that the first apply's handler (line 589-593) had correctly. Theoretically unreachable (same args succeed on first apply → succeed on second), but robustness fix applied.
- **b7303ee** `revert(T-300)`: T-300 (incremental SPR rescore) had a consistent `diff=-3` bug. Root cause found (analytically): overlapping ancestor paths in `fitch_incremental_downpass(nz) + fitch_incremental_downpass(nx)` when regraft edge is above nz. Chain 2 subtracts already-updated `above` local_cost (from chain 1) instead of original. Reverted to `full_rescore(tree, ds)` at the accept path in `tbr_search` (`src/ts_tbr.cpp:1138`).
- **b186e801** `fix(progress)`: Guards `R_FlushConsole()` behind `extern Rboolean R_Interactive` in the parallel progress loop (`src/ts_parallel.cpp:435`). Root cause of 6h GHA hangs.
- **44d929a8** `perf(sector)`: `copy flat_blocks` and `all_weight_one` in `build_reduced_dataset` — landed and passing.

## Pending jobs

| Type | ID / ref | Status | ETA | On completion |
|------|----------|--------|-----|---------------|
| GHA  | 26078069315 | in_progress | ~1–2 h | R-CMD-check on T-300 push (b67db1a1).  Verify zero DEBUG_RESCORE assertions across all platforms + that the suite is green. |
| GHA  | 26078069318 | in_progress | ~2–3 h | gcc-ASAN on T-300 push.  Watch for stack-smash / heap UB (prior T-300 attempt smashed after ~13 corrected iters). |
| PR   | #244 (T-302) | open, CI pending | — | Requeue after cpp-search CI green |
| PR   | #242 (T-298) | open | — | Windows DLL crash at `library("TreeSearch")` — GCOV instrumentation failure; needs rebase + diagnosis |

## Open items / next steps

1. **Monitor 26078069315 + 26078069318** — DEBUG_RESCORE assertion must show zero diffs across all platforms.  If clean: file a follow-up commit removing `#define DEBUG_RESCORE` (line 21 of `src/ts_tbr.cpp`).  If any diffs: read the values from the GHA log; the dirty-set algorithm is correct by construction, so a diff would indicate either (a) a node missing from the dirty walk that should be there, or (b) postorder not reflecting post-move topology — verify `build_postorder_prealloc` ran first at line 1136.
2. **NA variant of T-300** — Zhu2013 (the main benchmark) has inapplicable chars and won't benefit until T-300 covers NA.  Mirror `fitch_dirty_downpass` for the three-pass NA scoring.  Expected to recover the 18.2 % of DLL time spent in `fitch_na_score` via `full_rescore`.
3. **PR #244 (T-302)** — once cpp-search CI is green, push a no-op commit to requeue checks.
4. **PR #242 (T-298 quartet-concordance)** — rebase `feature/quartet-concordance-opt` from main; diagnose Windows DLL crash (GCOV instrumentation).
5. **Red-team area #5** (Data pipeline & simplification) — next in rotation after area #4 (Parallelism & RNG, clean).
6. **Profiling area #6** (`any_hit_reduce_avx2` — 9.6 % of DLL time, secondary T-300-independent candidate identified by round 2).

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

- **T-300 overlapping-chain hypothesis (PRIOR ANALYSIS WRONG)**: Earlier sessions claimed the −3 came from chain-2 subtracting chain-1's already-updated `above` local_cost.  Detailed accounting (Opus 4.7, 2026-05-19) shows the chain deltas DO sum correctly: chain1 = wrong − original, chain2 = correct − wrong, sum = correct − original.  The real source of the −3 is unidentified.  The dirty-set algorithm (commit f531bbcd) sidesteps the question by visiting each node exactly once.
- **Pre-hoc constraint screening via `impose_constraint`**: Move rejection is O(1) post-move; pre-screening only needed for NNI perturbation, Fuse, and Drift — not TBR.
- **R_FlushConsole removal entirely**: Would lose real-time progress display in interactive sessions. Guard on R_Interactive is the right approach.

## Suggested first action

```bash
gh run view 26078069315 --json status,conclusion
# Search the log for DEBUG_RESCORE lines (zero expected):
gh run view 26078069315 --log | grep DEBUG_RESCORE
```

If clean: remove `#define DEBUG_RESCORE` at `src/ts_tbr.cpp:21` in a follow-up commit, then start the NA variant of `fitch_dirty_downpass`.
