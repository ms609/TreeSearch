# Agent A Progress Log

## Current Task
**IDLE** — S-PROF round 6 complete, T-204 fix dispatched.

## S-PROF Round 6 + T-204 fix + T-266 tidy (2026-03-27 ~10:50 GMT)

### T-266 tidy-up
- Deleted from to-do.md; added to completed-tasks.md (2026-03-27 section).
- Removed feature/prune-reinsert branch (local + origin).
- TS-PruneRI directory orphaned (git metadata already gone); manual cleanup needed.

### T-204 fix (GHA 23641482723 → 23643078732)
- GHA 23641482723 failed: spelling ERROR ('cleanup'/'phyDat') + deprecated-fn
  warnings in examples. Root: T-204 added `.Deprecated()` to `PhyDat2Morphy`/
  `UnloadMorphy`; examples calling those (PhyDat2Morphy.Rd, MorphyWeights.Rd,
  GapHandler.Rd, SingleCharMorphy.Rd, Morphy.R constraint example) now emit
  warnings in `R CMD check`.
- Fix (eb21c588 on feature/native-search): WORDLIST + `\donttest{}`/
  `suppressWarnings()` wrappers. Re-dispatched as GHA 23643078732.

### S-PROF Round 6: thorough-preset phase distribution (Zhu2013, 75t)
Built cpp-search HEAD (post T-261+T-262+T-263). Ran MaximizeParsimony with
verbosity=2 to capture phase timings for thorough preset.

**Phase distribution (3 reps, ~11.2s/rep):**
| Ratchet 46.3% | NNI-perturb 34.3% | RSS 7.4% | CSS 4.4% | XSS 3.2% | TBR 3.2% |

**Key finding:** NNI-perturb = 34% of time with 14% hit rate (1 step/hit).
TBR is negligible (3%), confirming T-261+T-262+T-263 effectiveness.
Filed T-274 (P2): benchmark nniPerturbCycles=0 vs 5 at thorough-preset scale.

## T-204 + T-266 fixes (2026-03-27 ~10:15 GMT)
- T-204 (PR #216): GHA 23495097795 was a timing issue — docs commit `f59a193c` landed after the run was dispatched. Current HEAD (11622e90) has correct Rd files. Re-dispatched as 23641482723.
- T-266 (PR #235): Standard CI (R-CMD-check + gcc-ASAN) failed after PR opened. R CMD check failure: spelling ERROR — 'warmup' (from T-270 vignette) and 'config' not in WORDLIST for R 4.1 hunspell. Fixed in `de9e5210` (TS-PruneRI). Re-dispatched agent-check as 23641870390. gcc-ASAN/devel failures are infrastructure (rlang compile error), not package issues.

## S-PR + S-RED focus 1 (2026-03-27 ~10:00 GMT)
- S-PR: Updated T-204 to-do entry with GHA run 23495097795 failure details — undocumented `CleanNativeData`/`NativeBootstrap`/`NativeLength`/`PrepareNativeData`, codoc mismatches in `Jackknife.Rd`/`Ratchet.Rd`/`TreeSearch.Rd`. B needs to regenerate Rd files and add roxygen2 docs.
- S-RED focus 1: Reviewed ts_fitch.h/.cpp, ts_fitch_na.h, ts_fitch_na_incr.h, ts_simd.h. Focus on commits since 2026-03-19 (AVX2 dispatch, FlatBlock flat indirect, XFORM integration). No bugs found. AVX2 ops bit-identical to scalar; flat functions infrastructure only; XFORM no double-count (weight=0 removes hierarchy chars from Fitch blocks); incremental downpass/uppass stopping conditions correct.

## S-RED Focus 10 + S-PR + T-270 (2026-03-27 ~09:30–10:00 GMT)
- S-RED focus 10: reviewed ts_fitch.cpp IW/Profile paths. BUG FIXED: precompute_profile_delta old_cost=0 when s>info_max_steps. 15 tests pass. commit 7cff7870.
- S-PR: merged cpp-search into #235 (prune-reinsert, clean) and #216 (native-search, clean). #213 (cid-consensus) has ts_tbr.cpp conflict (CID vs T-263 snapshot) — aborted, needs E/human. Closed stale PR #178 (T-272 done).

## S-COORD Round 31 + T-270 (2026-03-27 ~09:20 GMT)
- T-266 PR #235 opened (GHA passed).
- T-150 GHA 23636944848 FAILED — InfoConsensus.Rd codoc mismatch. Updated T-150 row in to-do.md.
- Filed T-270 (vignette docs for T-257), T-272 (close PR #178).
- Completed T-270: updated vignettes/search-algorithm.Rmd (new pipeline step 5a, post-ratchet sectorial subsection, fixed stale consensusStableReps docs); updated AGENTS.md pipeline. commit d8f3c769.
- u.005.claimed-F: skipped (claimed by F).


### Session: 2026-03-27

Implemented taxon pruning-reinsertion (T-266): a perturbation strategy that
drops ~10% of leaves, TBR-optimizes the reduced backbone, then greedily
re-adds the dropped taxa via Wagner insertion + TBR polish. Complements the
ratchet (weight-space) and NNI-perturbation (topology-space).

**Commit:** `afbf531f` on `feature/prune-reinsert`

**Files added:**
- `src/ts_prune_reinsert.h/.cpp` — core algorithm (random + instability-weighted tip selection)
- `tests/testthat/test-ts-prune-reinsert.R` — 44 assertions (Tier 2)

**Files modified:**
- `src/ts_driven.h/.cpp` — pipeline phase 5c, timing, outer-cycle division
- `src/ts_wagner.h/.cpp` — exposed 3 helpers for reuse
- `src/ts_rcpp.cpp` — param unpacking + timing output
- `R/SearchControl.R` — 3 new params (pruneReinsertCycles/Drop/Selection)
- `R/ts-driven-compat.R` — backward-compat wrapper

**Local validation:** Build clean, 44/44 prune-reinsert tests pass,
234 related tests (driven/nni-perturb/wagner) pass with no regressions.

**GHA runs:**
- Run 23634563604: FAIL — `INT_MAX` undeclared on Linux/ARM (missing `<climits>`)
- Run 23635469688: FAIL — Codoc mismatch (SearchControl.Rd not regenerated)
- Run 23636145497: PASS — PR #235 opened to cpp-search

---

## Session: 2026-03-26 — S-RED focus 9 review

### Completed: S-RED standing task — focus area 9 (Wagner & addition trees)

Reviewed `ts_wagner.h/.cpp` (595 lines) and `ts_constraint.h/.cpp` (736+144 lines).

**Key findings:**
- No bugs found in Wagner tree construction (incremental scoring, constraint mapping, 3-taxon base case, biased addition, random constrained tree)
- Latent stale-reference issue in `impose_one_pass()` (best_node relocated when move_out_root is a direct child) — negligible severity, mitigated by retry loops and TBR enforcement
- `regraft_violates_constraint()` DFS timestamp logic verified correct
- `classify_clip_constraints()` bit masking and FORBIDDEN classification correct
- 902 constraint-related tests pass; 80/80 adversarial tests pass

No new bugs filed.

### Earlier: T-242 investigation (closed)
Confirmed T-242 was a display bug (ThreadSafePool::extract_into() resetting hits_to_best), not a search quality regression. Actual IW hit rate ~60-67%.

## Session: 2026-03-25 (evening) — Summary

### Completed: T-208 + T-211 → PR #229

Implemented `random_constrained_tree()` and fixed three `impose_constraint()`
bugs on `feature/random-constrained-tree` (worktree `TS-RCT`).

**GHA run 23557186264:** 0 FAIL, 10927 PASS on both Ubuntu and Windows.

**PR #229** created to cpp-search.
