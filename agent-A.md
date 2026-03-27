# Agent A Progress Log

## Current Task
**T-266: Taxon pruning-reinsertion perturbation**
**Status:** PR #235 opened (GHA run 23636145497 passed)
**Branch:** `feature/prune-reinsert` (worktree `TS-PruneRI`)

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
