# Agent E — Progress Log

## Current Task
- **Status:** Working on ASAN vector OOB bug investigation.

### ASAN vector OOB fix — DONE (2026-03-26)
- Root cause: `total_words == 0` when all characters are parsimony-uninformative
  (each state has <2 tips). Simplification removes all patterns, leaving
  `prelim`/`final_`/etc vectors empty. `compute_collapsed_flags`,
  `save_node_state`, and other code accessed `vector[0]` on empty vectors.
- Fix: early returns in TBR, SPR, NNI, drift, ratchet, and collapsed-flags
  when `total_words == 0`.
- Test fix: updated test-ts-collapsed.R datasets to have >=2 tips per state.
- Commit: 6505803f on cpp-search.
- Verified locally with `_GLIBCXX_ASSERTIONS` — 73 tests pass (char-ordering + collapsed).
- Bug was pre-existing on cpp-search, NOT specific to AVX2 branch (PR #233).

### S-COORD Round 27 — DONE
- Fixed R 4.1 `%||%` compat bug in `test-ts-anneal.R` (58fc2552).

### T-265 — RESOLVED (scoring method confound)
### Previous: S-RED Focus 8, T-261+T-262, T-255, T-260
