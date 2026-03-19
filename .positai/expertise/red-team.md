# Red-Team Expertise — TreeSearch

## Purpose

Red-teaming reviews code for (i) bugs and (ii) performance issues.
Fix trivial issues directly; add non-trivial ones to `to-do.md`.

## Focused rotation system

Each S-RED invocation targets **one focus area** from the rotation below.
The agent reads `last_focus` (bottom of this file) to determine which area
was reviewed last, then picks the **next** area in sequence. After
completing the review, update `last_focus`.

### Focus areas

| # | Area | Scope | Key questions |
|---|------|-------|---------------|
| 1 | **Fitch scoring correctness** | `ts_fitch.h/.cpp`, `ts_fitch_na.h`, `ts_fitch_na_incr.h` | Does incremental scoring match full `score_tree()`? Bounded variants bail correctly? NA three-pass edge cases? Write a targeted test if you find a gap. |
| 2 | **Search topology invariants** | `ts_tbr.cpp`, `ts_drift.cpp`, `ts_search.cpp` | After every rejected move, is topology fully restored? Undo stack correct? No stale `postorder`? Symmetry-breaking hash collisions? |
| 3 | **Ratchet & perturbation** | `ts_ratchet.cpp`, `ts_sector.cpp`, `ts_fuse.cpp` | `active_mask`/`upweight_mask` fully restored after perturbation? Sectorial reinsertion reverts on worse score? Fuse exchange handles tied scores? |
| 4 | **Parallelism & RNG** | `ts_parallel.cpp`, `ts_rng.h/.cpp`, `ts_driven.cpp` | Thread-local RNG set before any search call? No R API calls from worker threads? Pool mutex correct? Atomic stop flag races? Seeds generated from R RNG before spawning? |
| 5 | **Data pipeline & simplification** | `ts_data.h/.cpp`, `ts_simplify.h/.cpp`, `ts_constraint.h/.cpp` | `build_dataset` handles edge cases (all-ambiguous, single-state, zero-weight)? `build_reduced_dataset` copies all fields? Constraint column-major indexing correct? |
| 6 | **R ↔ C++ interface** | `ts_rcpp.cpp`, `TreeSearch-init.c`, `R/RcppExports.R`, `R/MaximizeParsimony.R` | Arg counts match? Concavity sentinel translated correctly? Edge matrix conventions? Return value attributes set? Parameter validation in R layer? |
| 7 | **Shiny module wiring** | `inst/Parsimony/server.R`, `server/mod_*.R`, `server/events.R` | Forward-ref callbacks resolve? Cross-module `updateXxxInput` targets correct namespace? Reactive graph has no orphaned observers? `isolate()` used correctly in result observers? |
| 8 | **Test suite health** | `tests/testthat/test-ts-*.R`, `tests/testthat/helper-ts.R` | Tier guards correct? Any tests that always pass (vacuous)? Missing `TreeSearch:::` prefixes? Edge-case coverage gaps (3-tip, single-char, all-NA)? Flaky tests? |
| 9 | **Wagner & addition trees** | `ts_wagner.h/.cpp` | NA-incremental scoring staleness acceptable? Constraint mapping (LCA-based) correct? Retry loop fires when needed? 3-taxon base case handles all orderings? |
| 10 | **Profile & IW scoring** | `ts_fitch.cpp` (IW/profile paths), `ts_data.cpp` (precompute) | `e/(k+e)` delta correct? Profile `info_amounts` lookup matches? `concavity = 1.0` sentinel activates weighted path? `precompute_profile_delta` includes `precomputed_steps` offset? |

### Workflow for a focused review

1. Read `last_focus` below → pick next area in rotation.
2. **Read the target files** thoroughly (not just skimming). Understand the
   logic before looking for bugs.
3. **Construct a specific adversarial scenario** — e.g. "what happens if
   the ratchet upweights a character that's already at max weight?" — and
   trace the code path.
4. **Write or run a targeted test** that exercises the scenario. If the
   test passes, note it. If it fails, file a bug.
5. Check for any **recent commits** touching the focus files (`git log`).
6. **Build and run** at minimum the test files most relevant to the focus
   area (not necessarily the entire suite).
7. Update `last_focus`, report findings.

Total time budget: aim for depth over breadth. A focused review that finds
one real bug is worth more than a broad sweep that confirms "all green."

## Bug patterns (reference — from past rounds)

| Pattern | Where to check |
|---------|---------------|
| Missing `GetRNGstate()`/`PutRNGstate()` around `unif_rand()` | Any .cpp using randomness |
| `std::random_device{}()` ignoring `set.seed()` | Seeding of `std::mt19937` |
| GCC-only builtins (`__builtin_popcountll`, etc.) | All .cpp/.h files |
| `.inc` file changes not triggering recompilation | `ts_fitch_na.inc`, `ts_fitch_na_incr.inc` |
| Missing `TreeSearch:::` prefix in tests | `tests/testthat/test-ts-*.R` |
| Arg count mismatch in `TreeSearch-init.c` | After adding/removing Rcpp params |
| `R_PosInf` in Rcpp defaults | `R/RcppExports.R` after `compileAttributes()` |
| No revert on worsening move | Sectorial reinsertion, fuse exchange |
| `active_mask`/`upweight_mask` not cleaned up | Ratchet perturbation restore paths |

## Performance patterns (reference)

| Pattern | Where to check |
|---------|---------------|
| Unbounded indirect scoring (missing `_bounded` variant) | Search inner loops |
| Full `score_tree()` where incremental would suffice | After clip/regraft |
| `build_postorder()` called unnecessarily | After unclip or snapshot restore |
| Full-tree copy where save/restore suffices | Fuse, sectorial |
| Missing early termination in loops | Block iteration in scoring |

## Known fragile areas (reference)

1. `ts_rcpp.cpp` + `TreeSearch-init.c`: append-only, check arg counts.
2. `RcppExports.R/.cpp`: concavity `Inf` → `-1.0` sentinel after regen.
3. `.inc` files: `touch src/ts_fitch.cpp` after changes.
4. Parallel: `ts_rng.h` thread_local must be set before search.
5. `init_from_edge`: first child → left convention.

## Reporting format

```
| T-NNN | P1/P2 | OPEN | — | [Bug/Perf] Brief description | Found by S-RED focus #N. File:line. Details. |
```

---

## last_focus

area: 10
reviewed_by: A
date: 2026-03-19
notes: Profile & IW scoring. Thorough review of ts_fitch.cpp (IW/profile paths), ts_data.cpp (build pipeline, precompute). Findings: (1) BUG FIXED — `precompute_iw_delta` overestimated delta when `divided_steps[p] < min_steps[p]` (clamp to e=0 treated below-minimum as at-minimum). Reconnection to min or below costs 0 in IW; old code computed freq/(k+1). Fixed with early return of 0.0 when raw e < 0. Conservative bug: caused missed improving moves in TBR screening, never affected final scores (full_rescore authoritative). All 1397 ts-* tests pass. (2) Concavity sentinel (1.0 for profile) correctly activates weighted pipeline in all 5 search modules (TBR, drift, ratchet, SPR, rcpp). (3) Profile `precomputed_steps` offset applied correctly in both `compute_profile` and `precompute_profile_delta`. Profile delta NOT affected by below-minimum issue because `info_amounts[min_steps, p] = 0` by construction. (4) Known design: `extract_divided_steps` uses stale `local_cost` for NA blocks (always zeros) — documented as intentional heuristic (ts_tbr.cpp lines 39-41). (5) `compute_iw`/`compute_weighted_score`/`precompute_weighted_delta` dispatch chain correct. (6) Ratchet IW perturbation: `pattern_freq` doubled per unique pattern (each pattern_index is unique across all blocks) — correct. (7) EW ew_offset correctly excluded from IW/profile paths; per-pattern precomputed_steps used instead.
