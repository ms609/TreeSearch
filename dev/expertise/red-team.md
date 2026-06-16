# Red-Team Expertise — TreeSearch (durable wisdom)

> **This file is the curated "lessons learned" — not the operational machinery.**
> The live rotation now lives under `dev/red-team/` (per the `/red-team` skill):
> - `dev/red-team/focus-areas.md` — the 10-area rotation table + per-area `start_tier`
> - `dev/red-team/log.md` — append-only round log + `last_focus:` pointer
> - `dev/red-team/findings.md` — verified OPEN findings
>
> Run a round with `/red-team`. Keep this file for the *patterns* that recur across
> rounds; keep the *records* in `dev/red-team/`. (Migrated 2026-06-16 from
> `.positai/expertise/red-team.md`, which is now a stub.)

## Purpose

Red-teaming reviews code for (i) correctness bugs and (ii) performance issues. Fix
trivial issues inline (and note them in the round log); file non-trivial, *verified*
findings in `dev/red-team/findings.md` (and the dispatcher queue `to-do.md`).

The goal is **issues fixed per token spent**, not issues found in the abstract. Depth
over breadth: one focused review that finds a real bug beats a broad "all green" sweep.

## Bug patterns (reference — distilled from past rounds)

| Pattern | Where to check | Seen in |
|---------|---------------|---------|
| Missing `GetRNGstate()`/`PutRNGstate()` around `unif_rand()` | Any `.cpp` using randomness | — |
| **R RNG API (`unif_rand`/`Get/PutRNGstate`) reached from a worker thread** | Serial helpers (`driven_search`, `resample_search`) called *from* a parallel worker — the worker path may re-enter a serial entry point that brackets R's RNG | T-309 (parallel `Resample()` → `ts_driven.cpp:690` via `ts_resample.cpp`) |
| `std::random_device{}()` ignoring `set.seed()` | Seeding of `std::mt19937`; use `ts::make_rng()` | — |
| `1u << s` undefined behaviour at `s >= 32` | `build_dataset`/`simplify_patterns` token bitmasks (`uint32_t`) when `n_states == 32` | DAT-001 (UBSAN); B's 2026-03-19 `MAX_STATES` guard |
| Division by an unguarded count → `NaN`/`Inf` that defeats a later guard | XPIWE `f = 1 + r*missing/obs` when `obs == 0`; NaN passes `> 0` and `<= floor` comparisons silently | DAT-002; T-311 (LS RSS=0 garbage) |
| `NaN`/`Inf` clamped to a benign value, hiding failure | `return rss > 0 ? rss : 0` turns NaN into a "valid" 0; validate finiteness instead | T-311 (LeastSquares) |
| `build_reduced_dataset` doesn't copy a needed field | Sector / prune-reinsert reduced datasets: `flat_blocks`, `all_weight_one`, `inapp_state`, HSJ (`hierarchy_blocks`/`tip_labels`) & Sankoff fields. Two copies of this function — keep them in sync | T-275, T-303; B 2026-03-19; 44d929a8 |
| Frozen-API return type drift (`as.integer` vs `as.logical`) | `R/SearchControl.R` boolean fields documented as logical | T-310 (`pruneReinsertNni`) |
| GCC-only builtins (`__builtin_popcountll`, etc.) | All `.cpp`/`.h` (need MSVC fallback) | — |
| `.inc`/`.h` changes not triggering recompilation | `ts_fitch_na.inc`, `ts_fitch_na_incr.inc`; `touch src/ts_fitch.cpp` after edits; no `Makevars` header-dep tracking | — |
| Missing `TreeSearch:::` prefix / no `set.seed()` before `sample()` | `tests/testthat/test-ts-*.R` | 2026-05-26 area 8 |
| Arg-count mismatch in `TreeSearch-init.c`; `R_PosInf` in Rcpp defaults | After adding/removing Rcpp params; `R/RcppExports.R` after `compileAttributes()` | — |
| Stale state arrays after a rejected move | Restore only covers the clip-to-root path; regraft-to-root nodes keep regrafted values — conservative (screening only) but real | T-235 (SPR); A 2026-03-19 |
| Stale metadata after **equal-score / tabu** rejection | Constraint `map_constraint_nodes`/`compute_dfs_timestamps` re-synced on the `!accepted` path but not the tabu path | T-316 (possible P1) |
| Out-of-range index from R reaching the kernel as `NA_INTEGER` (INT_MIN) | `AdditionTree(sequence=)` → `addition_order[i]-1` underflows → OOB write; validate in R | WGN-01 (P1) |
| No revert on worsening move; mask not cleaned up | Sectorial reinsertion, fuse exchange; ratchet perturbation restore paths | B 2026-03-19 |
| `pattern_freq *= 2` per char → exponential blowup when patterns shared | Ratchet `perturb_upweight`/`perturb_mixed`; use `+= 1` | B 2026-03-19 |
| `-Inf − (-Inf) = NaN` in log-space convolution | `.LogCumSumExp` (`R/pp_info_extra_step.r`) when both terms `-Inf` | B 2026-03-20 |

### Shiny (area 7) — async lifecycle patterns (immature seam, 2026-06-16)
- **Stamp the hash of the dataset that was *prepared*, not the current one**, when an
  async task completes — otherwise a mid-task data swap scores against the wrong data
  (T-309).
- **Re-entrancy:** `shinyjs::disable()` is an async browser round-trip; guard handlers
  with an explicit `searchInProgress` flag, not the disabled button (T-310).
- **`onStop` must write the worker's cancel signal** and `unlink` *all* temp-file
  prefixes the module actually creates (`ts_cancel_*`, `ts_progress_*`, `ts_profile_*`),
  or workers orphan and `tempdir()` grows (T-311, T-312; see `shiny-app.md` "Issue 6").
- **Topology dedup must strip branch lengths** before `write.tree()` (it serialises BLs),
  or BL-bearing user trees won't dedupe against parsimony trees (T-313).

## Performance patterns (reference)

| Pattern | Where to check |
|---------|---------------|
| Unbounded indirect scoring (missing `_bounded` variant) | Search inner loops |
| Full `score_tree()` where incremental / dirty-set would suffice | After clip/regraft; NNI accept calls full `fitch_uppass` (O(n)) |
| `build_postorder()` called unnecessarily | After unclip or snapshot restore |
| Full-tree `TreeState` copy where save/restore (`copy_topology`) suffices | Fuse, sectorial, NNI-perturb |
| `SankoffData` rebuilt every `score_tree()` | XFORM scoring hot path (could cache in DataSet) |
| Long mutex hold across `tree_fuse()` | `ThreadSafePool::fuse_round` (correct, but serialises) |

## Known fragile areas (reference)

1. `ts_rcpp.cpp` + `TreeSearch-init.c`: append-only, check arg counts every time.
2. `RcppExports.R/.cpp`: concavity `Inf` → `-1.0`/`HUGE_VAL` sentinel after regen.
3. `.inc`/`.h` files: `touch src/ts_fitch.cpp` after changes (no header-dep tracking).
4. Parallel: `ts_rng.h` `thread_local` must be set before any search call; **no R API
   from workers** — beware serial entry points re-entered from a worker.
5. `init_from_edge`: first child → left convention; 1-based R edge → 0-based internal.
6. **Two `build_reduced_dataset` functions** (`ts_sector.cpp`, `ts_prune_reinsert.cpp`) —
   asymmetric field-copy footgun; keep in sync, guard HSJ/XFORM where fields are absent.
7. `FlatBlock.active_mask` duplicates `CharBlock.active_mask` — any writer to one must
   sync the other (ratchet does; future writers must too).
8. Incremental/dirty-set Fitch rescore (`ts_fitch_na_dirty.h`, `ts_tbr.cpp:1138-1180`):
   the crown-jewel correctness risk (T-300 was a systematic delta=−3). It has **no
   enduring regression test** beyond T-304 — treat with suspicion after any edit.

## Reporting format (for `dev/red-team/findings.md`)

```
| T-NNN | P1/P2/P3 | <area #> | **Title.** | `path:line` | Detail + fix + verifier verdict. |
```

Severity: **P1** wrong user-visible result / crash / desk-reject · **P2** wrong on edge
input / frozen-shape inconsistency / search-quality · **P3** robustness / polish.
