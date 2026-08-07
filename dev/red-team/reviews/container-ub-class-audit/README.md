# Re-audit: degenerate-container undefined behaviour

Ordered 2026-08-07, after the July sign-off on this class proved wrong twice
in a week (`agent-issues/TreeSearch#124`, `#151`). Both misses shared a
cause: the earlier pass reasoned about reachability *transitively* ("callers
cannot produce zero words") instead of checking, and one of its false claims
was written into memory, where it then actively hid `#124`.

## The class, and why it takes two detectors

| Sub-class | Example | Detected by | Blind to it |
|---|---|---|---|
| 1. Null pointer to a `nonnull` parameter | `memcpy(v.data(), w.data(), 0)` where `v` is empty | UBSan `nonnull-attribute` (gcc-ASAN leg) | `_GLIBCXX_ASSERTIONS`; plain builds |
| 2. Out-of-range container address *formation* | `&v[0]`, `v.back()` on an empty `v` — no load or store | `_GLIBCXX_ASSERTIONS` (`glibcxx-assertions` leg, added by `#60`/PR #133) | ASan, which watches accesses, not address arithmetic |

Neither detector sees the other's sub-class, so any statement of the form
"the sanitizers are clean" is meaningful only once both legs have executed
the shape in question. That is the whole finding: coverage here is bounded
by *input* coverage, not by code reading.

## What was done

**Static pass, sub-class 1 — all 69 `memcpy`/`memmove`/`memset` sites.**
Classified per site as *guarded* / *provably non-empty at the site* /
*needs a guard*; a proof counts only if it is visible at the call, never as
a claim about callers.

| Group | Sites | Verdict |
|---|---|---|
| `MaddisonSlatkin.cpp`, `ts_data.cpp` `_pad` | 19 | Fixed-size C arrays; `sizeof(data)` cannot be zero |
| `ts_tbr.cpp:930` | 1 | Two scalars (`&kbits`, `&ds.concavity`) |
| `TbrSnapshot::save`/`restore` | 12 | Guarded — PR #134 |
| `ts_bench_tbr_phases` | 12 | Guarded at entry — PR #174 |
| `ts_tree.cpp` `load_tip_states`, `save_node_state` | 8 | Guarded on `tip_bytes > 0` / `total_words == 0` |
| `ts_tree.cpp` `restore_prealloc_undo` | 5 | Unreachable at zero words: the loop is `while (u.count > 0)`, and `save_node_state` returns before incrementing `count` |
| `ts_prune_reinsert.cpp`, `ts_sector.cpp`, `ts_collapsed.cpp` | 3 | Guarded on `tw > 0` / an early return at `total_words == 0` |
| `ts_tbr.cpp` L3b + vroot + NA clip | 8 | Gated by `l3b_active`/`use_directional` (both require `total_words > 0`) and by `has_na` (which implies at least one block) |
| `ts_splits.cpp` | 2 | `wps = (n_tip + 63) / 64 >= 1`; the emit loop's filters are identical to the count loop's, so `idx < n_splits` |

No unguarded site. The sibling patterns `.front()`, `.back()` and
`.data() + i` were swept too (45 + 11 sites); every `.back()` is inside a
`while (!stack.empty())` or on `postorder`, which is never empty.

**Dynamic pass, sub-class 2 — 1704 runs under a local
`-D_GLIBCXX_ASSERTIONS` build**, in three batteries:

| Battery | Runs | Axis |
|---|---:|---|
| 1 | 663 | 13 degenerate shapes (all-constant, all-autapomorphic, all-`?`, all-`-`, single character, 3 and 4 tips, …) x ~30 entry points x {EW, IW} |
| 2 | 382 | Hierarchy/HSJ/XFORM modes, constraints, profile parsimony, resampling API, 24- and 40-tip trees, 12-state characters |
| 3 | 659 | 37 `TS_*` environment knobs, each switching on an alternative kernel, plus the >=150-tip regime that activates L3b incremental edge sets |

Zero aborts. Battery 3 matters most: `l3b_active` requires `n_tip >= 150`
unless `TS_L3B_INCREMENTAL` is set, so eight `ts_tbr.cpp` sites had never
been executed by any test at any point.

**Positive control.** "No aborts" is uninterpretable without proof that the
harness can abort. Reverting the `#151` entry guard and rebuilding produced
`stl_vector.h:1130: Assertion '__n < this->size()' failed` on the very first
degenerate call, `rc=127`; restoring it returned the run to green. The same
control was then run against the new test file, which likewise aborts
without the guard — so that file is not tautological.

## What landed

`tests/testthat/test-ts-degenerate-shapes.R` (Tier 2, 329 expectations,
~5 s). Its expectations are contract checks; its purpose is to put the
degenerate shapes in front of both CI legs on every dispatch, so this class
is watched continuously rather than re-audited by hand after each escape.
It includes an `TS_L3B_INCREMENTAL`-forced block, which buys the L3b sites
on a 12-tip tree instead of 150.

## Residual risk — what this audit does NOT establish

- `_GLIBCXX_ASSERTIONS` hardens libstdc++ containers only. Raw arrays, and
  raw pointers *derived* from a container (`const uint64_t* bits = &v[i];`
  then `bits[w]`), are unchecked; those reads are ASan's job.
- Rcpp vector indexing is not hardened by either flag.
- The dynamic pass is bounded by the shapes imagined. Every historical
  instance was a zero-Fitch-word dataset, so that axis is now covered
  thoroughly and others less so.
- `p + 0` on a null `.data()` (`ts_rcpp.cpp:1750`) is UB before C++20 and
  well-defined from C++20; it is not reported by either detector and was
  left alone.

**Reopening condition.** Any new UBSan `nonnull-attribute` report, any
libstdc++ assertion abort, or any new entry point that indexes per-word
state without a `total_words == 0` guard. A new entry point should be added
to `test-ts-degenerate-shapes.R` in the same commit that introduces it.
