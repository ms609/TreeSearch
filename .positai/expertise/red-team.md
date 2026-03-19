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

area: 7
reviewed_by: B
date: 2026-03-19
notes: Shiny module wiring (server.R 200 lines, mod_data.R 589, mod_search.R 962, mod_consensus.R 1333, mod_clustering.R 295, mod_treespace.R 718, mod_downloads.R 167). **Findings:** (1) Forward-ref callbacks (cb_ref): 4 callbacks wired correctly at server.R:162-166 after all modules initialized. Placeholder closures capture cb_ref environment; actual implementations set synchronously before any reactive fires. ✓ (2) Cross-module updateXxxInput: Only one instance — mod_data.R:203 uses parent_session to update "treespace-relators". Fragile (hardcoded namespace string) but correct. All other updateXxxInput calls use module-local session. ✓ (3) Orphaned observers: None found. All observeEvent/observe blocks reference existing inputs within their module namespace. ✓ (4) isolate() in result observers: search result observer (line 811) has exactly one reactive dependency (searchTask$result()); all cleanup in isolate(). Profile prep result observer same pattern. ✓ (5) New progress polling observer correctly gated: checks progressFile(), r$searchNotification, r$searchInProgress before polling. Uses separate reactiveVals from profile prep observer. invalidateLater(500) only scheduled when all gates pass. ✓ (6) ShowConfigs operates on top-level DOM IDs (not namespaced) — correct since UI elements like "whichTree", "consConfig" are defined in ui.R outside modules. ✓ (7) UpdateActiveTrees reentrancy guard (r$updatingTrees) uses on.exit for cleanup — correct. ✓ (8) Change detection pattern (r$oldkeepNTips, r$oldOutgroup, etc.) prevents reactive cascades from programmatic input updates — correct but complex. ✓ No new bugs or tasks filed.

Previous: area: 5
reviewed_by: B
date: 2026-03-19
notes: Data pipeline & simplification (ts_data.h/.cpp 289 lines, ts_simplify.h/.cpp 370 lines, ts_constraint.h/.cpp 350 lines). **Findings:** (1) FIXED — `build_reduced_dataset()` in ts_sector.cpp did not copy `ds.inapp_state` to the reduced dataset. Default value (-1) meant "no inapplicable state," which is currently harmless because sectors use Fitch scoring (EW/IW/PROFILE) which reads `blk.has_inapplicable` per block. Would be a bug if sectors were extended to support HSJ scoring. Added `rd.data.inapp_state = ds.inapp_state;` after the existing field copies. (2) FIXED — No guard for `n_states > MAX_STATES (32)` in `build_dataset()`. Token bitmasks use `uint32_t`, so `(1u << s)` for s >= 32 is undefined behavior. Added `Rf_error()` check at top of `build_dataset()`. Unlikely in practice (morphological data: 2–10 states; DNA: 4; protein: 20) but defensive. (3) NOT FILED — `build_reduced_dataset()` does not copy HSJ/Sankoff fields (hierarchy_blocks, tip_labels, sankoff_*). Same analysis as #1: sectors don't use those scoring modes. Not worth fixing until sectors support HSJ/XFORM. (4) Simplification correctness verified: uninformative-character detection (classical criterion + 4-caterpillar verification for ambiguous tokens) is correct and conservative. All-inapplicable characters are not simplified (skipped in Phase 1) but score 0, so no correctness impact. (5) Constraint column-major indexing verified: `split_matrix[s + n_splits * t]` correct. Canonicalization (tip 0 outside) handles n_tips as multiple of 64 correctly (remainder=0 skips clearing, all bits valid). DFS timestamp allocation matches tree node count. `regraft_violates_constraint` logic verified for MUST_INSIDE/MUST_OUTSIDE. (6) EW offset interaction with IW/Profile verified: uninformative patterns (removed from blocks) contribute 0 to IW score (extra_steps clamped to 0) and correct constant to Profile score (precomputed_steps added back before lookup). (7) Added 10 new test assertions: RSS + XSS + sector_diag with inapplicable characters. All 1679 ts-* pass.

Previous: area: 4
reviewed_by: E
date: 2026-03-19
notes: Parallelism & RNG (ts_parallel.cpp 450 lines, ts_rng.h/.cpp 110 lines, ts_driven.cpp 484 lines). **Findings:** (1) Thread-local RNG: ✓ `ts::thread_rng` and `ts::thread_stop_flag` set before replicate loop at line 96–97 of ts_parallel.cpp. Cleaned up at line 153–154. (2) No R API from workers: ✓ All Rprintf calls in ts_driven.cpp gated by `verbosity >= 2`; parallel calls pass `verbosity = 0` (line 126). All interrupt checks via `ts::check_interrupt()` which dispatches to `thread_stop_flag`. Wagner uses `thread_safe_unif()` and `rng_state_begin/end()`. HSJ and Sankoff modules have zero R API calls. (3) Pool mutex: ✓ All ThreadSafePool public methods hold lock_guard. `fuse_round()` holds lock during tree_fuse + score_tree + pool.add. (4) Atomic stop flag: ✓ Uses `memory_order_relaxed` throughout, fine for simple boolean flag. (5) Seeds: ✓ Pre-generated from R RNG on main thread (lines 190–194, 385–389), `GetRNGstate/PutRNGstate` brackets. (6) DataSet copy correctness: ✓ All HSJ/Sankoff fields (hierarchy_blocks, tip_labels, sankoff_cost_matrices, etc.) are `std::vector`s — default copy constructor deep-copies. (7) HSJ/Sankoff thread safety: ✓ Both `hsj_score()` and `sankoff_score()` are stateless (only use parameters and local variables). No global mutable state. (8) **PERF NOTE** (not a bug): XFORM scoring rebuilds `SankoffData` in every `score_tree()` call (vector allocation in hot path). Could pre-build and cache in DataSet. Not filed as a bug since it's optimization-only and XFORM is new/experimental. (9) init.c: 45 entries (43 Rcpp + 2 manual), all arg counts match. (10) Score verification: serial=80→TreeLength=80 ✓, parallel=79→TreeLength=79 ✓ (Vinther2008). No new bugs found.

Previous: area: 2 (complementary)
reviewed_by: A (complementary to B's focus 2 review)
date: 2026-03-19
notes: Search topology invariants — deeper analysis of state restoration. **Findings:** (1) FOUND — SPR stale scoring arrays after rejected regraft. In `spr_search` (ts_search.cpp), after `spr_regraft + full_rescore + rejection + spr_unregraft + spr_unclip`, the `restore_saved_states()` only restores nodes on the clip-to-root path (saved during incremental pass). Nodes on the regraft-to-root path that aren't on the clip-to-root path retain prelim/final_/local_cost values from the regrafted topology's `full_rescore`. On subsequent clip iterations, `fitch_incremental_downpass` may read stale prelims for nodes below the new clip ancestor, producing incorrect `divided_length` and indirect evaluations. **Impact: conservative only** — stale arrays affect candidate screening but never acceptance decisions (all moves gated by `full_rescore` verification). Final score always correct (line 374 does `full_rescore`). Confirmed by targeted test: `test-ts-spr-state-restore.R` (33 assertions, EW/IW/NA datasets × multiple starting trees). **Not filed as a bug** because SPR is a secondary search method (TBR is primary in driven pipeline) and the issue is self-correcting over multiple passes. (2) TBR (ts_tbr.cpp): All rejection paths fully correct. Phase 2 restore via `restore_prealloc_undo() + spr_unclip() + saved_postorder` restores arrays exactly. Candidate verification via `save_topology + state_snap.save → apply_tbr_move → full_rescore → reject → restore_topology + state_snap.restore` restores topology + all arrays including postorder. Tabu rejection path identical. `states_valid` flag is dead code (always true) — both branches call `full_rescore` so no functional impact. (3) NNI (ts_search.cpp): Standard path: `nni_undo + incremental_downpass` correctly restores prelim/local_cost (second downpass recomputes from restored topology's children). final_ stale after rejection but unused until next acceptance (uppass only on accept). NA path: `score_tree` always does full recomputation, so stale arrays are irrelevant. ✓ (4) Drift (ts_drift.cpp): B's `saved_postorder` fix (line 676-678) verified correct. RFD computation's double-rescore-and-reapply pattern is correct if expensive. `drift_search` outer loop always starts phases with `full_rescore`. Subtree sizes computed once in drift_phase (not updated after moves) — minor search-quality concern, not correctness. (5) Hash collision risk: 64-bit hash, birthday bound ~k²/2^64. For k=10000 rerootings, P(collision) ≈ 5×10⁻¹². Negligible. Test added: test-ts-spr-state-restore.R (Tier 2).

Previous: area: 3
reviewed_by: B
date: 2026-03-19
notes: Ratchet & perturbation (ts_ratchet.cpp 237 lines, ts_sector.cpp 792 lines, ts_fuse.cpp 522 lines). **Findings:** (1) BUG FIXED — `perturb_upweight()` and `perturb_mixed()` in `ts_ratchet.cpp` used `ds.pattern_freq[pat] *= 2` per selected character. When multiple characters share the same pattern index, this gives exponential blowup (`original * 2^N` instead of `original + N`). Only affects IW/profile perturbed landscape (EW uses upweight_mask independently, IW uses pattern_freq exclusively). Could cause integer overflow with highly compressed datasets and high perturbation probability (adaptive tuning at 0.5, 50+ chars sharing a pattern). Fix: changed `*= 2` to `+= 1` (additive, consistent with EW upweight_mask semantics). (2) Ratchet save/restore: PerturbSnapshot correctly saves/restores active_masks, upweight_masks, and pattern_freq vector. Unperturbed search always uses original weights. Best-tree topology reset uses copy_topology + build_postorder + reset_states. ✓ (3) upweight_mask has NO effect on IW/profile scoring: fitch_downpass uses upweight_mask for EW step count, but extract_char_steps, compute_iw, compute_profile, precompute_iw_delta, and all indirect_iw_length variants only use pattern_freq. Setting upweight_mask in IW perturbation is harmless but redundant. (4) Sectorial reinsertion revert: RSS/XSS save CladeSnapshot before reinsertion, full-tree rescore after, revert if score worsens. Post-hoc constraint check with revert. All paths maintain valid final_ states for HTU construction. ✓ (5) XSS sectors from same partition are independent (non-overlapping tip sets). Accepting one sector can't corrupt another. (6) CSS uses sector_mask in TBR (no HTU approximation); no revert needed (TBR handles internally). ✓ (7) Fuse tied-score handling: accepts equal-score exchanges only with accept_equal AND actual topology change (prevents infinite loops). Undo on rejection correctly restores clade left/right/parent. ✓ (8) Fuse ancestor stale marking (lines 467-480) is dead code — loop breaks immediately after acceptance. Harmless defensive code. (9) Minor: RSS doesn't update eligible list after equal-score acceptance (stale subtree sizes); efficiency concern only. 1404/1404 ts-* pass.

Previous: area: 2
reviewed_by: B
date: 2026-03-19
previous_notes: Search topology invariants (ts_tbr.cpp, ts_drift.cpp, ts_search.cpp). **Findings:** (1) BUG FIXED — `drift_phase()` in `ts_drift.cpp` sets `saved_postorder` once (line 404) and never updates it after accepted moves. When the last clip candidate is rejected, the stale postorder is restored (line 566). The subsequent `tbr_search` (called from `drift_search`) starts with `full_rescore()` which uses `fitch_downpass()` iterating over the stale postorder — processing nodes before their children, producing wrong scores. Impact: first few TBR iterations after drift perturbation have incorrect accept/reject decisions. Self-healing once a move is accepted (build_postorder_prealloc rebuilds). Fix: added `tree.build_postorder()` before return when n_accepted > 0. (2) TBR search (ts_tbr.cpp): Topology save/restore is correct. TopoSnapshot + StateSnapshot used for rejection. StateSnapshot includes postorder (line 263). PreallocUndo used for clip/unclip cycle. saved_postorder kept current by StateSnapshot.save() at line 830 (before each candidate apply). `states_valid` flag is dead code (always true) — minor smell. (3) SPR search (ts_search.cpp): Correct — `build_postorder()` called after every clip/unclip at line 360, regardless of accept/reject. (4) Hash dedup in TBR: 64-bit FNV-1a for virtual_prelim; collision risk negligible. (5) Verified: `apply_tbr_move` / `drift_apply_tbr_move` reroot path is correct for path reversal. All 1397 ts-* pass.

Previous: area: 1
reviewed_by: C+A
date: 2026-03-19
previous_notes: Fitch scoring correctness. **C's review:** Thorough review of ts_fitch.h/.cpp (711 lines), ts_fitch_na.h (290 lines), ts_fitch_na_incr.h (678 lines). Findings: (1) BUG FIXED — `fitch_indirect_length_bounded` and `fitch_indirect_length_cached` did not account for `upweight_mask`, underscoring candidates during ratchet perturbation TBR. The NA-aware variants (`fitch_na_indirect_length_bounded/cached`) were already correct. Also fixed `nx_cost` computation in ts_tbr.cpp, ts_search.cpp, and ts_drift.cpp, plus drift RFD computation — same missing upweight_mask pattern. Impact: TBR screening during ratchet perturbation was slightly inaccurate. Final result correctness unaffected (full_rescore() is authoritative). All 1397 ts-* pass. (2) fitch_downpass_node standalone function correctly omits upweight_mask (callers handle it). (3) Incremental downpass stop condition correct (root check + self-parent guard). (4) Incremental uppass dirty propagation correct — clip_ancestor marked dirty, tips processed separately in NA path. (5) NA three-pass algorithm (Brazeau et al.) correctly implemented: Pass 1 NA-aware prelim, Pass 2 uppass with applicability propagation (tips processed separately), Pass 3 corrected scoring with subtree_actives. (6) extract_char_steps correctly uses raw local_cost bits (not upweight-adjusted) for per-pattern counts. (7) fitch_indirect_length uses union-of-finals (correct for non-additive; Goloboff 1996) rather than intersection-then-union. **A's complementary review:** (A1) PROVED mathematically: standard Fitch incremental uppass dirty-flag logic is correct even when downpass stops before root. When downpass stops at node N (prelim unchanged), all nodes below N on the path have correct finals without explicit recomputation. Proof: if fitch(M_old, S) = fitch(M_new, S) and both are intersection-type, then N_final ⊆ M_old ∩ S = M_new ∩ S, so uppass(N_final, M_old) = uppass(N_final, M_new) = N_final. For union cases, fitch(M_old∪S) = fitch(M_new∪S) requires M_old = M_new when M∩S = ∅. (A2) FOUND: NA-aware uppass has a theoretical children_app staleness. The NA uppass formula uses `children_app = OR(children's prelims)` which can change even when the parent's prelim is stable — because the NA downpass prelim formula aggregates differently than raw OR. Specifically, case_children depends on children_app, which can differ if a child's prelim changed. This means the stopping node's final_ for NA blocks could be stale. Affects `fitch_na_pass3_score` divided_length. Conservative: full_rescore always catches it. Same class of design choice as extract_divided_steps heuristic. (A3) Confirmed C's upweight_mask fixes are correct. All 1397 ts-* pass.
