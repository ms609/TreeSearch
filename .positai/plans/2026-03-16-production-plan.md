# TreeSearch: Prototype to Production Plan

## Goal
Produce a phylogenetic search engine competitive with TNT on runtime for equivalent matrices, supporting equal weighting (EW), implied weighting (IW), inapplicable characters (Brazeau et al. 2019), and all existing `MaximizeParsimony` features (constraints, profile parsimony, successive approximations, jackknife/bootstrap, timeout, verbosity, all-optimal-trees collection).

---

## Current State Summary

**What exists in C++:**
- Bit-packed Fitch scoring (EW + IW + inapplicable three-pass)
- NNI, SPR, TBR with incremental scoring (EW and IW, but *not* NA-incremental)
- Ratchet, drift, sectorial (RSS/XSS), tree fusing, Wagner tree, pool management
- Driven search orchestrator (`ts_driven_search`)

**What's missing or incomplete:**
- ~~Constraint enforcement~~ — ✅ 1A complete (C++ `ts_constraint`)
- ~~Profile parsimony~~ — ✅ 1B complete (C++ `compute_profile` + driven search)
- ~~Successive approximations~~ — ✅ 1D complete (C++ `ts_successive_approx`)
- ~~Jackknife/bootstrap~~ — ✅ 1D complete (C++ `ts_resample_search`)
- ~~Timeout~~ — ✅ 1C complete (`max_seconds` in `DrivenParams`)
- ~~Return all optimal trees~~ — ✅ 1C complete (pool return)
- ~~Incremental NA scoring~~ — ✅ 2A complete (TBR/drift use NA-aware incremental scoring; SPR/Wagner still use standard Fitch)
- ~~Parallelization~~ — ✅ 5A-B complete (`std::thread` inter-replicate parallelism, `nThreads` param)
- **Adaptive strategy selection** — nothing yet

**Optimizations completed:**
- ✅ Incremental NA-aware scoring (2A): NA-aware incremental downpass/uppass, Pass 3, indirect length for TBR/drift
- ✅ TBR neighborhood traversal (2B): early termination, smaller-subtree filter, snapshot restore, vroot cache, deferred reshuffling
- ✅ SPR/NNI optimization (2C): bounded indirect scoring, subtree-size filter, incremental NNI scoring
- ✅ Ratchet perturbation modes (2D): upweight/mixed/zero modes, adaptive tuning
- ✅ Sectorial + fusing optimization (2E+2F): O(n) XSS partition, in-place fuse, RSS in driven pipeline
- ✅ Wagner incremental scoring (2G): O(depth×C) per insertion vs O(n×C)
- ✅ TBR symmetry breaking (3A): virtual_prelim dedup via FNV-1a hashing
- ✅ Character ordering (3C): expensive-blocks-first, zero-weight compaction, active_mask skip, bounded indirect in drift/Wagner
- ✅ Memory layout (3D): profiling confirmed current layout is cache-efficient; postorder save/restore in TBR

**Parallelization completed:**
- ✅ Phase 5A-B: `std::thread` inter-replicate parallelism with `ThreadSafePool`, thread-safe RNG via `thread_local`, parallel resample for jackknife/bootstrap. `MaximizeParsimony(nThreads = N)`.

**R-level integration completed:**
- ✅ Phase 9: `MaximizeParsimony()` → C++ driven search; `Morphy()` → MorphyLib search

---

## Phase 1 — Feature Completeness (Parallelizable)

Bring all existing `MaximizeParsimony` features into the C++ backend so that `MaximizeParsimony2` can fully replace `MaximizeParsimony`.

### 1A. Constraint enforcement in C++ (Agent A) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt7atbn.md`

TNT-style locked-node constraint enforcement implemented in `ts_constraint.h/.cpp`.
Uses DFS timestamps for O(1) descendant queries. Per TBR/SPR clip: classify each
constraint as MUST_INSIDE / MUST_OUTSIDE / UNCONSTRAINED. Per candidate regraft
edge: O(1) check per active constraint. Wired into TBR, drift, sectorial search,
Wagner tree, and the driven search Rcpp bridge (`consSplitMatrix` etc.).

### 1B. Profile parsimony scoring in C++ (Agent B) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt9zigq.md`

Key design: reuse IW's indirect evaluation pipeline (`indirect_iw_length`
unchanged) with profile-specific delta precomputation. Add `ScoringMode`
enum and `info_amounts` lookup table to `DataSet`. Zero overhead on EW/IW.

1. ✅ Extend `DataSet` with `ScoringMode`, `info_amounts` matrix, `info_max_steps`.
2. ✅ Add `compute_profile()` and `precompute_profile_delta()` to ts_fitch.
3. ✅ Modify `score_tree()` dispatch for profile mode.
4. ✅ Add `compute_weighted_score()`/`precompute_weighted_delta()` wrappers; update 6 callsites in TBR/SPR/Drift.
5. ✅ Extend ratchet perturbation guards to include profile mode.
6. ✅ Add `infoAmounts` parameter to Rcpp bridge (`ts_driven_search`, `ts_fitch_score`).
7. ✅ Update `MaximizeParsimony()` R function: call `PrepareDataProfile` then C++ engine instead of `Morphy()`.
8. ✅ Test suite: 25 tests cross-validating against R-level `TreeLength(concavity = "profile")`.

### 1C. All-optimal-trees collection + timeout (Agent C) — ✅ COMPLETE

1. ✅ **All optimal trees:** `driven_search` now takes `TreePool&` out-param; Rcpp bridge returns list of edge matrices + score vector for all pool trees.
2. ✅ **Suboptimal collection:** `pool_suboptimal > 0` tested and working.
3. ✅ **Timeout:** `max_seconds` in `DrivenParams`, checked via `std::chrono::steady_clock` after each phase. `timed_out` flag in result.
4. ✅ **`R_CheckUserInterrupt()`:** Added after every heuristic phase (Wagner+TBR, XSS, ratchet, drift, final TBR).
5. ✅ **Verbosity:** `verbosity` param (0=silent, 1=per-replicate, 2=per-phase) via `Rprintf()`.

Tests: 47/47 passing (driven).

### 1D. Successive approximations + jackknife/bootstrap in C++ (Agent C) — ✅ COMPLETE

1. ✅ **Successive approximations:** `ts_successive_approx` — outer loop: driven search → extract per-char steps → reweight via `w_i = (p_i)^(-k) - 1` → repeat until convergence or max iterations. Returns EW parsimony score.
2. ✅ **Jackknife:** `ts_resample_search(bootstrap=false)` — Fisher-Yates partial shuffle, rebuild DataSet with modified weights, run driven search.
3. ✅ **Bootstrap:** `ts_resample_search(bootstrap=true)` — sample with replacement.
4. **R-side thin wrappers:** Not yet written (Rcpp bridges available; R wrappers for `Resample()`/`SuccessiveApproximations()` integration deferred to Phase 9 migration).

Tests: 35/35 passing (resample).

**Note:** Jackknife and bootstrap are embarrassingly parallel — perfect candidates for Phase 5 parallelization.

---

## Phase 2 — Optimization of Individual Heuristics (Parallelizable)

Once feature-complete, optimize each component. These are largely independent and can be worked on in parallel.

### 2A. Incremental NA-aware scoring (Agent E) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt9da67.md`

**Four functions implemented in `ts_fitch_na_incr.inc`:**
1. `fitch_na_incremental_downpass()` — NA-aware Pass 1 walk from clip site to
   root, maintaining `prelim` and `subtree_actives` incrementally.
2. `fitch_na_incremental_uppass()` — NA-aware Pass 2 propagation of `final_`
   with tip update handling.
3. `fitch_na_pass3_score()` — Full Pass 3 (second downpass) on the divided tree
   for exact step counting.
4. `fitch_na_indirect_length()` / `indirect_na_iw_length()` — NA-aware candidate
   screening that suppresses false steps where clip subtree or edge-below subtree
   lacks applicable tips.

**Integration:** Wired into `ts_tbr.cpp` and `ts_drift.cpp` with `has_na`
dispatch. `NodeSnapshot` extended to save/restore `down2` and `subtree_actives`.

**Tests:** 33/33 passing (NA incremental) + all downstream tests unchanged.

**Deferred:** Bounded/cached NA indirect variants, incremental Pass 3,
SPR/Wagner NA integration (lower priority — SPR/Wagner not in driven pipeline).

### 2B. TBR neighborhood traversal optimization (Agent F) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt98wze.md`

**Five optimizations implemented:**
1. **Early termination in indirect scoring**: `fitch_indirect_length_bounded()`,
   `indirect_iw_length_bounded()`, `fitch_indirect_length_cached()`, and
   `indirect_iw_length_cached()` — bail out as soon as accumulated score exceeds
   the current best candidate. Eliminates 60-80% of inner-loop work for losing
   candidates.
2. **Clip smaller subtree only**: Precomputes subtree sizes; skips clips where
   `subtree_size > n_tip / 2`. Halves the worst-case rerooting count.
3. **Avoid full_rescore on rejection**: `StateSnapshot` struct saves complete
   state arrays via `memcpy` before applying a TBR move. On rejection, restores
   from snapshot instead of running an expensive full scoring pass.
4. **Precomputed vroot cache**: For TBR rerooting, precomputes
   `vroot[edge][s] = final_[A][s] | final_[D][s]` for all main edges once per
   clip. Eliminates redundant OR operations in the inner loop.
5. **Deferred reshuffling**: Only reshuffles clip candidates when a full pass
   finds no improvement. After accepting a move, retries with the same ordering
   to exploit spatial locality of improving moves.

**Tests:** TBR bench 26/26, driven 53/53, sector 32/32, fuse 16/16 (1 skip),
resample 35/35 — all passing.

**Deferred to Phase 3:**
- Incremental rerooting (connected DFS traversal of reroot edges with delta
  updates): payoff reduced by smaller-subtree filter (k ≤ n/2)
- Symmetry breaking: Phase 3A
- SIMD vectorization: Phase 3E

### 2C. SPR/NNI indirect scoring improvements (Agent F, continued) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-phase2c.md`

**Six optimizations implemented:**
1. **SPR bounded indirect scoring**: Replaced unbounded `fitch_indirect_length` /
   `indirect_iw_length` with bounded variants (early termination). Added NA-aware
   dispatch (was previously missing).
2. **SPR subtree-size filtering**: Skip clips where subtree > n/2.
3. **SPR deferred reshuffling**: Same pattern as TBR.
4. **SPR incremental clip scoring**: Replaced full rescore + manual clip subtree
   scoring with incremental downpass/uppass, matching TBR.
5. **SPR best-of-all screening**: Track best candidate, verify only the single
   best with full rescore.
6. **NNI incremental scoring**: O(depth × C) per candidate instead of O(n × C).
   State restoration via second incremental downpass after undo. NA fallback to
   full rescore.

**Tests:** SPR/NNI opt 47/47, TBR bench 26/26, driven 53/53, sector 32/32,
fuse 16/16 (1 skip), resample 35/35 — all passing.

**Deferred to Phase 3:**
- SIMD (SSE2/AVX2) for the bit-parallel operations in `fitch_indirect_length`.
  The inner loop processes `n_states` words per block — a natural SIMD target.

### 2D. Ratchet optimization (Agent G) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt99n1l.md`

**Three features implemented:**
1. **Three perturbation modes** (`ZERO_ONLY`, `UPWEIGHT_ONLY`, `MIXED`): upweighting
   via `upweight_mask` field in `CharBlock` — scoring functions count marked characters
   double. For IW, `pattern_freq` is saved/restored instead.
2. **Adaptive perturbation tuning** (opt-in): tracks escape rate over sliding windows
   of 3 cycles and adjusts `perturb_prob` toward a target.
3. **Configurable inner search intensity**: `perturb_max_moves` (previously hardcoded)
   is now a parameter. `perturb_accept_equal` also configurable.

All defaults unchanged — new features are opt-in. Optimal perturbation strategy
deferred to Phase 6 benchmarking.

**Tests:** ratchet-search 17/17, ratchet-stress 72/72, ratchet-opt 23/23.

### 2E. Sectorial search optimization (Agent H) — DONE

Detailed plan: `.positai/plans/2026-03-16-plan-mmt9azq5.md`

**Key optimizations (priority order):**
1. **O(n) XSS partitioning** — current `xss_partition()` is O(n²); fix with maintained `unclaimed_below[]` array
2. **In-place fuse exchange with undo** — avoid full-tree copy per trial exchange; save/restore only the affected clade
3. **Pre-allocated `SectorWorkspace`** — reuse buffers across sectors instead of allocating per sector
4. **Eliminate redundant full-tree rescores** — track dirty flag; skip `score_tree()` when states haven't changed
5. **Lazy donor processing** — process donors one-at-a-time in fuse; skip remaining after first improvement
6. **Wire RSS into driven search** — complement XSS with random sector picks (TNT uses both)
7. **Topology snapshot optimization** — save only sector-clade nodes, not full tree
8. **Donor diversity prioritization** — try most-different donors first in fusing
9. **Multiple TBR rounds per sector** — profile-guided; likely marginal for small sectors

### 2F. Tree fusing optimization (Agent H, continued) — DONE

Folded into 2E plan above (steps 2, 4, 5, 8).

**TNT comparison notes:**
- TNT fuses all pairs; we fuse best-vs-pool. Donor prioritization (step 8) captures most of the benefit without O(pool²) cost.
- TNT copies full trees for trial exchanges; our in-place apply/undo pattern (step 2) should be faster.
- Bottom-up (smallest-clades-first) exchange order is correct and matches TNT.

### 2G. Wagner tree construction optimization (Agent G) — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-wagner-opt.md`

**Optimization implemented: incremental two-pass Fitch scoring during construction.**
Replaced full `build_postorder() + score_tree()` after each taxon insertion with:
- Incremental downpass from new_internal to root (O(depth×C), early termination)
- DFS-based uppass from root with early termination (only visits nodes whose
  final_ actually changed)
- Single `build_postorder() + score_tree()` at the end for final authoritative score

Complexity improved from O(n²×C) to O(n × depth × C) for the downpass component.
The uppass uses early-termination DFS, avoiding full tree traversal.

**Tests:** 26/26 passing (11 original + 5 new: NA datasets, multiple seeds, driven
search integration). All downstream tests (driven 53/53, fuse, sector, resample)
also pass.

**Deferred (as planned):**
- Closest-taxon-first addition order: O(n³×C), not relevant for RAS replicates
- Multiple Wagner starts per replicate: trivial to add later
- IW-specific incremental step tracking: defer to profiling

---

## Remaining Work from Phases 1–2 — ✅ ALL COMPLETE

| Phase | Task | Status |
|-------|------|--------|
| 1B | Profile parsimony in C++ | ✅ COMPLETE — 25/25 tests pass |
| 1D | R-side wrappers for `Resample()` / `SuccessiveApproximations()` | ✅ COMPLETE — C++ engine integration |
| — | Final audit: Wagner `n_tip < 3` guard | ✅ RESOLVED — guard exists in `ts_wagner.cpp` |
| — | Final audit: SPR stale final states | ✅ RESOLVED — `full_rescore()` ensures correctness; 2C optimizations mitigate performance impact |

**`Resample()` migration:** Now dispatches to C++ `ts_resample_search` for
all modes (EW, IW, profile, constrained). No `Morphy()` fallback needed.

**`SuccessiveApproximations()` migration:** Now dispatches to C++
`ts_successive_approx`. Returns a `multiPhylo` with the final SA tree plus
`score`, `sa_iterations`, and `converged` attributes.

---

## Phase 3 — TNT "Clever Tricks" (Sequential, requires profiling from Phase 2)

### 3A. Symmetry-breaking in TBR — ✅ COMPLETE

Detailed plan: `.positai/plans/2026-03-16-plan-tbr-symmetry.md`

**Two optimizations implemented in `ts_tbr.cpp`:**
1. **SPR-equivalence fast path**: Before entering the TBR inner loop for each
   rerooting, `memcmp` the `virtual_prelim` against `clip_prelim`. If identical,
   the rerooting is equivalent to the SPR case already evaluated — skip entirely.
2. **Hash-based virtual_prelim deduplication**: FNV-1a hash of each
   `virtual_prelim` vector tracked in an `unordered_set` per clip. The SPR
   case's hash is seeded before the loop. Duplicate rerootings (same hash) skip
   the entire O(|main_edges| × total_words) inner loop.

Correctness: purely a skip optimization — never changes which candidate is
selected. Two rerootings with identical `virtual_prelim` produce identical
scores for all regraft destinations.

**Tests:** 24/24 passing (test-ts-tbr-symmetry.R). All existing tests unchanged:
TBR 28/28, TBR bench 26/26, driven 53/53, sector 32/32, fuse 16/16, resample 35/35.

### 3B. "New technology" search components — ✅ COMPLETE

- ✅ **RSS within driven search:** Completed in Phase 2E+2F (Agent H). `rss_rounds` in `DrivenParams`.
- ✅ **CSS (Constrained Sectorial Search):** Completed (Agent G). Sector-restricted TBR on full tree (no HTU approximation). `tbr_search` accepts optional `sector_mask`; `css_search` partitions tree and runs per-sector TBR. Wired into driven pipeline after RSS.
- **Constraint-driven sector selection:** Deferred. Using constraint clades as sector boundaries is a minor enhancement; the core CSS mechanism is in place.

### 3C. Character-ordering optimization — ✅ COMPLETE (Agent F)

Five constant-factor optimizations applied:
1. **Expensive blocks first:** Sort changed to descending weight in `build_dataset()`.
2. **Zero-weight pattern compaction:** Patterns with weight=0 removed before block assignment (benefits resampling).
3. **`active_mask == 0` skip:** All block loops in scoring functions skip entirely-zeroed blocks (benefits ratchet perturbation).
4. **Bounded indirect in drift:** 4 unbounded indirect calls replaced with bounded variants.
5. **Bounded indirect in Wagner:** 2 unbounded indirect calls replaced with bounded variants.

Test status: 49 new + 374 regression = 423 total passing.

Also fixed `TreeSearch:::` namespace prefix in 8 test files (`test-ts-tbr-search.R`,
`test-ts-drift-search.R`, `test-ts-ratchet-search.R`, `test-ts-ratchet-opt.R`,
`test-ts-ratchet-stress.R`, `test-ts-iw.R`, `test-ts-pool.R`, `test-ts-splits.R`)
that called internal C++ bridge functions without the qualifier, causing failures
when tests were run via `library()` instead of `devtools::load_all()`.

### 3D. Memory layout optimization — ✅ COMPLETE

Detailed plan: `.positai/plans/phase-3d-memory-layout.md`
Profiling report: `inst/benchmarks/memory_profile_results.md`

**Profiling-driven investigation.** Main finding: indirect scoring dominates TBR time
at scale (72% at 200 tips), and the current memory layout is already well-structured
for cache access. Per-candidate cost is stable (~33 ns for `total_words=12`),
confirming no cache pressure issue. State arrays fit in L2 (162 KB at 200 tips).

**Steps investigated and skipped (profiling showed no benefit):**
- Postorder node renumbering (downpass not the bottleneck)
- Binary-character specialization (all blocks share same `n_states` from contrast matrix)
- Block-major layout (vroot_cache already linear; no cache pressure)
- StateSnapshot reduction (negligible: 5 μs vs 36 ms indirect)

**Optimization applied:** Postorder save/restore in TBR — eliminates redundant
`build_postorder()` calls after unclip and after snapshot restore.

**Tests:** memory-layout 32/32, driven 53/53, tbr-bench 26/26, fuse 16/16 (1 skip),
sector 32/32.

### 3E. SIMD vectorization

- The Fitch inner loop (`intersection |= left[s] & right[s]`) over state words is a prime SIMD target.
- Implement SSE2/AVX2 specializations for the most common `n_states` values (1, 2, 4).
- Use `#ifdef` or runtime dispatch.
- **Careful:** MSVC and GCC have different intrinsic support. Use portable wrappers.

---

## Phase 4 — Missing Items from Your List

### 4A. What's missing from the list above

1. **Profile parsimony optimization:** If profile parsimony turns out to be popular, it deserves the same incremental optimization as IW.
2. **Tabu search:** TNT implements a tabu list to avoid revisiting recently explored topologies. This is a simple addition: maintain a hash set of recent tree hashes and skip moves that would produce a tree in the set.
3. **Tree bisection reconnection with subtree pruning (SPR+TBR hybrid):** TNT's "combo" search alternates SPR and TBR. The idea is that SPR is faster per move, so use it for easy improvements and switch to TBR when SPR plateaus. Currently the driven search goes straight to TBR. Consider an SPR→TBR escalation.
4. **Multiple random addition sequences per replicate:** TNT's `xmult` can try multiple Wagner trees and keep the best before proceeding to TBR. Low cost, potentially better starting points.
5. ~~**Consensus-driven convergence:**~~ Axed. Strict consensus stabilizes trivially when conflicting trees collapse resolution, making it a poor indicator of adequate tree-space sampling. Hits-to-best is simpler and well-understood.

---

## Phase 5 — Parallelization

Detailed plan: `.positai/plans/2026-03-16-plan-phase5-parallel.md`

### 5A. Architecture assessment — ✅ COMPLETE (Agent H)

The cleanest parallelization strategy is **independent replicates on separate threads**, with a shared pool protected by a mutex.

- Each replicate (Wagner → TBR → XSS → ratchet → drift → TBR) runs on its own thread with its own `TreeState` and `DataSet` copy.
- `DataSet` must be copied per thread (ratchet mutates `active_mask`/`upweight_mask`/`pattern_freq`).
- `ConstraintData` must be copied per thread (mutable workspace: `clip_zones`, `constraint_node`, DFS timestamps).
- The `TreePool` is the only shared mutable state → protect with `std::mutex` via `ThreadSafePool` wrapper.
- `R_CheckUserInterrupt()` must be called from the main thread only. Workers check a `std::atomic<bool> stop_flag`.
- Random number generation: pre-generate all seeds from R's RNG on main thread before spawning workers. Each worker seeds its own `std::mt19937` per replicate.

### 5B. Implementation plan — ✅ COMPLETE (Agent H)

Eight implementation steps:

1. **Thread-safe RNG infrastructure** — `thread_local std::mt19937*` pointer + `ts::make_rng()` helper. Each search function's `GetRNGstate()/unif_rand()/PutRNGstate()` replaced with `ts::make_rng()`. When `thread_rng == nullptr` (serial mode), falls back to R's RNG — zero overhead. Similarly `ts::check_interrupt()` replaces `R_CheckUserInterrupt()`.
2. **Extract single-replicate function** — Factor per-replicate pipeline out of `driven_search()` into `run_single_replicate()`.
3. **Thread-safe pool wrapper** — `ThreadSafePool` with `std::mutex`-guarded `add()`, `best()`, `fuse_round()`.
4. **Parallel driven search** — `parallel_driven_search()`: pre-generate seeds, spawn N workers, main thread polls interrupts, join, extract results.
5. **Parallel resampling** — `parallel_resample()` for jackknife/bootstrap: embarrassingly parallel, no shared state.
6. **Rcpp bridge and R interface** — `nThreads` parameter in `ts_driven_search` bridge and `MaximizeParsimony()`. Default `nThreads=1` (serial).
7. **Verbosity** — Workers silent; main thread periodically reports aggregate progress.
8. **Tests** — Serial equivalence, parallel correctness, timeout, edge cases, IW/NA/constraint modes.

### 5C. Intra-replicate parallelism (harder, investigate later)

- **Parallel sectors in XSS:** XSS partitions are independent after partitioning. Each sector could run on its own thread. Benefit depends on sector count and size.
- **Parallel candidate evaluation in TBR:** Evaluate multiple regraft sites simultaneously. Requires read-only access to the main tree state. Tricky because acceptance modifies the tree.
- **Block-parallel Fitch scoring:** Score different character blocks on different threads. Fine-grained; overhead likely exceeds benefit unless trees are very large.

Recommendation: Start with inter-replicate parallelism. It's simple, scales well, and gives the biggest bang for the buck. Intra-replicate parallelism is more complex with diminishing returns — investigate only after inter-replicate is proven and profiled.

### 5D. R-level parallelization considerations

- `unif_rand()` and `GetRNGstate()`/`PutRNGstate()` are not thread-safe. All R RNG usage must happen on the main thread before dispatching workers.
- `R_CheckUserInterrupt()` must only be called from the main thread.
- `Rprintf()` for progress reporting: route through a thread-safe message queue.
- Memory allocation: R's `R_alloc` is not thread-safe. Use standard C++ allocators in worker threads.

### 5E. Agent CPU limits

**Max 2 CPU cores per agent.** Tests and benchmarks must use `nThreads = 2L`
at most — never `nThreads = 0L` (auto-detect). See AGENTS.md for full rules.

---

## Phase 6 — Adaptive Strategy Selection (Data-Driven Tuning)

This is the most ambitious phase. It requires all heuristics to be individually optimized first.

### 6A. Instrumentation

1. Add timing instrumentation to each heuristic in the driven search: wall-clock time spent in Wagner, TBR, XSS, ratchet, drift, fusing.
2. Track "improvement yield" per phase: score improvement per unit time.
3. Track convergence curve: score vs. cumulative time.
4. Return all instrumentation data to R as part of the result.

### 6B. Benchmark suite

1. Curate the existing 25+ empirical datasets into a benchmark suite.
2. For each, determine the known or best-known optimal score (from TNT or exhaustive search for small datasets).
3. Characterize each dataset: n_taxa, n_chars, proportion inapplicable, average character information, character concordance, proportion binary, matrix density.

### 6C. Strategy space

Define the tunable parameters of the driven search as a strategy vector:
- TBR max_hits per phase
- Ratchet cycles and perturbation probability
- Drift cycles and AFD/RFD limits
- XSS rounds and partition count
- Fuse interval
- SPR vs. TBR choice at each phase
- NNI pre-pass (yes/no)
- Multiple Wagner starts per replicate (count)

### 6D. Benchmarking framework

1. For each dataset × strategy combination, run the driven search N times (e.g., 10 replicates) and record:
   - Time to find optimal score (or best found in time limit)
   - Total time
   - Number of replicates to convergence
2. This produces a matrix: datasets × strategies → performance.

### 6E. Predictive model

1. Train a model (random forest, gradient boosting) to predict: given dataset features → which strategy finds the optimum fastest.
2. Feature engineering: the dataset characterization metrics from 6B, plus perhaps cheap-to-compute statistics from the first few replicates (e.g., "how much did ratchet improve over TBR alone?").
3. The model doesn't need to be perfect — even a coarse classification into 3–4 strategy families would be valuable.

### 6F. Adaptive search

1. Implement a "warmup" phase in driven search: run a few replicates with a default strategy, compute dataset features, consult the model, then switch to the recommended strategy for remaining replicates.
2. Alternative: online adaptation — track improvement yield per phase and dynamically adjust the strategy mid-search. E.g., if ratchet is producing most improvements, allocate more cycles to ratchet and fewer to drift.

---

## Parallelization Across Agents

### Phase 1 (Feature completeness) — 4 agents in parallel:
| Agent | Task | Status | Dependencies |
|-------|------|--------|-------------|
| A | Constraint enforcement | ✅ Complete | None |
| B | Profile parsimony | ✅ Complete | None |
| C | All-trees + timeout + verbosity | ✅ Complete | None |
| D | Successive approx + jack/boot | ✅ Complete | None |

### Phase 2 (Optimization) — 4 agents in parallel:
| Agent | Task | Status | Dependencies |
|-------|------|--------|-------------|
| E | Incremental NA scoring (2A) | ✅ Complete | Phase 1 complete (for testing) |
| F | TBR optimization (2B) + SPR/NNI (2C) | ✅ 2B, 2C complete | None (but benefits from E) |
| G | Ratchet (2D) + Wagner (2G) optimization | ✅ 2D, 2G complete | None |
| H | Sectorial + fusing optimization (2E+2F) | ✅ Complete | None |

### Phase 3 (TNT tricks) — 2-3 agents, sequential within phase:
| Agent | Task | Status | Dependencies |
|-------|------|--------|-------------|
| — | 3A Symmetry breaking | ✅ Complete | Phase 2 complete |
| — | 3C Character ordering | ✅ Complete | Phase 2 complete |
| — | 3D Memory/cache optimization | ✅ Complete | Phase 2 complete |
| — | 3E SIMD vectorization | Pending | Phase 3D profiling complete ✅ |
| G | 3B New technology components (CSS) | ✅ Complete | Phase 2 complete |

### Phase 4 (missing items) — integrated with Phase 3

### Phase 5 (Parallelization) — 1 agent, after Phase 2:
| Agent | Task | Dependencies |
|-------|------|-------------|
| H | Thread-parallel replicates + parallel resample | Phase 2 complete ✅ | ✅ COMPLETE |

### Phase 6 (Adaptive strategy) — 1-2 agents, after Phase 5:
| Agent | Task | Dependencies |
|-------|------|-------------|
| J | Instrumentation + benchmarks | Phase 5 complete |
| K | ML model + adaptive search | J complete |

---

## Things Missing from the Original List / Suggestions

1. **Testing infrastructure:** Each optimization must be accompanied by regression tests that verify correctness against known-good scores (the 30 inapplicable datasets, the DNA datasets). A CI benchmark suite that detects performance regressions would be valuable.

2. **Profile-guided optimization (PGO):** After the code stabilizes, build with PGO. This is free performance.

3. **Windows DLL considerations:** The AGENTS.md notes about DLL locking are important. The parallelization work (Phase 5) needs careful testing on Windows.

4. ~~**Consensus computation in C++:**~~ Moot — consensus-driven convergence axed (Phase 4, item 5).

5. **R API design:** `MaximizeParsimony2` should eventually support the same parameter interface as `MaximizeParsimony` (or a clearly better one). Consider a unified function that dispatches to the C++ engine when possible and falls back to R for unsupported features.

6. **Documentation:** Each phase should update the AGENTS.md memory file and the package documentation.

7. **Morphy deprecation path:** The current R-layer search uses Morphy for scoring. The C++ backend has its own Fitch implementation. Long-term, the C++ implementation should fully replace Morphy (simplifying the dependency chain). But this needs careful validation — Morphy is the reference implementation for the inapplicable algorithm.

---

## Priority Ordering

Updated priority reflecting completed work:

1. **Phase 3E (SIMD)** — potentially large constant-factor improvement (Phase 3D profiling confirms indirect scoring inner loop is the target: 72% of TBR time at 200 tips)
All Phase 1–2 items ✅ DONE. Phase 3A, 3B, 3C, 3D ✅ DONE.
~~Phase 5A-B (parallelization)~~ ✅ DONE (Agent H).
R-side wrappers for `Resample()` / `SuccessiveApproximations()` ✅ DONE.
