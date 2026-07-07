# Architecture Reference

Load this when: editing `src/ts_*.cpp`/`.h`, adding Rcpp exports, reading
the R-level API, or reviewing design decisions.

---

## R-level API

| Function | Engine | Purpose |
|----------|--------|---------|
| `MaximizeParsimony()` | C++ driven search | Primary search (EW, IW, profile, constraints) |
| `Resample()` | C++ | Jackknife/bootstrap resampling |
| `SuccessiveApproximations()` | C++ | Successive approximations weighting |
| `TreeLength()` | C++ `ts_fitch_score` | Score one or more trees |
| `FastCharacterLength()` | C++ `ts_char_steps` | Per-character step counts |
| `AdditionTree()` | C++ `ts_wagner_tree` | Wagner tree construction |
| `RandomTreeScore()` | C++ (phyDat) or MorphyLib (morphyPtr) | Score a random tree |
| `TaxonInfluence()` | C++ via `MaximizeParsimony()` | Per-taxon search |
| `SearchControl()` | — | Expert parameter constructor for `MaximizeParsimony()` |
| `ParsSim()` | Pure R | Simulate datasets under parsimony (EW/IW/profile) |

---

## C++ module map

| Module | Header/Source | Purpose |
|--------|--------------|---------|
| Fitch scoring | `ts_fitch.h/.cpp` | Downpass, uppass, incremental, indirect |
| NA scoring | `ts_fitch_na.h` | Three-pass inapplicable algorithm (Brazeau et al. 2019) |
| NA incremental | `ts_fitch_na_incr.h` | Incremental NA-aware scoring for TBR/drift |
| SIMD | `ts_simd.h` | SSE2/NEON portability layer for bit-parallel ops |
| Data | `ts_data.h/.cpp` | `DataSet`, `CharBlock`, `build_dataset`, simplification |
| Tree | `ts_tree.h/.cpp` | `TreeState`, topology manipulation, `PreallocUndo` |
| Constraint | `ts_constraint.h/.cpp` | Topological constraint enforcement |
| TBR | `ts_tbr.h/.cpp` | TBR search (with sector_mask for CSS) |
| SPR/NNI | `ts_search.h/.cpp` | SPR and NNI search (standalone, not in driven pipeline) |
| Ratchet | `ts_ratchet.h/.cpp` | Perturbation (zero/upweight/mixed, adaptive) |
| Drift | `ts_drift.h/.cpp` | Accept suboptimal moves within AFD/RFD limits |
| Wagner | `ts_wagner.h/.cpp` | Greedy addition tree (incremental scoring, NA-aware) |
| Sectorial | `ts_sector.h/.cpp` | RSS (conflict-guided), XSS, CSS; from-above HTU |
| Fuse | `ts_fuse.h/.cpp` | Tree fusing (in-place exchange) |
| Pool | `ts_pool.h/.cpp` | Dedup, eviction, consensus hash, split frequency table |
| Splits | `ts_splits.h/.cpp` | Bipartition computation, comparison, `hash_single_split()` |
| Driven | `ts_driven.h/.cpp` | Multi-replicate orchestrator |
| Resample | `ts_resample.h/.cpp` | Jackknife, bootstrap, successive approximations |
| Parallel | `ts_parallel.h/.cpp` | `std::thread` inter-replicate parallelism |
| RNG | `ts_rng.h/.cpp` | Thread-safe RNG (`thread_local` dispatch) |
| Simplify | `ts_simplify.h/.cpp` | Character compression and uninformativeness checks |
| Collapsed | `ts_collapsed.h/.cpp` | Zero-length edge detection for clip skipping |
| NNI perturb | `ts_nni_perturb.h/.cpp` | Stochastic NNI-perturbation (IQ-TREE-style topology escape) |
| HSJ scoring | `ts_hsj.h/.cpp` | Hopkins & St. John hierarchy scoring |
| Sankoff | `ts_sankoff.h/.cpp` | Sankoff step-matrix scoring (x-transform) |
| Rcpp bridge | `ts_rcpp.cpp` | All Rcpp-exported functions |

---

## Scoring modes

`ScoringMode` enum in `ts_data.h`: `EW`, `IW`, `PROFILE`, `XFORM`.
- **EW**: standard Fitch parsimony
- **IW**: implied weights via `e/(k+e)` where `e = steps - min_steps`
- **PROFILE**: lookup in `info_amounts` table (structurally identical to IW pipeline)
- **XFORM**: Fitch(non-hierarchy) + Sankoff(recoded composite characters)

Profile mode sets `ds.concavity = 1.0` (finite sentinel) so existing
`isfinite()` checks activate the weighted pipeline without code duplication.

---

## Parallelism design

- `std::thread` (not OpenMP) to avoid R memory allocator conflicts
- Per-thread: `DataSet` copy, `ConstraintData` copy, `std::mt19937` RNG
- Shared: `ThreadSafePool` (mutex-guarded), atomic stop flag
- Main thread: pre-generates seeds from R's RNG, polls
  `R_CheckUserInterrupt()` and timeout every 200ms
- Worker threads make no R API calls — `ts_rng.h` provides `thread_local`
  dispatch (null → R API for serial; set → thread-local for parallel)

---

## Scoring notes

- `.h` file changes (`ts_fitch_na.h`, `ts_fitch_na_incr.h`) may require
  `touch src/ts_fitch.cpp` before rebuild if the build system doesn't track
  header dependencies.
- Incremental scoring is a **screening heuristic** for candidate selection;
  `full_rescore()` / `score_tree()` is always authoritative.
- See `.positai/expertise/fitch-scoring.md` for detailed invariants:
  uppass correctness proof, NA staleness analysis, `upweight_mask` audit.

---

## Constraint enforcement

- `build_constraint()` reads R split matrix with **column-major** indexing:
  `split_matrix[s + n_splits * t]`.
- Wagner uses LCA-based constraint mapping (`wagner_map_constraint_nodes`)
  since splits aren't fully present during incremental construction.
- Wagner has a posthoc retry loop (up to 100 random addition orders) as a
  safety net for edge cases.

---

## Exported Rcpp functions

All registered in `ts_rcpp.cpp` and `TreeSearch-init.c`. Run
`Rscript check_init.R` to verify consistency.

| Function | Module | Purpose |
|----------|--------|---------|
| `ts_fitch_score` | ts_fitch | Score a tree |
| `ts_char_steps` | ts_rcpp | Per-pattern step counts (with simplification offsets) |
| `ts_na_debug_char` | ts_fitch_na | Per-node debug for a single pattern |
| `ts_na_char_steps` | ts_fitch_na | Per-pattern step counts (raw, no offsets) |
| `ts_debug_clip` | ts_fitch | Debug SPR clip/regraft |
| `ts_test_indirect` | ts_fitch | Debug indirect length |
| `ts_nni_search` | ts_search | NNI hill-climbing |
| `ts_spr_search` | ts_search | SPR hill-climbing |
| `ts_tbr_search` | ts_tbr | TBR with plateau exploration |
| `ts_ratchet_search` | ts_ratchet | Ratchet perturbation |
| `ts_drift_search` | ts_drift | Drift search |
| `ts_wagner_tree` | ts_wagner | Wagner tree (specified addition order) |
| `ts_random_wagner_tree` | ts_wagner | Wagner tree (random order) |
| `ts_compute_splits` | ts_splits | Bipartition splits from edge matrix |
| `ts_trees_equal` | ts_splits | Compare two trees |
| `ts_pool_test` | ts_pool | Pool deduplication test |
| `ts_tree_fuse` | ts_fuse | Fuse two trees |
| `ts_sector_diag` | ts_sector | Sectorial search diagnostics |
| `ts_rss_search` | ts_sector | Random Sectorial Search |
| `ts_xss_search` | ts_sector | Exclusive Sectorial Search |
| `ts_driven_search` | ts_driven | Full driven search |
| `ts_resample_search` | ts_resample | One jackknife/bootstrap replicate |
| `ts_successive_approx` | ts_resample | Successive approximations |
| `ts_parallel_resample` | ts_parallel | Batch resample with parallelism |
| `ts_bench_tbr_phases` | ts_rcpp | TBR phase timing diagnostic |
| `ts_hsj_score` | ts_hsj | HSJ hierarchy scoring |

---

## Key design decisions

1. **PreallocUndo** (`ts_tree.h`): Pre-allocated flat buffers for TBR/drift
   undo stack. Uses `grow()` to dynamically expand when capacity exceeded
   (NA uppass saves both internal nodes and tips). Initial capacity `3 * n_node`.

2. **TBR symmetry breaking** (`ts_tbr.cpp`): FNV-1a hash deduplication of
   `virtual_prelim` vectors to skip redundant rerooting evaluations.

3. **Bounded indirect scoring**: All search modules use `_bounded` variants
   that bail out when accumulated score exceeds best candidate.

4. **Profile parsimony**: Reuses IW indirect pipeline unchanged; only delta
   precomputation differs. `ds.concavity = 1.0` sentinel activates weighted
   path. Max 2 informative states per character; inapplicable → ambiguous.

5. **MPT enumeration**: Post-search TBR plateau walk from all pool seeds.
   `tbr_search()` accepts optional `TreePool* collect_pool` parameter.

6. **All-ambiguous phyDat guard**: `TreeLength()` and `MaximizeParsimony()`
   check for `levels = NULL` / 0-column contrast matrix before calling C++.

7. **From-above HTU for sectorial search** (`ts_sector.cpp`):
   `compute_from_above_for_sector()` computes `from_above[sector_root]` —
   the Fitch state-set the rest of the tree sends *down* to the sector
   boundary, excluding the sector's own contribution. Used instead of
   `final_[parent]` in `build_reduced_dataset()`. O(depth × total_words).

8. **Split frequency table** (`ts_pool.h/.cpp`): `SplitFrequencyTable` maps
   per-split FNV-1a hash → occurrence count across best-score pool trees.
   Used by conflict-guided RSS to weight sector selection. The same FNV-1a
   hash (`hash_single_split()` in `ts_splits.h`) is used by consensus
   hashing and split frequency counting — must stay consistent.

9. **Consensus-stability hash** (`ts_pool.cpp`): XOR of FNV-1a hashes of
   splits present in ALL best-score trees. Updated after each replicate.
   Hash collision false-matches are conservative (over-count stability).

10. **Diversity-aware pool eviction** (`ts_pool.cpp`): When the pool is full
    and a new tree ties the worst score, the entry most similar to the new
    tree (most shared splits, counted via per-split FNV-1a hash set
    membership) is evicted. Falls back to arbitrary worst entry when the
    new tree is strictly better.

11. **Cross-replicate consensus constraint tightening** (`ts_driven.cpp`):
    When `consensus_constrain = true` and no user constraint is supplied,
    after ≥5 replicates, unanimous pool splits are extracted and enforced
    as topological constraints via `build_constraint_from_bitsets()`. The
    TBR/SPR search then avoids breaking established consensus clades.
    Constraints are cleared and rebuilt whenever the best score changes.
    Sector/fuse operations do not enforce auto-constraints.
