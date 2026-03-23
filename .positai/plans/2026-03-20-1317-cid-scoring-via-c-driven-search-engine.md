# Plan: CID Scoring via C++ Driven Search Engine

**Date:** 2026-03-20
**Task:** Replace R-level `Ratchet()`/`TreeSearch()` backend for `CIDConsensus()`
with the C++ driven search engine (`ts_driven_search`).

---

## Problem

The R-level search is too slow due to:
1. R dispatch overhead per candidate (phylo construction, S3 dispatch for
   `as.Splits()`, R-loop over input trees)
2. No incremental scoring (full CID recomputed for every TBR candidate)
3. No access to driven search machinery (sectorial search, drift, fusing,
   multi-replicate parallelism, pool management)

## Architecture: MRP Screening + CID Verification

**Key insight:** The driven search already separates *screening* (incremental
Fitch scoring for candidate selection) from *verification* (`score_tree()`
for acceptance). HSJ and XFORM modes exploit this — screening uses Fitch on
standard characters; verification dispatches to the full composite scorer.

We apply the same pattern:

1. **MRP (Matrix Representation with Parsimony) characters** encode input tree
   splits as binary Fitch characters. Each input tree's non-trivial splits
   become binary characters (tip in smaller partition → state 1, else → 0).
   With `n_trees` input trees × ~(`n_tips` - 3) splits each, we get a
   synthetic phyDat of ~`n_trees * (n_tips - 3)` binary characters.

2. **Incremental Fitch screening** on MRP characters provides a fast proxy:
   MRP score = total incompatible splits across all input trees ≈ sum of
   Robinson–Foulds distances / 2. This correlates strongly with CID and
   supports the existing TBR indirect-scoring machinery.

3. **CID verification** via `score_tree()` with `ScoringMode::CID` computes
   the actual clustering information distance using TreeDist's C++ MCI
   algorithm. Only called for the best TBR candidate per clip edge (same
   pattern as HSJ/XFORM).

### Benefits

All driven search infrastructure works unchanged:

| Component | CID mode behavior |
|-----------|-------------------|
| Wagner starting trees | MRP-optimal addition trees (diverse, fully resolved) |
| TBR hill-climbing | Incremental MRP screening + CID verification |
| Ratchet | Perturb MRP character weights + CID tree weights (see below) |
| Drift | Accept suboptimal CID moves within AFD/RFD limits |
| Sectorial search | XSS/RSS on MRP characters (CSS if beneficial) |
| Tree fusing | Exchange subtrees between pool members |
| Pool management | Dedup via split hashing, score-based eviction on CID score |
| Multi-replicate parallelism | `std::thread`, independent MRP datasets per thread |

### Ratchet perturbation

MRP characters are organized so that all characters from input tree `i`
occupy consecutive positions across one or more `CharBlock`s. During ratchet
perturbation:

- Standard `perturb_*()` modifies `active_mask`/`upweight_mask` on MRP
  blocks, affecting incremental screening.
- Additionally, `CidData::tree_weights[i]` is set to 0 (zeroed tree) or
  2.0 (upweighted tree) in tandem, so that `cid_score()` uses a perturbed
  CID objective consistent with the MRP perturbation.
- `save_perturb_state()` / `restore_perturb_state()` extended to save/restore
  `CidData::tree_weights`.

This makes CID ratchet equivalent to the R-level `CIDBootstrap` (resample
input trees with replacement).

---

## Phase 1: TreeDist header extraction

TreeDist's MCI computation (`mutual_clustering()` in `tree_distances.cpp`)
and LAP solver (`lap.cpp`) are not currently exposed for `LinkingTo`.

### New files in `../TreeDist/inst/include/TreeDist/`

| File | Contents | Source |
|------|----------|--------|
| `types.h` | `cost`, `lap_dim`, `splitbit`, `int16`, `int32`, constants | `ints.h`, `lap.h` |
| `cost_matrix.h` | `CostMatrix` class | `lap.h` |
| `lap_scratch.h` | `LapScratch` struct | `lap.h` |
| `lap.h` | `lap()` declaration + convenience overload | `lap.h` |
| `mutual_clustering.h` | `mutual_clustering_score()` taking raw split arrays (not Rcpp types) | `tree_distances.cpp` lines 382–548 |
| `clustering_entropy.h` | `clustering_entropy()` from split sizes | inline, ~15 lines |

The extracted headers must be **Rcpp-free** (pure C++ with `<cstdint>`,
`<vector>`, `<algorithm>`) so that TreeSearch can include them without
pulling in Rcpp headers in its own TUs.

**LAP implementation**: TreeSearch compiles its own copy of the LAP solver.
Options:
- (a) `inst/include/TreeDist/lap_impl.h` — header-only, guarded by
  `TREEDIST_LAP_IMPLEMENTATION` define, included in exactly one TreeSearch TU.
- (b) TreeSearch copies `lap.cpp` into `src/ts_lap.cpp`.

Option (a) is cleaner; option (b) is a fallback if alignment sensitivity
causes problems.

**Lookup tables** (`lg2[]`, `lg2_rooted[]`, etc.): These use
`__attribute__((constructor))` for initialization. In the header, declare as
`extern`; provide initialization in a single TU (TreeSearch's `ts_cid.cpp`).

### Changes to existing TreeDist source

- `src/lap.h` → include from `<TreeDist/lap.h>` etc. (thin wrapper)
- `src/tree_distances.cpp` → call shared `mutual_clustering_score()`
- `src/tree_distances.h` → `#include <TreeDist/clustering_entropy.h>` etc.

---

## Phase 2: TreeSearch C++ changes

### 2a. New `ScoringMode::CID`

In `ts_data.h`:
```cpp
enum class ScoringMode { EW, IW, PROFILE, HSJ, XFORM, CID };
```

### 2b. `CidData` struct (`ts_cid.h`)

```cpp
struct CidSplitSet {
  int n_splits;
  int n_bins;
  std::vector<uint64_t> data;    // n_splits * n_bins, contiguous
  std::vector<int> in_split;     // popcount per split

  const uint64_t* split(int i) const {
    return &data[static_cast<size_t>(i) * n_bins];
  }
};

struct CidData {
  int n_trees;
  int n_tips;
  int n_bins;                    // ceil(n_tips / 64)

  // Precomputed input tree splits
  std::vector<CidSplitSet> tree_splits;

  // Precomputed per-tree clustering entropies
  std::vector<double> tree_ce;
  double mean_tree_ce;

  // Per-tree weights (1.0 normally; modified during ratchet)
  std::vector<double> tree_weights;
  double weight_sum;             // sum of tree_weights (for weighted mean)

  // Normalized scoring mode
  bool normalize;

  // Block boundaries: mrp_tree_start[i] = first CharBlock index for tree i
  // Used to synchronize MRP perturbation with tree_weights.
  std::vector<int> mrp_tree_start;

  // LAP scratch (reused across cid_score calls to avoid allocation)
  mutable LapScratch lap_scratch;
};
```

### 2c. `ts_cid.cpp`

Key functions:

| Function | Purpose |
|----------|---------|
| `cid_score(TreeState& tree, const CidData& cd)` | Compute mean CID (or normalized) of candidate against input trees, using precomputed input splits and TreeDist's MCI algorithm |
| `build_mrp_dataset(const CidData& cd, int n_tips)` → `DataSet` | Convert input tree splits to MRP binary characters, packed into `CharBlock`s |
| `build_cid_data(split_matrices, n_tips, normalize)` → `CidData` | Construct `CidData` from R-side split matrices |
| `perturb_cid_weights(CidData& cd, DataSet& ds, mode, prob, rng)` | Synchronize MRP block masks with CidData tree weights during ratchet |
| `save_cid_weights(CidData& cd)` / `restore_cid_weights(CidData& cd)` | Save/restore tree_weights for ratchet |

`cid_score()` implementation:
1. Call `ts::compute_splits(tree)` to get candidate splits
2. Compute candidate clustering entropy
3. For each input tree (weighted by `tree_weights[i]`):
   - Compute MCI via TreeDist's `mutual_clustering_score()` using raw arrays
4. Return `(weighted_mean_CID)` or `1 - weighted_mean(MCI_i / CE_i)` if normalized

### 2d. `score_tree()` dispatch

In `ts_fitch.cpp`:
```cpp
double score_tree(TreeState& tree, const DataSet& ds) {
  if (ds.scoring_mode == ScoringMode::CID) {
    return cid_score(tree, *ds.cid_data);
  }
  // ... existing HSJ, XFORM, EW/IW/PROFILE dispatch ...
}
```

### 2e. DataSet extension

Add to `struct DataSet`:
```cpp
// CID scoring data (populated when scoring_mode == CID).
// Owned by the caller; DataSet holds a non-owning pointer.
CidData* cid_data = nullptr;
```

### 2f. Ratchet integration

Extend `save_perturb_state()` / `restore_perturb_state()` in `ts_ratchet.cpp`:
```cpp
struct PerturbSnapshot {
  // ... existing active_mask / upweight_mask vectors ...
  std::vector<double> cid_tree_weights;  // saved tree_weights
  double cid_weight_sum;
};
```

In `perturb_*()`, after modifying MRP block masks, synchronize
`CidData::tree_weights` by inspecting which tree-blocks were zeroed/upweighted.

### 2g. Rcpp bridge (`ts_rcpp.cpp`)

New exported function:
```cpp
// [[Rcpp::export]]
Rcpp::List ts_cid_consensus(
    Rcpp::List split_matrices,  // list of RawMatrix (one per input tree)
    Rcpp::IntegerVector nTip,
    Rcpp::LogicalVector normalize,
    // DrivenParams fields (same as ts_driven_search):
    int max_replicates, int target_hits, double max_seconds,
    int n_threads, int verbosity,
    // SearchControl fields:
    Rcpp::List control
) {
  // 1. Build CidData from split_matrices
  // 2. Build MRP DataSet from CidData
  // 3. Set ds.scoring_mode = ScoringMode::CID; ds.cid_data = &cid_data;
  // 4. Call ts_driven_search(pool, ds, params)
  // 5. Return pool trees + scores as R list
}
```

### 2h. Registration

Add `ts_cid_consensus` to `TreeSearch-init.c` and run `Rscript check_init.R`.

---

## Phase 3: R-level changes

### `R/CIDConsensus.R` rewrite

Replace the `Ratchet()`/`TreeSearch()` backend with a call to
`ts_cid_consensus()`. The R function becomes a thin wrapper:

```r
CIDConsensus <- function(trees, metric = ClusteringInfoDistance,
                         start = NULL, normalize = FALSE,
                         maxReplicates = 100L, targetHits = 10L,
                         maxSeconds = 0, nThreads = 1L,
                         collapse = TRUE, neverDrop = FALSE,
                         maxDrop = ceiling(NTip(trees[[1]]) / 10),
                         verbosity = 1L, control = SearchControl(),
                         ...) {
  # Input validation (unchanged)
  # Build split matrices from input trees
  tipLabels <- trees[[1]]$tip.label
  splitMats <- lapply(trees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })

  # Call C++ driven search
  result <- ts_cid_consensus(splitMats, NTip(trees[[1]]),
                              normalize, maxReplicates, targetHits,
                              maxSeconds, nThreads, verbosity, control)

  # Post-process: convert edge matrix back to phylo
  bestTree <- result$tree
  attr(bestTree, "score") <- result$score

  # Phase 2: Collapse/resolve (R-level, unchanged)
  if (collapse) {
    cidData <- .MakeCIDData(trees, metric, tipLabels, normalize)
    bestTree <- .CollapseRefine(bestTree, cidData, verbosity)
  }

  # Phase 3: Rogue dropping (R-level, unchanged)
  if (!isTRUE(neverDrop)) {
    # ... existing rogue logic ...
  }

  bestTree
}
```

### API changes

| Old parameter | New parameter | Notes |
|---------------|---------------|-------|
| `method` | Removed | Always driven search |
| `ratchIter`, `ratchHits` | `maxReplicates`, `targetHits` | Match `MaximizeParsimony()` API |
| `searchIter`, `searchHits` | `control = SearchControl()` | Expert tuning via `SearchControl()` |
| — | `maxSeconds` | New: timeout |
| — | `nThreads` | New: parallelism |
| — | `control` | New: expert `SearchControl()` for strategy tuning |
| `metric` | `metric` | Retained for collapse/resolve and rogue phases |

The `metric` parameter is retained for the R-level collapse/resolve and rogue
phases, but the core search always uses CID (via the C++ engine). The `metric`
parameter's scope narrows to post-search refinement only.

### DESCRIPTION

- Add `LinkingTo: TreeDist`
- Version bump if needed

---

## Phase 4: Testing

### C++ unit tests (`tests/testthat/test-ts-cid.R`)

Tier 2. Tests:

1. `ts_cid_consensus()` returns valid tree with score attribute
2. Score matches R-level CID computation (`mean(ClusteringInfoDistance(...))`)
3. Search improves over majority-rule consensus starting score
4. Multi-replicate produces equal-or-better score than single replicate
5. Normalized scoring returns value in [0, 1]
6. `nThreads > 1` produces valid result (correctness, not speed)

### Updated `test-CIDConsensus.R`

Existing tests adapted to new API (parameter name changes). Core assertions
unchanged: search improves score, collapse helps, rogue dropping works,
custom metric accepted.

### Regression test

Verify that on the POC's 20-tip test case, the C++ engine achieves CID ≤ 2.2
(within tolerance of the R-level result of 2.16).

---

## Phase 5: Future work (not in this plan)

- **Collapse/resolve in C++**: Move `CollapseRefine` into the driven search
  as a post-TBR phase.
- **Rogue dropping in C++**: Move the drop/restore loop into C++.
- **Incremental CID**: Track which splits changed during TBR and
  incrementally update MCI (requires LAP warm-starting — research needed).
- **Non-binary search in C++**: Native polytomy handling in the driven search.

---

## Implementation order

1. TreeDist: create `inst/include/TreeDist/` headers
2. TreeSearch: `ts_cid.h/.cpp` with `CidData`, `cid_score()`, `build_mrp_dataset()`
3. TreeSearch: `ScoringMode::CID` + `score_tree()` dispatch
4. TreeSearch: ratchet integration (perturb/save/restore CID weights)
5. TreeSearch: Rcpp bridge `ts_cid_consensus()`
6. TreeSearch: R-level `CIDConsensus()` rewrite
7. Tests + regression check
8. Update `TreeSearch-init.c`, `DESCRIPTION`, `NAMESPACE`

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| MRP screening poorly correlated with CID for some tree distributions | MRP = RF proxy, which is well-correlated with CID in practice. Ratchet/drift escape any screening-induced local optima. Can add CID-based candidate re-ranking if needed. |
| LAP solver alignment sensitivity when compiled in TreeSearch | Use option (b) (copy lap.cpp) if header-only causes regression. Benchmark. |
| CID verification too slow for large trees (>100 tips, >500 input trees) | LAP dimension is small when consensus is close to inputs (many exact matches). LapScratch reuse avoids allocation. Candidate set is already filtered by MRP screening. Profile and optimize if needed. |
| TreeDist header refactor affects TreeDist's own tests | Run TreeDist's full test suite after extraction. Headers are pure refactoring. |
| Per-tree weight synchronization in ratchet is complex | Unit test: verify that after perturb + restore, tree_weights are identical to pre-perturb state. |
