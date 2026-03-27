# CID Consensus: Scaling with Input Tree Count

**Date:** 2026-03-26
**Hardware:** Hamilton 8 HPC (AMD EPYC 7702, 1 core, 8 GB)
**TreeSearch version:** 2.0.0 (commit dc28a696, `feature/cid-consensus`)
**Dataset:** Mammals bootstrap (Lemoine et al. 2018), 1000 trees × 1449 tips,
via Consense package `LoadDataset("mammals_full")`

## Summary

CID verification cost per candidate scales as O(T × n^2.4), dominated by the
O(n^3) LAP solver within `mutual_clustering_info()`.  MRP split deduplication
(97% compression at T=1000) makes the Fitch screening layer essentially
insensitive to T, so the bottleneck is CID verification alone.

Consensus quality (full-set MCI) saturates at T_sub ≈ 50–100 trees and
degrades at higher T_sub under fixed time budgets.  This motivated the
`treeSample` auto-selection in `InfoConsensus()`.

## Experiment 1: Quality vs. T_sub (tree subsampling)

**Design:** For each T_sub value (10–1000), run `InfoConsensus()` on a random
subsample of T_sub trees, then evaluate the result against the full 1000-tree
set.  3 random seeds per condition.

**Tip subsamples:** 50 tips (60s budget) and 100 tips (120s budget), drawn
from the 1439 shared tips via `SubsampleTips()`.

**Data:** `cid_tsub_50tips_20260326_1017.csv`, `cid_tsub_100tips_20260326_1227.csv`

### Key metrics

| Metric | Description |
|--------|-------------|
| `score` | MCI of result against the T_sub subsample (internal objective) |
| `full_mci` | MCI of result against all 1000 trees (true quality) |
| `cid_to_ref` | CID to NCBI reference tree (external validation) |

### Results

**50 tips:** Quality plateaus at T_sub ≈ 20–30.  MCI within 1% of maximum
from T_sub=30 onward.  Wall time barely increases (54s → 62s) due to MRP
dedup (47,000 total splits → 1,217 unique, 97.4% compression).

**100 tips:** Quality peaks at T_sub ≈ 30–70 (median MCI 16.1 bits) and
*degrades* to 15.4 bits at T_sub=1000 (95.2% of peak).  Wall time increases
from 109s to 143s.  The degradation occurs because additional CID cost
displaces search replicates under the fixed time budget.

### Baselines

| Scenario | InfoConsensus MCI | Majority-rule MCI | Strict MCI |
|----------|:-----------------:|:-----------------:|:----------:|
| 50 tips  | 10.7              | 2.27              | 0.40       |
| 100 tips | 16.1              | 3.13              | 0.00       |

InfoConsensus achieves 4–5× the MCI of majority-rule consensus.

### Auto-selection heuristic

Based on these results, `InfoConsensus()` now defaults to
`treeSample = "auto"`, which selects `min(T, max(50, 2 * n_tip))` trees
for Phase 1.  Phases 2 (collapse/resolve) and 3 (rogue dropping) always
use the full tree set.

## Experiment 2: Per-candidate CID verification cost

**Design:** Generate 50 random candidate trees per condition, batch-score
via `ts_cid_score_trees()` against subsamples of the input set.  3 reps
per condition.  This isolates the CID scoring cost from search dynamics.

**Data:** `cid_profile_20260326_1353.csv`

### Per-candidate cost (median, microseconds)

| T_sub |  50 tips | 100 tips |  200 tips |
|------:|---------:|---------:|----------:|
|    10 |    1,420 |    7,300 |    40,740 |
|    20 |    2,820 |   14,060 |    79,960 |
|    50 |    7,100 |   35,280 |   199,540 |
|   100 |   14,320 |   70,760 |   404,780 |
|   200 |   28,480 |  140,920 |   978,840 |
|   500 |   71,180 |  352,740 | 1,998,020 |
| 1,000 |  142,020 |  705,540 | 4,030,760 |

### Scaling analysis

**Linear in T:** Per-tree marginal cost is constant within each tip count
(R² = 1.000):

| Tips | Per-tree cost | Intercept |
|-----:|--------------:|----------:|
|   50 |       142 µs  |     42 µs |
|  100 |       706 µs  |     53 µs |
|  200 |     4,031 µs  | 27,385 µs |

**Power law in n:** Fitting `cost ~ n^α` across tip counts gives
α ≈ 2.41 (R² = 0.999).  The dominant cost is the Jonker–Volgenant LAP
solver within `mutual_clustering_info()`, which is O(n^3) worst-case.
The observed exponent < 3 reflects early exit and hash-based split matching.

### Practical impact on TBR throughput

Each TBR pass evaluates approximately (n_tip − 3) candidates.
CID wall time per TBR pass (seconds):

| T_sub | 50 tips | 100 tips | 200 tips |
|------:|--------:|---------:|---------:|
|    50 |     0.3 |      3.4 |     39.7 |
|   100 |     0.7 |      6.8 |     79.4 |
|   500 |     3.3 |     34.2 |    397   |
| 1,000 |     6.7 |     68.5 |    794   |

At 200 tips with 1000 trees, a single TBR pass takes ~13 minutes in CID
verification alone.  This confirms that subsampling is essential for
large-scale CID consensus.

## Where the time goes

| Component | Cost | Notes |
|-----------|------|-------|
| Input split construction (as.Splits) | ~0.5 ms/tree | One-time; negligible |
| MRP split dedup + Fitch screening | ~0.02 ms/candidate | O(unique_splits × n) |
| Candidate split computation | O(n) per candidate | ~microseconds |
| **CID verification (LAP solver)** | **O(T × n^2.4)/candidate** | **Bottleneck** |

The LAP solver accounts for >99% of per-candidate evaluation cost.
Optimizing input processing (e.g. TreeTools C++ headers for split
construction) would not measurably improve throughput.

## Implications for future optimization

1. **Subsampling (implemented):** T_sub ≈ 50–100 is sufficient for
   quality; auto-selection via `treeSample = "auto"` addresses this.

2. **LAP solver optimization:** The O(n^3) Jonker–Volgenant solver is the
   fundamental bottleneck.  Possible approaches:
   - Approximate MCI via sampling or truncation
   - Sparse cost matrix exploitation (many split pairs have zero overlap)
   - SIMD-accelerated LAP
   - Bounded LAP (early exit when assignment can't beat budget) — already
     partially implemented via `score_budget`

3. **Incremental CID:** Currently every candidate requires full CID
   re-evaluation.  An incremental variant that updates only the affected
   splits after a TBR move would reduce per-candidate cost from O(T × n^3)
   to O(T × k^3) where k = number of changed splits (~4 per TBR move).
   This is the highest-impact optimization path.

4. **Parallel CID:** The T tree evaluations within `cid_score()` are
   independent and could be parallelized across threads.  Currently
   parallelism is only at the inter-replicate level.

## Note on TreeTools C++ headers

TreeSearch already `LinkingTo: TreeTools`, making several C++ headers
available.  Relevance to the CID pipeline:

### `SplitList` (split bitvector construction from RawMatrix)

Structurally identical to our `CidSplitSet` unpack (both are byte → 64-bit
word with popcount).  Could replace our manual unpack in
`ts_cid_consensus()` lines 2727–2749, but this is one-time setup (~0.5 ms
per tree), not the bottleneck.

### `preorder_edges_and_nodes()` (fast preorder renumbering)

Our C++ `cid_tree_from_edge()` + `build_postorder()` serves the same
purpose.  The TreeTools version renumbers nodes in a single pass with
stack-allocated frames.  Could save a few microseconds per
`ts_cid_score_trees()` call, but again not on the critical path.

### `ClusterTable` (Day 1985)

Day's cluster table represents a tree as a list of clusters ⟨L, R⟩ where
each cluster is an interval of leaf internal-labels.  This enables O(n)
shared-cluster counting between two trees (`SETSW` + `SHARED`).

**Potential relevance to CID:**  Phase 1 of `mutual_clustering_info()`
already does O(n) exact-match split identification via hash lookup
(lines 210–236).  ClusterTable could replace this with an O(n) interval
check.  However:

- Our hash-based path is already O(n) amortised.  ClusterTable's advantage
  is constant factors and no hash collisions, but exact matching is a small
  fraction of total cost.
- The LAP phase (which handles non-matching splits) dominates.  ClusterTable
  cannot help here because it only answers containment queries, not the
  pairwise overlap counts (`|a ∩ b|`) needed for the LAP cost matrix.

**Verdict:** ClusterTable would marginally speed up the exact-match phase
of MCI, but the LAP dominates (~95% of cost for dissimilar trees).
Worth flagging for if exact matches become a larger fraction (e.g. with
subsampled bootstrap trees that share substantial structure).  Not a
priority for the current bottleneck.
