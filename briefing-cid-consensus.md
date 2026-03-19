# Briefing: CID-Optimal Consensus Trees via TreeSearch

**Date:** 2026-03-19  
**From:** TreeDist dev cycle  
**Status:** Proof of concept validated; ready for design work

---

## Summary

TreeSearch's SPR/TBR/NNI rearrangement machinery and Ratchet metaheuristic
can be used to find consensus trees that minimise Clustering Information
Distance (CID) to a set of input trees.  A proof-of-concept plugging a CID
scorer into `TreeSearch()` via the existing custom-scorer interface works
out of the box.  This briefing describes the opportunity and implementation
considerations.

---

## Background

The **median tree** problem: given N input trees (e.g., bootstrap replicates
or MCMC posterior samples), find the tree C that minimises the mean distance
to all inputs.  The majority-rule (MR) consensus solves this for
Robinson-Foulds distance but produces poorly-resolved trees when signal is
weak.

Takazawa et al. (2025, bioRxiv 10.64898/2026.03.16.712085) proposed
**transfer-distance consensus** — a greedy algorithm that constructs more
resolved trees by optimising a finer-grained split dissimilarity.  Their
algorithm is implemented in TreeDist as `TransferConsensus()` (pure R,
working prototype).

**The idea here:** go further and optimise TreeDist's CID (or PID, MSID,
etc.) directly, using TreeSearch's tree-rearrangement infrastructure.  CID
is an information-theoretic metric based on bipartite matching (solved via
the Linear Assignment Problem), so it captures richer structure than either
RF or transfer distance.

---

## Proof of Concept (validated)

A custom CID scorer was plugged into `TreeSearch()` via the
`TreeScorer` / `InitializeData` / `CleanUpData` interface:

```r
CIDScorer <- function(parent, child, dataset, ...) {
  tr <- structure(list(
    edge = cbind(parent, child),
    tip.label = names(dataset),
    Nnode = length(unique(parent))
  ), class = "phylo")
  mean(ClusteringInfoDistance(tr, attr(dataset, "trees_input")))
}
```

This works with `SPRSwap`, `TBRSwap`, `NNISwap`, and `Ratchet()`.

On a test case (20-tip trees, competing topologies):

| Method | Mean CID |
|---|---|
| Majority-rule | 3.03 |
| Transfer consensus | 2.85 |
| SPR hill-climb from TC | **2.16** |

The Ratchet bootstrapper adapts naturally: instead of resampling
characters, resample input trees with replacement.

---

## What TreeSearch Would Need

### 1. Additional rearrangement moves: collapse & resolve

SPR/TBR/NNI operate on **bifurcating** trees.  But the CID-optimal
consensus may be partially unresolved — collapsing a weakly-supported split
can reduce CID by removing a bad match from the LAP.

Two new move types are needed:

- **Collapse**: contract an internal edge to create a polytomy.  Equivalent
  to removing one split from the tree.  Neighbourhood size: O(n).
- **Resolve**: pick a polytomy node and resolve it into a binary split.
  The number of binary resolutions of a k-way polytomy is (2k−5)!! — but
  restricting to single-step resolutions (splitting one multifurcation into
  two groups) gives O(k²) candidates per polytomy.

These moves allow the search to explore non-binary trees, which is essential
for consensus tree optimisation.  They complement SPR/NNI/TBR, which handle
the binary rearrangement landscape.

A natural search strategy:

1. Start from the transfer consensus (or MR).
2. SPR/TBR within the binary sub-landscape.
3. Periodically try collapse moves on each internal edge.
4. After collapsing, try resolve moves to explore alternative resolutions.
5. Ratchet for escaping local optima.

### 2. Performance: batch CID scoring

The CID scorer is currently called once per candidate tree, computing
`mean(ClusteringInfoDistance(candidate, all_input_trees))`.  This
involves N LAP solutions (one per input tree).  For 50-tip trees with
100 inputs:

- Per-candidate cost: ~100 × 12 μs = 1.2 ms
- SPR neighbourhood: ~4,700 candidates → ~5.6 s per iteration
- Full search (100 iterations): ~10 min

This is usable but not fast.  Potential optimisations:

- **Incremental LAP update**: when an SPR changes only a few splits, most
  of the cost matrix is unchanged.  A warm-started LAP (seeded with the
  previous solution) could skip most of the Dijkstra augmentation phase.
  This is a research question — LAPJV doesn't natively support warm starts.
- **Neighbourhood pruning**: evaluate a random subset of SPR moves per
  iteration rather than the full neighbourhood (stochastic hill-climbing).
  TreeSearch already supports this via `EdgeSwapper` returning single moves.
- **Parallel candidate evaluation**: each candidate tree is independent;
  OpenMP over the candidate set (orthogonal to the existing per-pair
  parallelism).
- **C++ scorer**: currently the scorer goes R → C++ → R per candidate.
  A C++ `TreeScorer` that avoids R dispatch overhead would be faster.
  TreeSearch's existing C++ scoring interface (via Morphy) provides a
  template.

### 3. Interface design

A natural R interface:

```r
CIDConsensus(trees,
             start = c("transfer", "majority", "random"),
             method = c("spr", "ratchet", "tbr"),
             metric = ClusteringInfoDistance,
             maxIter = 1000,
             ...)
```

- `start`: how to initialise the search tree.
- `method`: search strategy.
- `metric`: which TreeDist metric to optimise (CID recommended, but PID,
  MSID, IRF could also work).
- Returns: a `phylo` tree (possibly non-binary) with the optimised score
  as an attribute.

### 4. Metrics where this applies

Any TreeDist metric can serve as the objective function.  The decomposition
properties affect which optimisations are available:

| Metric | LAP? | Exact-match shortcut? | Notes |
|---|---|---|---|
| CID | Yes | Yes | Recommended primary metric |
| PID (= SPI) | Yes | **No** (see below) | Exact-match detection incorrect for SPI |
| MSID | Yes | **No** (same issue) | |
| MSD | Yes | Yes | Simpler (no IC terms) |
| IRF | No | Yes | No LAP — just per-split info sum |
| RF | No | N/A | MR already optimal for RF |

**Important correctness note for SPI/MSI**: exact-match detection (matching
identical splits before solving the LAP) is *incorrect* for SPI and MSI.
`spi_overlap(A, B)` where B ⊇ A can exceed `spi_overlap(A, A)`, so the
full LAP may find a better assignment by NOT matching identical splits.
This bug was found and fixed in the TreeDist C++ batch path during this
dev cycle.

---

## Hill-Climbing Experiment Results (from TreeDist)

Tested iterative single-split add/remove (not SPR) across three scenarios:

| Scenario | MR | TC | HC from TC | HC from MR |
|---|---|---|---|---|
| Competing topologies (20-tip) | 3.03 | 2.85 | **2.16** | 2.84 |
| Mixed random (30-tip) | **11.54** | 13.91 | 12.46 | 11.54 |
| Heterogeneous noise (40-tip) | 5.06 | 5.06 | 5.06 | 5.06 |

Key findings:
- CID hill-climbing can substantially improve over both TC and MR.
- **Starting point matters**: different starts reach different local optima.
  Best-of-multiple-starts is advisable.
- When TC ≈ MR, the consensus is already near a CID local optimum.

---

## Relationship to Existing TreeSearch Architecture

The existing `EdgeListSearch` → `TreeScorer` → `EdgeSwapper` loop is the
right abstraction.  CID scoring slots in as a custom `TreeScorer`.  The
main gaps are:

1. **Non-binary tree support**: current `EdgeSwapper` functions assume
   binary trees.  Collapse/resolve moves need new swappers.
2. **Ratchet bootstrapper**: `CIDBootstrap` (resample input trees) replaces
   `MorphyBootstrap` (resample characters).  Straightforward.
3. **Score interpretation**: TreeSearch minimises parsimony (integer score).
   CID is a real-valued distance.  The search loop uses `<` comparison,
   so this works, but `stopAtScore` and `suboptimal` parameters need
   sensible defaults for CID-scale values.

---

## References

- Takazawa Y, Takeda A, Hayamizu M, Gascuel O (2025). "Outperforming the
  majority-rule consensus tree using fine-grained dissimilarity measures."
  bioRxiv. doi:10.64898/2026.03.16.712085
- Smith MR (2020). "Information-theoretic Generalized Robinson-Foulds
  metrics for comparing phylogenetic trees." Bioinformatics, 36(20),
  5007–5013.
