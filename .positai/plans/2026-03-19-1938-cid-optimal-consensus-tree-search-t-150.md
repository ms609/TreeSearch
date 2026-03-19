# Plan: CID-Optimal Consensus Tree Search (T-150)

**Date:** 2026-03-19
**Task:** T-150 — CID-optimal consensus tree search
**Branch:** `feature/cid-consensus`

---

## Goal

Implement `CIDConsensus()`: a function that finds a consensus tree minimizing
mean Clustering Information Distance (CID) to a set of input trees, using
TreeSearch's existing R-level search infrastructure (`TreeSearch()`, `Ratchet()`,
`EdgeListSearch()`) with a custom CID scorer.

The proof of concept (briefing-cid-consensus.md) showed CID hill-climbing from
a transfer consensus reduces mean CID from 3.03 (MR) / 2.85 (TC) to **2.16**
on a 20-tip test case.

---

## Architecture

Plug into the existing `TreeScorer` / `EdgeSwapper` / `Bootstrapper` interfaces:

```
CIDConsensus(trees, ...)
  → Ratchet(startTree, cidData,
      InitializeData = identity,
      CleanUpData    = .NoOp,
      TreeScorer     = .CIDScorer,
      Bootstrapper   = .CIDBootstrap,
      swappers       = list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))
```

**Dataset representation:** An environment (reference semantics) containing:
- `$trees` — the input `multiPhylo`
- `$tipLabels` — shared tip labels
- `$nTip` — number of tips
- `$metric` — distance function (default: `ClusteringInfoDistance`)

Using an environment allows `CIDBootstrap` to temporarily swap the tree list
(resample with replacement) without copying large objects.

---

## Phases

### Phase 1: Core CID consensus (binary trees)

Use the existing `Ratchet()` / `TreeSearch()` on **binary** trees with CID
scoring. This alone delivered the POC's best result (2.16 CID via SPR from TC).

#### New file: `R/CIDConsensus.R`

**Exported:**

| Function | Purpose |
|----------|---------|
| `CIDConsensus(trees, ...)` | Main user-facing function |

**Internal:**

| Function | Purpose |
|----------|---------|
| `.CIDScorer(parent, child, dataset, ...)` | TreeScorer: build phylo, compute `mean(metric(candidate, trees))` |
| `.MakeCIDData(trees, metric)` | Create env with trees, tipLabels, metric |
| `.CIDBootstrap(edgeList, cidData, ...)` | Bootstrapper: resample input trees, search, restore |
| `.EdgeListToPhylo(parent, child, tipLabels)` | Helper: parent/child → phylo object |

**`CIDConsensus()` signature:**

```r
CIDConsensus <- function(
  trees,
  metric = ClusteringInfoDistance,
  start = NULL,                    # phylo, or NULL → majority rule consensus
  method = c("ratchet", "spr", "tbr", "nni"),
  ratchIter = 100L,
  ratchHits = 10L,
  searchIter = 500L,
  searchHits = 20L,
  verbosity = 1L,
  ...
)
```

**Design decisions:**

1. **`start = NULL`** defaults to `multi2di(Consensus(trees, p = 0.5))`.
   User can pass any phylo (e.g., a transfer consensus from TreeDist dev).
   Starting tree is always resolved with `multi2di()` for Phase 1.

2. **`method = "ratchet"`** is the default. Delegates to `Ratchet()` with
   `CIDBootstrap` and `list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)`.
   Other methods delegate to `TreeSearch()` with the corresponding swapper.

3. **Lower default `searchIter`/`searchHits`** than parsimony (500/20 vs
   4000/42) because CID scoring is ~100× slower per candidate than Fitch.

4. **`metric` parameter** allows swapping in `MutualClusteringInfo`,
   `SharedPhylogeneticInfo`, or any function with signature `f(tree1, tree2)`.
   Default is `ClusteringInfoDistance`.

5. **Return value:** A `phylo` tree with attributes `"score"` (mean CID)
   and `"hits"`.

6. **`InitializeData = identity`**, **`CleanUpData`** = no-op function.
   The dataset is the env created by `.MakeCIDData()`, passed directly.

**`.CIDBootstrap()` design:**

```r
.CIDBootstrap <- function(edgeList, cidData, EdgeSwapper, maxIter,
                          maxHits, verbosity, stopAtPeak, stopAtPlateau, ...) {
  origTrees <- cidData$trees
  nTree <- length(origTrees)
  cidData$trees <- origTrees[sample.int(nTree, replace = TRUE)]
  on.exit(cidData$trees <- origTrees)
  res <- EdgeListSearch(edgeList[1:2], cidData,
                        TreeScorer = .CIDScorer,
                        EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        verbosity = verbosity,
                        stopAtPeak = stopAtPeak,
                        stopAtPlateau = stopAtPlateau, ...)
  res[1:2]
}
```

Resampling input trees with replacement is the CID analogue of character
bootstrapping — it perturbs the objective function to escape local optima.

#### New file: `tests/testthat/test-CIDConsensus.R`

Tier 2 (skip on CRAN). Tests:

1. `.CIDScorer()` returns correct mean CID for a known tree/tree-set pair.
2. `.CIDBootstrap()` returns a valid edgeList (2 elements, valid topology).
3. `CIDConsensus()` with `method = "spr"` improves or equals starting score.
4. `CIDConsensus()` with `method = "ratchet"` runs without error on a small
   case (10 tips, 20 trees, `ratchIter = 3`).
5. `CIDConsensus()` accepts a custom `metric` (e.g., `MutualClusteringInfo`).
6. `CIDConsensus()` rejects non-multiPhylo input.
7. Starting from user-supplied tree works.
8. Score attribute is set on returned tree.

#### File modifications

| File | Change |
|------|--------|
| `DESCRIPTION` | Add `CIDConsensus.R` to Collate field |
| `NAMESPACE` | `export(CIDConsensus)` + any new TreeDist imports |
| `R/CIDConsensus.R` | New file |
| `tests/testthat/test-CIDConsensus.R` | New file |

### Phase 2: Collapse and Resolve moves (non-binary trees)

Add EdgeSwapper functions that operate on potentially non-binary trees,
enabling the search to find partially-resolved consensus optima.

#### New functions in `R/CIDConsensus.R`:

| Function | Purpose |
|----------|---------|
| `.CollapseSwap(parent, child, nTip)` | Contract a random internal edge → polytomy |
| `.ResolveSwap(parent, child, nTip)` | Resolve a random polytomy → new binary split |
| `.CollapseAllSwap(parent, child, nTip)` | Return list of all single-collapse candidates |
| `.ResolveAllSwap(parent, child, nTip)` | Return list of all single-resolve candidates |

**Collapse algorithm:**
1. Find internal edges (both endpoints > nTip, excluding root edge).
2. Pick one (random, or enumerate all with `edgeToBreak = -1`).
3. Reparent all children of the collapsed child node to its parent.
4. Remove the child node, renumber, return new parent/child.

**Resolve algorithm:**
1. Find polytomy nodes (degree > 2 in parent vector).
2. Pick one polytomy node and two or more of its children.
3. Insert a new internal node as intermediate parent of the selected children.
4. Renumber, return new parent/child.

**Search strategy with non-binary trees:**

A new composite swapper `ConsensusSwap` tries three move types per iteration:
1. With probability p₁: collapse a random edge
2. With probability p₂: resolve a random polytomy
3. With probability p₃: SPR on the fully-resolved version
   (temporarily resolve all polytomies, SPR, then re-collapse original polytomies)

OR, simpler: use `Ratchet()` with a mixed swapper list:
```r
swappers = list(
  .ResolveAndSPRSwap,  # resolve → SPR (coarse)
  .CollapseSwap,       # collapse weak edges
  .ResolveSwap         # refine polytomies
)
```

The key insight: `RearrangeEdges()` and `EdgeListSearch()` don't enforce
bifurcation — only the existing `*Swap` functions do. New swappers that
handle polytomies plug in cleanly.

**Updated `CIDConsensus()` behavior:**
When `start` is a non-binary tree (or when the search reaches a non-binary
optimum via collapse), the mixed swapper list is used automatically.

#### Additional tests

9. `.CollapseSwap()` produces valid non-binary topology.
10. `.ResolveSwap()` produces valid topology with one fewer polytomy degree.
11. Collapse then resolve is reversible (same split count).
12. CIDConsensus with collapse/resolve finds equal-or-better score than
    binary-only search on a case where optimal consensus is non-binary.

### Phase 3: Performance optimizations (stretch)

Not blocking initial delivery; document as future work.

- **Precompute per-tree split entropies** and cache in cidData to avoid
  recomputation across candidates.
- **Parallel candidate evaluation**: `parallel::mclapply` or `future`
  over candidates within `RearrangeEdges`.
- **C++ CID scorer**: Avoid R dispatch overhead by computing CID entirely
  in C++ (would require porting LAP solver).

---

## Implementation order

1. Create `R/CIDConsensus.R` with Phase 1 functions.
2. Write `tests/testthat/test-CIDConsensus.R` (Tier 2).
3. Update `DESCRIPTION` (Collate) and `NAMESPACE`.
4. Run tests, verify on briefing's competing-topology example.
5. Add Phase 2 collapse/resolve functions.
6. Add Phase 2 tests.
7. Verify full search on non-binary case.
8. Update documentation/vignette.

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| CID scoring too slow for large trees (>50 tips, >200 input trees) | Lower default `searchIter`; document performance expectations; Phase 3 optimizations |
| SPR local optima on binary trees miss non-binary optimum | Phase 2 adds collapse/resolve; multi-start with `nSearch` parameter |
| `Ratchet()` bifurcation check (line 89-90) rejects non-binary trees | Phase 2: bypass by calling `EdgeListSearch` directly for non-binary case, or add a `bifurcating = TRUE` parameter |
| TreeDist `TransferConsensus` not yet on CRAN | Default start = majority rule; user can pass any starting tree |
| Edge renumbering after collapse/resolve may break node numbering conventions | Use `ape::collapse.singles()` and `TreeTools::Renumber()` for canonicalization |
