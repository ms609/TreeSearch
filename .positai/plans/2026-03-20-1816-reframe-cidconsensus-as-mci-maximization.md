# Plan: Reframe CIDConsensus as MCI Maximization

## Motivation

Currently `CIDConsensus()` minimizes mean CID (Clustering Information Distance).
This creates a normalization headache for rogue taxon dropping: raw CID
naturally decreases when tips are dropped (the CE terms shrink), so dropping
tips looks beneficial even when it isn't. The current workaround is a
`normalize` parameter that switches to `1 - mean(MCI_i / CE_i)`.

The cleaner framing: **maximize mean Mutual Clustering Information (MCI)**.
MCI measures shared clustering information between the consensus and each
input tree. It doesn't have the CE-deflation bias — dropping a tip removes
splits from both sides, so MCI only improves if the tip was genuinely causing
disagreement. This eliminates the need for the `normalize` parameter entirely.

## Design

### Internal representation: negated MCI

The search infrastructure (TreePool, TBR, ratchet, drift, fuse) all assume
**minimization** (lower score = better). Rather than flipping every comparison
operator in both C++ and R, we negate: the internal score is `-mean(MCI)`.
Lower (more negative) = higher MCI = better. On output, we negate back to
present the user with a positive MCI value (higher = better).

### What changes

#### C++: `ts_cid.cpp` — `cid_score()`

Replace both the normalize and non-normalize branches with a single
computation:

```cpp
double cid_score(TreeState& tree, const CidData& cd) {
  SplitSet ss = compute_splits(tree);
  CidSplitSet cand = splitset_to_cid(ss, cd.n_tips);
  // No need for cand_ce — we only need MCI

  double mci_sum = 0.0;
  for (int i = 0; i < cd.n_trees; ++i) {
    if (cd.tree_weights[i] <= 0.0) continue;
    double mci = mutual_clustering_info(cand, cd.tree_splits[i],
                                         cd.n_tips, cd.lap_scratch);
    mci_sum += cd.tree_weights[i] * mci;
  }
  return -mci_sum / cd.weight_sum;  // negated for minimization
}
```

The `normalize` field in `CidData` becomes unused. Keep the field for now
(avoid struct layout churn); mark with a comment. `mean_tree_ce` also becomes
unused for scoring but can be kept for diagnostics.

#### C++: `ts_cid.h`

- Update `cid_score()` doc comment.
- Add comment noting `normalize` and `mean_tree_ce` are vestigial.

#### C++: `ts_rcpp.cpp` — `ts_cid_consensus()`

- Keep the `normalize` parameter in the Rcpp signature for backward
  compatibility but ignore it (or emit a deprecation message).
- The `mean_tree_ce` computation can stay (harmless).
- **The returned `best_score` and per-tree scores are now negated MCI.**
  R-side will negate on receipt.

#### R: `CIDConsensus.R`

1. **`CIDConsensus()` signature**:
   - Remove `normalize` parameter (or keep with deprecation warning if
     passed).
   - Remove the normalization warning (lines 135–141).
   - Update roxygen docs: "maximizes mean MCI", score is "mean MCI
     (higher is better)".

2. **`.CIDScoreFast()`**:
   Simplify to compute `-mean(MCI)`:
   ```r
   .CIDScoreFast <- function(parent, child, dataset) {
     nTip <- dataset$nTip
     candidate <- .EdgeListToPhylo(parent, child, dataset$tipLabels)
     candSp <- as.Splits(candidate, dataset$tipLabels)
     nTree <- length(dataset$inputSplitsRaw)
     mciSum <- 0
     for (i in seq_len(nTree)) {
       mciSum <- mciSum + MutualClusteringInfoSplits(
         candSp, dataset$inputSplitsRaw[[i]], nTip)
     }
     -mciSum / nTree
   }
   ```
   Remove the normalize branch and CE computation entirely.

3. **`.MakeCIDData()`**:
   - Remove `normalize` parameter.
   - Keep CE precomputation (used if someone passes a custom metric; harmless
     otherwise).

4. **`.CIDDrivenSearch()`**:
   - Pass `normalize = FALSE` to C++ (value is ignored but signature expects it).
   - Negate `result[["best_score"]]` before attaching as attribute:
     `attr(tree, "score") <- -result[["best_score"]]`
   - Same for hits.

5. **`.TopologySearch()`**:
   - Remove `normalize` from the `.NullOr()` call.
   - Negate returned score.

6. **`.CollapseRefine()`**:
   - No comparison changes needed (all `<` comparisons still work with
     negated MCI).
   - Negate the final score before attaching to the output tree.

7. **`.RogueRefine()`**:
   - Remove `originalNormalize` variable and its propagation.
   - No comparison changes needed.
   - Negate scores on output.

8. **`.PrescreenMarginalNID()`**, **`.BestInsertion()`**:
   - Remove `normalize` parameter passing to `.MakeCIDData()`.

9. **`.CIDBootstrap()`**:
   - No changes needed (operates on internal scores).

10. **Output score convention**:
    All `attr(result, "score")` attachments must negate the internal
    score to produce positive MCI for the user. Audit every code path
    that sets this attribute.

#### R: `CIDConsensus.R` — messaging

Update `message()` calls in verbosity output to say "MCI" instead of
"score" where appropriate, and note that higher is better.

#### Tests: `test-CIDConsensus.R`

- Remove or update tests that check `score >= 0` (MCI scores are positive
  but the bound is `[0, max_CE]`, not `[0, 1]`).
- Remove or update normalized-scoring tests (lines 386–407): the
  `normalize` parameter goes away.
- Remove the "Warning fires for rogue dropping without normalization"
  test (lines 516–524).
- Update the rogue-dropping tests to not pass `normalize = TRUE`.
- Update collapse-comparison test (line 381): direction is now `>=`
  (higher MCI = better), or keep `<=` if comparing internal negated
  scores — depends on whether the test uses `attr(result, "score")`
  (which will be positive MCI).
- Add a test verifying that `attr(result, "score")` matches
  `mean(MutualClusteringInfo(result, trees))`.

#### Tests: `test-ts-cid.R`

- Update `CIDConsensus` wrapper tests analogously.

### What doesn't change

- `mutual_clustering_info()` — unchanged.
- `clustering_entropy()` — unchanged (still used by MRP builder, diagnostics).
- `build_mrp_dataset()` — unchanged.
- All search infrastructure (TreePool, TBR, ratchet, drift, fuse) — unchanged
  (they see negated MCI, which is still a "lower is better" score).
- The `metric` parameter on `CIDConsensus()` — kept for non-CID metrics.
  The generic path (non-CID metric) still minimizes `mean(metric(...))`.
  This is a separate concern; leave as-is for now.
- `sync_cid_weights_from_mrp()` — unchanged (weight sync is metric-agnostic).

### Risk assessment

- **Low risk**: The core algorithm (LAP, MCI, MRP screening) is untouched.
  The change is purely in how the objective is aggregated and presented.
- **Medium risk**: Score comparison direction. All internal comparisons
  use `<` (lower = better) and we preserve this via negation. The only
  places where direction matters are user-facing: `attr(result, "score")`
  and verbosity messages. These need careful auditing.
- **Test coverage**: Existing tests cover scoring, collapse, rogue dropping,
  and end-to-end. Many will need updating but the coverage itself is good.

## Task sequence

1. C++ changes: `ts_cid.cpp`, `ts_cid.h` (simplify `cid_score()`)
2. C++ bridge: `ts_rcpp.cpp` (ignore `normalize`, keep param for compat)
3. R scoring: `.CIDScoreFast()`, `.MakeCIDData()` (remove normalize logic)
4. R pipeline: `.CIDDrivenSearch()`, `.TopologySearch()`, `.CollapseRefine()`,
   `.RogueRefine()`, helpers (remove normalize plumbing, negate on output)
5. R API: `CIDConsensus()` signature and docs
6. Tests: update `test-CIDConsensus.R` and `test-ts-cid.R`
7. Build and run full test suite
