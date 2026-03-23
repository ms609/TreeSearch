# Implementation Plan: Reframe CIDConsensus as MCI Maximization

Based on the design in `2026-03-20-1816-reframe-cidconsensus-as-mci-maximization.md`.

## Summary

Replace CID minimization with MCI maximization. Internal score becomes
`-mean(MCI)` (negated for the minimization infrastructure). User-facing
scores are positive MCI (higher = better). The `normalize` parameter is
removed entirely.

## Step-by-step

### 1. C++: `src/ts_cid.cpp` — simplify `cid_score()`

Replace the two-branch `if (cd.normalize) { ... } else { ... }` body
(lines 518–560) with a single branch that computes `-weighted_mean(MCI)`:

```cpp
double cid_score(TreeState& tree, const CidData& cd) {
  compute_splits_cid(tree, cd.cand_tip_bits, cd.cand_buf);
  CidSplitSet& cand = cd.cand_buf;
  double budget = cd.score_budget;

  double mci_sum = 0.0;
  double weight_done = 0.0;
  for (int i = 0; i < cd.n_trees; ++i) {
    if (cd.tree_weights[i] <= 0.0) continue;
    double mci = mutual_clustering_info(cand, cd.tree_splits[i],
                                         cd.n_tips, cd.lap_scratch);
    mci_sum += cd.tree_weights[i] * mci;
    weight_done += cd.tree_weights[i];
    // Early termination: even with perfect MCI for remaining trees,
    // score can't beat budget
    if (budget < HUGE_VAL) {
      double remaining = cd.weight_sum - weight_done;
      // Upper bound on remaining MCI: each tree contributes at most cand_CE
      double cand_ce = clustering_entropy_fast(cand, cd.n_tips, cd.lg2_n);
      double best_possible = -(mci_sum + remaining * cand_ce) / cd.weight_sum;
      if (best_possible > budget) return best_possible;
    }
  }
  return -mci_sum / cd.weight_sum;
}
```

Note: `cand_ce` for the early-termination bound should be computed once
before the loop (not per iteration). Move it before the loop and only
use it in the budget check.

### 2. C++: `src/ts_cid.h` — mark `normalize` vestigial

- Add comment on `normalize` field: `// Vestigial; no longer used by cid_score()`
- Update `cid_score()` doc comment: "Returns negated mean MCI (lower = better consensus)."

### 3. C++: `src/ts_rcpp.cpp` — ignore `normalize` parameter

Keep the `normalize` parameter in the Rcpp signature (backward compat)
but ignore it. No code change needed in the bridge logic since
`cid_data.normalize` is still set but `cid_score()` no longer reads it.

### 4. R: `R/CIDConsensus.R` — remove normalize plumbing

**Functions to modify:**

| Function | Change |
|----------|--------|
| `CIDConsensus()` | Remove `normalize` param; remove rogue-without-normalize warning; update docs |
| `.CIDDrivenSearch()` | Remove `normalize` from args; pass `normalize = FALSE` to C++ (ignored) |
| `.CIDDrivenSearch()` | Negate `result[["best_score"]]` when attaching as `attr(tree, "score")` |
| `.TopologySearch()` | Remove `normalize` from `.NullOr()` call; negate returned score |
| `.CIDScoreFast()` | Remove normalize branch; return `-mean(MCI)` |
| `.MakeCIDData()` | Remove `normalize` parameter |
| `.CollapseRefine()` | Negate `bestScore` on output (internal comparisons unchanged) |
| `.RogueRefine()` | Remove `originalNormalize`; remove normalize propagation to helpers |
| `.PrescreenMarginalNID()` | Remove `normalize` arg from `.MakeCIDData()` calls |
| `.BestInsertion()` | Remove `normalize` arg from `.MakeCIDData()` call |

**Score output convention:**
- C++ returns negated MCI (lower = better). Negate back to positive MCI
  at every `attr(tree, "score") <-` assignment in `.CIDDrivenSearch()`,
  `.TopologySearch()`, `.CollapseRefine()`, `.RogueRefine()`.
- `.CollapseRefine()` and `.RogueRefine()` use `.CIDScorer()` internally
  which already returns negated MCI — so all `<` comparisons work unchanged.
  Only the final `attr(result, "score")` needs negation.
- `.ScoreTree()` returns the internal score (negated MCI) — used for
  internal comparisons only.

**Verbosity messages:** Update messages in `.CollapseRefine()` and
`.RogueRefine()` to negate scores before display (show positive MCI).

**Roxygen updates:**
- `@return` score description: "mean MCI (higher is better)"
- Remove `@param normalize`
- Remove normalize-related `@details`

### 5. Tests: `tests/testthat/test-CIDConsensus.R`

| Test | Change |
|------|--------|
| "CIDConsensus runs end-to-end" (line 155-157) | `score >= 0` still valid (MCI is non-negative) |
| "Collapse should be equal or better" (line 381) | Flip to `>=` (higher MCI = better) |
| "Normalized scoring returns value in [0, 1]" (line 388) | **Delete** (normalize removed) |
| "Normalized scoring differs from raw" (line 397) | **Delete** |
| Rogue tests (lines 461-497) | Remove `normalize = TRUE` |
| maxDrop test (line 506) | Remove `normalize = TRUE` |
| "Warning fires for rogue dropping" (line 516) | **Delete** |
| "CIDConsensus normalized mode" (line 529) | **Delete** |
| Add new test | Verify `attr(result, "score")` ≈ `mean(MutualClusteringInfo(result, trees))` |

### 6. Tests: `tests/testthat/test-ts-cid.R`

| Test | Change |
|------|--------|
| All `normalize = FALSE` args | Remove (default is now to ignore) |
| `result$best_score >= 0` (line 40) | Change to `<= 0` (negated MCI) |
| `attr(coll, "score") <= attr(noColl, "score")` (line 139) | Flip to `>=` |
| "Normalized CID score is in [0, 1]" (line 266) | **Delete** |
| Add new test | Verify negated MCI matches `-mean(MCI)` from R |

### 7. Build and test

- Tarball build to `.agent-build`
- Run `test-ts-cid.R` and `test-CIDConsensus.R`
- Verify CID score matches `mean(MutualClusteringInfo())` from TreeDist

### 8. Commit

Single commit: `refactor: reframe CIDConsensus as MCI maximization`

## Risk

Low. Core algorithm (LAP, MCI, MRP screening) is untouched. Only the
aggregation formula and presentation layer change. All internal `<`
comparisons remain valid via negation.
