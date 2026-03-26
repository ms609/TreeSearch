# Incremental CID Scoring — Feasibility Analysis

Date: 2026-03-26

## Problem

CID verification is the dominant cost in `InfoConsensus()`. Each accepted
TBR move triggers a full `cid_score()` call, which computes
`mutual_clustering_info()` for every input tree. Per-candidate cost is
O(T × n^2.4) due to the Jonker-Volgenant LAP solver.

Measured per-candidate costs (from cid-scaling benchmarks):

| n_tip | T=100 | T=1000 |
|:-----:|:-----:|:------:|
| 50 | 14 ms | 142 ms |
| 100 | 71 ms | 706 ms |
| 200 | 405 ms | 4,031 ms |

## Current Pipeline

```
TBR candidate evaluation (all O(n^2) candidates):
  → MRP/Fitch indirect scoring (fast, ~0.02 ms/candidate)
  → Best MRP candidate selected

If best MRP candidate not dominated:
  → Apply move to topology
  → cid_score() — FULL rescore (expensive):
      for each input tree i in 1..T:
        compute_splits_cid(candidate)           # O(n)
        mutual_clustering_info(cand, tree[i]):  # O(k^2.4)
          Phase 1: exact match via hash         # O(n) amortised
          Phase 2: LAP for unmatched splits     # O(k^3) where k = unmatched
      return -mean(MCI)
```

Only 1 CID call per accepted TBR move (not per candidate), so the issue
is cost per call, not number of calls.

## Incremental Approaches Considered

### A. Incremental Split Computation

**Idea:** After a TBR move, only ~4 splits change. Maintain the split set
incrementally instead of recomputing from scratch.

**Cost of `compute_splits_cid()` currently:** O(n) — negligible compared
to LAP. At 200 tips, split computation is ~microseconds; LAP is ~4 ms
per tree. Not worth optimizing.

### B. Per-Tree MCI Caching with Invalidation

**Idea:** Cache `mutual_clustering_info()` result per input tree. After a
TBR move, only recompute for trees whose matching is affected by the
changed splits.

**Problem:** A TBR move changes ~4 candidate splits, which potentially
affects the matching with ALL input trees. Whether the matching is
actually affected depends on the optimal LAP assignment, which we'd need
to check per tree — i.e., re-examine the matching anyway.

**Partial optimisation:** If a changed candidate split was not part of any
LAP matching (it was exactly matched), we only need to check whether its
replacement is also exactly matchable. If exact_match(old_split) and
exact_match(new_split), the MCI for that tree is unchanged except for the
local contribution of those two splits. This can be computed in O(1).

**Benefit:** Mainly helps when input and candidate trees are similar (many
exact matches), which is exactly the regime where subsampling already
makes T small. Net benefit is modest.

### C. Warm-Start LAP Solver

**Idea:** After a TBR move, the LAP cost matrix changes by ~4
rows/columns. Use the previous optimal assignment as a warm start and
re-solve only the affected portion.

**JV algorithm warm-start:** The Jonker-Volgenant SAP algorithm maintains
dual variables (u[], v[]) and a column assignment. In principle, one could
invalidate only the rows/columns corresponding to changed splits and
re-augment from there.

**Complexity:** Implementation is non-trivial. The JV algorithm's
augmenting-path phase assumes fresh initialization. Partial re-solve
requires careful bookkeeping of which rows/columns are "dirty" and
re-running augmentation only for those.

**Expected speedup:** At most 4× if 4 of k splits changed and the
re-solve is O(4 × k^2) instead of O(k^3). In practice, k (unmatched
splits) is often small (5-20), so the LAP is already fast. The constant
factor of implementing warm-start may negate the asymptotic benefit.

**Verdict:** High implementation cost, moderate benefit. Worth pursuing
only if the simpler approaches prove insufficient.

### D. Early Termination on Input Trees (already implemented)

The current code sets `score_budget = best_score + eps` before CID
rescoring. If accumulated MCI can't beat the budget (using a per-tree
upper bound from `clustering_entropy_fast()`), it returns early.

**Measured benefit:** When the TBR move is actually an improvement, the
budget is the old score — so all trees are evaluated. When the move is
not an improvement, early termination kicks in after ~T/2 trees on
average (depends on how bad the move is). Since the dominated check
already filters most bad moves, the remaining CID calls are mostly
improvements, so early termination helps less often.

**Enhancement:** Sort input trees by "expected discriminating power"
(e.g., similarity to the current candidate). Evaluate the most
discriminating trees first so early termination fires sooner for
non-improving moves.

### E. Batch Evaluation of Top-k Candidates

**Idea:** Instead of MRP-selecting the single best candidate, keep the
top-k MRP candidates (k=3-5) and CID-score all of them. Take the best
CID result.

**Current code:** Only the best MRP candidate gets CID-scored. If MRP
ranking disagrees with CID ranking (which it often must — they're
different objectives), we might miss CID-improving moves.

**Cost:** k× more CID calls per TBR pass. At 100 tips with T_sub=200
(auto), each CID call is ~14 ms. With k=5, the additional cost is ~56 ms
per TBR pass — probably acceptable.

**Benefit:** Better CID exploration. The MRP screening is a coarse proxy;
evaluating more candidates should find better CID moves per pass,
potentially converging in fewer TBR passes.

**Verdict:** Medium implementation cost, potentially high benefit for
search quality (not raw per-call speed). Worth benchmarking.

## Recommendation

Given the treeSample = "auto" subsampling already reduces T dramatically
(1000→100-200) and makes per-call CID cost manageable:

**Priority 1 (low effort, potential high impact):**
- **Approach E**: Batch top-k candidate evaluation. Simple to implement in
  `ts_tbr.cpp` — collect top-k MRP candidates instead of just the best,
  CID-score each. Expected to improve search quality.

**Priority 2 (medium effort, moderate impact):**
- **Approach D enhancement**: Sort input trees by discriminating power.
  Requires precomputing a "difficulty" score per tree during
  `prepare_cid_data()`. Modest code change.

**Priority 3 (high effort, uncertain impact):**
- **Approach C**: Warm-start LAP. Only worth pursuing if 200+ tip trees
  become a primary use case and treeSample can't reduce T enough.

**Not recommended:**
- **Approach A**: Split computation is not on the critical path.
- **Approach B**: Complexity doesn't justify benefit given subsampling.

## Quantitative Impact Estimates

For a typical TBR convergence at 100 tips with treeSample="auto" (T_sub=200):

| Metric | Current | With top-5 | With sorted trees |
|--------|:-------:|:----------:|:-----------------:|
| CID calls/pass | 1 | 5 | 1 |
| Time/CID call | 14 ms | 14 ms | 14 ms (avg) |
| CID time/pass | 14 ms | 70 ms | 14 ms (same, but earlier termination) |
| TBR passes to converge | ~50 | ~30 (better moves) | ~50 |
| Total CID time | 700 ms | 2,100 ms | 500 ms (est.) |
| Total replicate time | ~3 s | ~5 s | ~2.5 s |

The "fewer TBR passes" estimate for top-5 is speculative; actual benefit
depends on how often MRP and CID disagree on the best move.
