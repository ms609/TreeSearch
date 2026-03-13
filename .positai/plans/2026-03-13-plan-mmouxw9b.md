# Plan: Implied Weights (IW) Scoring

## Context

Under equal weights (EW), the tree score is `Σ weight_i × steps_i` — a linear sum.
Under implied weights (IW, Goloboff 1993), each character's contribution is nonlinear:

```
IW score = Σ freq_i × extra_i / (k + extra_i)
```

where `extra_i = steps_i - min_steps_i` and `k` is the concavity constant. Lower score = better tree. When `k = Inf`, this reduces to EW.

The key challenge: the current bit-packed scoring sums 64 characters at once via `popcount`. Under IW, each character feeds into a nonlinear function individually, so we can't use the same trick for final scoring — but we *can* still use bit-packing for state computation and for the indirect TBR calculation (with a marginal-delta trick).

## Architecture

The plan minimizes changes to search modules (TBR, ratchet, drift, sector, fuse) by encapsulating IW logic in the scoring layer. Search modules already use `double` for scores (forward-compat design). The main changes are:

1. **Scoring functions** return `double` and dispatch EW vs IW internally
2. **Indirect calculation** returns `double` under IW (using precomputed marginal deltas)
3. **Search modules** need minimal changes — just switching `int` to `double` in a few local variables

## Steps

### Step 1: DataSet changes

**File: `ts_data.h`**
- Add `double concavity` to `DataSet` (default `HUGE_VAL` = EW mode)

**File: `ts_data.cpp`**
- Accept `min_steps` from R in `build_dataset()` (new parameter: `const int* min_steps_r`)
- Populate `ds.min_steps[p]` from the input array

**File: `ts_rcpp.cpp`**
- Update `make_dataset()` helper to accept `IntegerVector min_steps` and `double concavity`
- Pass through to `build_dataset()`
- Update all callers of `make_dataset()` to pass the new parameters (default: empty min_steps, k=Inf for backward compat)

### Step 2: Per-character step extraction

**File: `ts_fitch.h` / `ts_fitch.cpp`**

Add a function to extract per-character step counts from `local_cost` masks after a downpass:

```cpp
void extract_char_steps(const TreeState& tree, const DataSet& ds,
                        std::vector<int>& char_steps);
```

This iterates over all internal nodes, scanning `local_cost` masks:
- For standard blocks: scan `local_cost[node * n_blocks + b]`
- For NA blocks: recompute from `down2` children (same logic as `ts_na_char_steps` debug function)

The result is `char_steps[pattern_idx]` = total steps for that pattern.

### Step 3: IW scoring function

**File: `ts_fitch.h` / `ts_fitch.cpp`**

```cpp
double score_tree(TreeState& tree, const DataSet& ds);
```

This is the new unified entry point:
- If `ds.concavity` is infinite → call `fitch_score()` (or `fitch_na_score()`), return as double
- If finite → run downpass/uppass (or NA three-pass), then `extract_char_steps()`, then apply IW formula

All search modules will call `score_tree()` instead of `fitch_score()` / `fitch_na_score()`.

### Step 4: IW indirect calculation

**File: `ts_fitch.h` / `ts_fitch.cpp`**

New function for IW-aware indirect evaluation:

```cpp
double indirect_iw_score(
    const uint64_t* clip_prelim,
    const TreeState& tree, const DataSet& ds,
    int node_a, int node_d,
    double base_iw,                    // precomputed IW of divided tree
    const std::vector<double>& iw_delta);  // marginal cost of +1 step per pattern
```

The approach:
1. **Once per clip** (in TBR loop): extract `divided_steps[]` per character, compute `base_iw`, and precompute `iw_delta[p]` = marginal IW cost if pattern `p` gains one more step:
   ```
   iw_delta[p] = freq[p] * ((e+1)/(k+e+1) - e/(k+e))
   where e = divided_steps[p] - min_steps[p]
   ```
2. **Per regraft candidate**: compute `needs_step` mask (same as current `fitch_indirect_length`), then scan set bits and sum `iw_delta[pattern_index[c]]`
3. Return `base_iw + Σ iw_delta` for the candidate

This is efficient: the per-regraft inner loop only touches characters that actually gain a step at that edge, which is typically a small fraction.

For EW mode, continue using the existing `fitch_indirect_length()` — no performance regression.

### Step 5: Update TBR search

**File: `ts_tbr.cpp`**

Changes:
- Replace `full_rescore()` calls with `score_tree()`
- In the clip loop, after incremental downpass:
  - If IW: extract `divided_steps`, compute `base_iw` and `iw_delta[]`
  - Replace `fitch_indirect_length()` calls with `indirect_iw_score()`
- Score comparisons already use `double`; add small epsilon for float equality

The `full_rescore` helper becomes:
```cpp
static double full_rescore(TreeState& tree, const DataSet& ds) {
  tree.reset_states(ds);
  return score_tree(tree, ds);
}
```

### Step 6: Update other search modules

**Ratchet** (`ts_ratchet.cpp`): calls `tbr_search()` internally. The TBR changes propagate. The `RatchetResult::best_score` is already `double`. The ratchet perturbs `active_mask` — this interacts correctly with IW since masked-out characters contribute zero to both EW and IW scores.

**Drift** (`ts_drift.cpp`): Has its own TBR loop copy. Needs the same changes as TBR (Step 5). The AFD/RFD calculations use step deltas — under IW, these become IW-score deltas rather than raw step deltas. The RFD computation using `local_cost` diffs needs adjustment: under IW, compute RFD from per-character fit changes rather than per-node step diffs.

**Sector** (`ts_sector.cpp`): Builds a reduced dataset and runs TBR on it. The reduced dataset inherits `concavity` and `min_steps`. After reinsertion, the full tree is rescored. The sector's `min_steps` for the reduced problem is an approximation — we can use the full-tree `min_steps` for the subset of characters in the sector, which is a valid lower bound.

**Fuse** (`ts_fuse.cpp`): Compares two full trees. Uses `fitch_score()` — switch to `score_tree()`.

### Step 7: Update R interface

**File: `ts_rcpp.cpp`**

- Update `ts_fitch_score` to accept optional `concavity` and `min_steps` parameters
- Update `ts_tbr_search`, `ts_ratchet_search`, `ts_drift_search`, `ts_rss_search`, `ts_xss_search` to accept `concavity` and `min_steps`
- The existing `morphy_iw()` R function provides the reference implementation for verification

### Step 8: Tests

**File: `tests/testthat/test-ts-iw.R`**

- Compare `ts_iw_score()` against `morphy_iw()` on the same trees for multiple datasets
- Test with multiple concavity values (k = 2, 3, 10, 100)
- Test with inapplicable datasets (IW + NA interaction)
- Test that IW TBR search finds trees at least as good as EW search on IW criterion
- Regression: verify EW scoring is unchanged when concavity = Inf

## Key Design Decisions

1. **Unified `score_tree()` entry point** rather than making callers branch. Keeps search modules clean.

2. **Marginal-delta trick for indirect calc** rather than per-character loop per candidate. Preserves the speed advantage of the indirect approach under IW.

3. **`min_steps` passed from R** rather than computed in C++. The existing `MinimumLength()` R function handles this correctly including inapplicable characters.

4. **No changes to bit-packing or state computation**. The Fitch downpass/uppass/NA algorithm is identical under IW — only the final score aggregation changes.

5. **Float comparison with epsilon** for score equality checks.

## Risk Assessment

- **Drift RFD**: The current RFD calculation uses per-node step diffs from `local_cost`. Under IW, the analogous quantity would be per-node *fit* diffs. May need a simplified approach (e.g., use score delta directly). Low risk — drift is a perturbation heuristic.

- **Sector reduced datasets**: The `min_steps` for the reduced dataset should be the full-tree `min_steps` for the included characters (a lower bound). Medium risk — needs careful validation.

- **NA + IW interaction**: The three-pass algorithm's step counts feed into IW the same way as standard Fitch steps. Low risk — already demonstrated in `ts_na_char_steps`.

## Estimated Scope

- ~6 files modified (ts_data.h/cpp, ts_fitch.h/cpp, ts_tbr.cpp, ts_rcpp.cpp)
- ~2 files with lighter changes (ts_drift.cpp, ts_fuse.cpp)
- 1 new test file (test-ts-iw.R)
- Total: ~300-400 lines of new/modified C++ code
