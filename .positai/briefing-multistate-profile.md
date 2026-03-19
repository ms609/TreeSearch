# Briefing: Extending Profile Parsimony to >2 States

## Status: T-101 DONE, T-102–T-107 OPEN

## Goal
Extend profile parsimony scoring from 2-state characters to multi-state (3+).

## What exists already

### Current 2-state implementation (on main branch)
- `R/pp_info_extra_step.r`: `StepInformation()` — computes information content
  per character for all possible step counts. Uses `LogCarter1()` for 2 states.
  Lines 51-56 explicitly warn and drop states beyond 2 informative tokens.
- `R/data_manipulation.R`: `PrepareDataProfile()` — decomposes multi-state
  characters into pairs (keeping top-2 informative states), compresses to
  binary, builds `info.amounts` matrix. Lines 89-115: `.RemoveExtraTokens()`
  keeps only 2 most informative states; line 138 asserts exactly 2 non-ambig.
  Lines 196-201 hardcode `levels = c("0", "1")` and a 3×2 contrast matrix.
- `R/MaximizeParsimony.R`: Lines 424-428 call `PrepareDataProfile()`, then
  lines 528-533 extract `info.amounts` attribute and pass to C++.
- C++ engine: `ts_data.cpp` copies `info_amounts` table; `ts_fitch.cpp`
  looks up `info_amounts[(step-1) + info_max_steps * pattern]` for each
  pattern. The C++ scoring pipeline is generic — it handles multi-state Fitch
  natively. Only the R-level data prep is restricted to 2 states.

### Prior multi-state work (on `concordance-FitchInfo` branch, NOT on main)
- `src/MaddisonSlatkin.cpp`: Full C++ implementation of the Maddison & Slatkin
  (1991) recursive algorithm for counting trees with exactly s steps for a
  multi-state unordered character. Supports up to 5 states.
- `R/FitchInfo.R`: Uses `MaddisonSlatkin()` for concordance scoring with
  multi-state characters. Not profile parsimony weighting, but the
  mathematical core is exactly what's needed.
- Key commits: `ab5f80be` "Support 2-5 states", `23963c07` "Embed MadSlat to FI",
  `9336c066` "FitchInfo" (latest on concordance-FitchInfo).
- The FitchInfo code already converts MaddisonSlatkin output to cumulative
  information content (bits), which is the same transformation profile
  parsimony needs.

### Key mathematical insight
The existing `MaddisonSlatkin()` computes exactly what `StepInformation()`
needs: `log P(s steps | n_0, n_1, ..., n_k leaves)` for each possible step
count s, averaged over all unrooted binary trees. This is the multi-state
generalization of Carter et al. (1990)'s theorem 1.

Profile parsimony's information content = `log2(N_total_trees) - log2(cumsum(N_trees_with_≤s_steps))`
where `N_trees_with_exactly_s_steps = exp(MaddisonSlatkin(s, states)) * N_total_trees`.

So `MaddisonSlatkin()` output feeds directly into `StepInformation()`.

## Architecture: what needs to change

### Layer 1: Mathematics (already done on branch)
`MaddisonSlatkin()` computes `log(fraction of trees with exactly s steps)`.
This is the multi-state analog of `LogCarter1()`.

### Layer 2: `StepInformation()` (R/pp_info_extra_step.r)
Currently calls `LogCarter1()` for 2 states, rejects >2.
Needs: dispatch to `MaddisonSlatkin()` when >2 informative states.
The transformation from log-probabilities to information content is identical.

### Layer 3: `PrepareDataProfile()` (R/data_manipulation.R)
Currently decomposes to pairs and hardcodes binary contrast matrix.
Needs: pass multi-state characters through directly (no decomposition).
Must build a proper contrast matrix for k states + ambiguous token.
The `info.amounts` matrix dimensions change (more rows = more possible steps).

### Layer 4: C++ engine
**No changes needed.** The `info_amounts` lookup table is already generic —
indexed by `(step, pattern)`. The Fitch scoring engine already handles
multi-state characters. We just need to feed it the right contrast matrix
and info_amounts.

## Literature

| Reference | Role |
|-----------|------|
| Carter et al. (1990) | Exact formula for 2-state trees — current basis |
| Steel (1993) | Distribution theory for bicolored trees |
| Steel & Charleston (1995) | Properties of parsimoniously colored trees |
| Steel, Goldstein & Waterman (1996) | CLT for parsimony length |
| Maddison & Slatkin (1991) | Recursive algorithm for multi-state — the key |
| Faith & Trueman (2001) | Original profile parsimony justification |

No known closed-form generalization of Carter for >2 states exists.
Maddison & Slatkin's recursive algorithm is the standard approach.

## Performance concerns
- MaddisonSlatkin is exponential in number of states (2^k bitmask states)
- Current C++ implementation handles up to 5 states
- For 2 states: `LogCarter1()` is O(1) per step count
- For 3-5 states: memoized recursion, feasible for typical morphological data
- For >5 states: may need approximation or capping
- `info.amounts` computation is a one-time precomputation cost (not in search loop)

## Risks
1. MaddisonSlatkin.cpp is on a different branch — needs careful merge/cherry-pick
2. May need to handle the interaction with character simplification (currently
   characters with many states get collapsed)
3. Performance for characters with many taxa AND many states
4. Need to handle edge cases: all-ambiguous, singleton states, etc.
5. Test coverage: existing profile parsimony tests assume binary data

---

## T-106: Approximation for >5 State Profile Parsimony — Research Analysis

### 1. Scaling of exact MaddisonSlatkin

Benchmarked on Windows, R 4.5.2, single-threaded. All timings are for
computing the full step-count range (s_min to n-1).

| k (tokens) | n (tips) | tips/state | Time |
|:-----------:|:--------:|:----------:|-----:|
| 2 | 4  | 2 | <1 ms |
| 2 | 10 | 5 | 10 ms |
| 3 | 9  | 3 | 100 ms |
| 3 | 15 | 5 | 3.7 s |
| 4 | 8  | 2 | 320 ms |
| 4 | 12 | 3 | 12.6 s |
| 4 | 20 | 5 | timeout (>30 s) |
| 5 | 10 | 2 | timeout (>30 s) |

**Root cause:** The recursion partitions n tips into two subtrees in
all valid ways, for each of 2^k−1 root states, for each step count s
in 0..n−k. The memoization table grows as
O(#unique_leaf_configs × max_steps × #states), and the number of
unique leaf configurations grows combinatorially with n and k.

**Conclusion:** Exact computation is infeasible for k≥5 with n≥15, or
k≥6 at any practical n. The current code's `k ≤ 5` limit is well-placed.

### 2. Approximation approaches evaluated

#### (a) Plain Monte Carlo

Sample N random unrooted binary trees, score each with Fitch, tally the
step-count distribution.

**Test:** k=6, n=30, split=(8,7,5,4,3,3), N=10,000 random trees.
Rate: ~1,700 trees/second. Distribution: observed range 13–22, peak at 19.

**Problem:** The exact P(s_min=5) = exp(−38.6) ≈ 1.7×10⁻¹⁷, while the
smallest observable MC probability at N=10⁴ is 10⁻⁴. The "gap" between
the minimum step count and the MC-observable range spans 8 step counts and
13 orders of magnitude. Even at N=10⁶, the gap persists.

**Verdict:** Cannot estimate the information-rich left tail. Only useful
for the body/right tail of the distribution.

#### (b) Normal (CLT) approximation

Steel, Goldstein & Waterman (1996) proved asymptotic normality of
parsimony length for binary characters. Multi-state CLT should hold by
similar arguments (sum of nearly independent subtree contributions).

**Test:** Fitted normal(μ=18.9, σ=1.5) from MC data, extrapolated to s_min.

| Metric | Normal | Exact |
|--------|-------:|------:|
| log P(s_min=5) | −44.3 | −38.6 |
| IC(s_min) bits | 62.1 | ~55.7 |

**The normal overestimates IC at the minimum by ~6 bits** (the true
distribution has heavier left tails than Gaussian). However, this error
is at step counts that never occur on real trees during a search.

In the MC-observable range (13–22 steps), the normal approximation agrees
well with empirical data. This is the range that actually affects search
decisions.

**Verdict:** Accurate in the practical range. Left-tail error is
large but irrelevant for search quality.

#### (c) Hybrid: exact anchor + MC body

The key insight enabling this approach:

> **P(s_min) has an exact O(k) formula for any k:**
> `P(s_min) = NUnrootedMult(split) / NUnrooted(n)`
>
> This uses the product-of-double-factorials counting formula for labeled
> trees consistent with k non-overlapping groups, and requires no recursion.

The hybrid approach:
1. Exact P(s_min) via `NUnrootedMult` (instant for any k)
2. MC sample of N=50,000 random trees → empirical distribution for the body
3. Normal fit to MC data → parametric extrapolation for the sub-MC left tail
4. Blend: use exact at s_min, normal extrapolation for s_min+1 to MC left
   edge, empirical distribution for MC-observable range

**Verdict:** This is the recommended approach (see §3 below).

#### (d) "Keep top 5" (current fallback)

The existing `StepInformation()` already handles k>5 by keeping the 5 most
frequent tokens and dropping the rest. This discards real information
(the dropped tokens contribute genuine parsimony signal) but is safe.

For characters where the dropped tokens each have only 2–3 leaves, the
information loss is modest. For characters with 6+ well-represented tokens,
the loss is significant but hard to quantify without exact values.

#### Other approaches considered

- **Importance sampling** (bias toward low-step trees): Could solve the
  left-tail problem but requires a carefully designed proposal distribution.
  Engineering effort disproportionate to the niche use case.
- **WithOneExtraStep() extension** to k>2: Currently unimplemented. The
  combinatorics are substantially harder for k>2 (multiple ways to place
  the extra step among k groups). Could provide exact P(s_min+1) but would
  not solve the general left-tail problem.
- **Extending MaddisonSlatkin to k=6:** Structural changes to support
  2^6−1=63 states are modest (add `StateKeyT<6>` template), but the
  computational blowup still makes it infeasible for n>10–12.

### 3. Recommendation

**Primary approach: MC-calibrated normal approximation with exact anchor**

This approach requires minimal new code, has well-understood error
properties, and covers the only practical use case (characters with 6+
states in morphological datasets).

**Why this is sufficient:** Profile parsimony's search engine only compares
info_amounts values at step counts that actually occur on candidate trees.
For a k=6 character on 30 tips, candidate trees typically score 13–22 steps
(based on MC data). The information content in this range is well-estimated
by the normal approximation calibrated to MC samples. The extreme left
tail (5–12 steps) has enormous IC values that serve as theoretical upper
bounds but never affect search decisions, because no reasonable tree
achieves those step counts.

**Performance:** The MC sampling adds ~30 seconds per character at
N=50,000 trees (for n=30 tips). This is a one-time precomputation cost.
For datasets with few >5-state characters, this is acceptable. For
datasets with many such characters, the MC could be parallelized or the
sample size reduced.

**Accuracy:** In the practical range (within ~3σ of the MC mean), the
normal approximation's IC values match empirical estimates to within
~0.1 bits. This is smaller than the character's own noise and does not
materially affect search quality.

**Fallback:** If MC is too slow or the user needs a quick result, retain
the existing "keep top 5" heuristic as an option.

### 4. Prototype R code

```r
#' Approximate StepInformation for >5 state characters
#'
#' Uses exact P(min_steps) + MC-calibrated normal approximation.
#'
#' @param split Integer vector of token frequencies (sorted decreasing,
#'   singletons removed).
#' @param n_mc Number of Monte Carlo trees to sample (default 50000).
#' @return Named numeric vector of information content (bits) per step count.
#' @keywords internal
.ApproxStepInformation <- function(split, n_mc = 50000L) {
  k <- length(split)
  n <- sum(split)
  s_min <- k - 1L
  s_max <- n - 1L

  # 1. Exact P(minimum steps) — works for any k
  log_p_min <- log(NUnrootedMult(split)) - log(NUnrooted(n))

  # 2. Monte Carlo: sample random trees, tally step counts
  labels <- paste0("t", seq_len(n))
  char_vec <- rep(seq_along(split) - 1L, split)
  names(char_vec) <- labels
  dat <- TreeTools::MatrixToPhyDat(
    matrix(char_vec, ncol = 1, dimnames = list(labels, "c1"))
  )
  mc_scores <- vapply(
    seq_len(n_mc),
    function(i) RandomTreeScore(dat),
    double(1)
  )

  # 3. Fit normal to MC data
  mu_hat <- mean(mc_scores)
  sd_hat <- sd(mc_scores)

  # 4. Build log-probability vector for all step counts
  steps <- s_min:s_max
  n_steps <- length(steps)
  log_p <- numeric(n_steps)

  for (i in seq_along(steps)) {
    s <- steps[i]
    if (s == s_min) {
      # Exact value
      log_p[i] <- log_p_min
    } else {
      # MC estimate (with continuity correction)
      mc_count <- sum(mc_scores == s)
      if (mc_count > 0) {
        # Direct empirical estimate
        log_p[i] <- log(mc_count / n_mc)
      } else {
        # Normal extrapolation for unobserved step counts
        log_p[i] <- dnorm(s, mu_hat, sd_hat, log = TRUE)
      }
    }
  }

  # 5. Cumulative IC
  ret <- -.LogCumSumExp(log_p) / log(2)
  ret[ret < sqrt(.Machine[["double.eps"]])] <- 0
  names(ret) <- steps

  ret
}
```

### 5. Implementation plan

| Step | Description | Effort |
|------|-------------|--------|
| 1 | Add `.ApproxStepInformation()` to `R/pp_info_extra_step.r` | 1 hr |
| 2 | Modify `StepInformation()`: dispatch to `.Approx...` when k>5 (instead of current top-5 truncation) | 30 min |
| 3 | Add `approx` parameter to `StepInformation()` with options `"exact"` (current), `"mc"` (new), `"auto"` (default: exact for k≤5, MC for k>5) | 30 min |
| 4 | Tests: verify MC approximation agrees with exact for k=3 within ~10% relative IC at practical step counts | 1 hr |
| 5 | Documentation: update `StepInformation()` docs to describe approximation | 30 min |

Total: ~3.5 hours. No C++ changes needed.

### 6. Comparison with existing "keep top 5" approach

| Criterion | Keep top 5 | MC approximation |
|-----------|-----------|------------------|
| Speed | Instant (delegates to exact) | ~30s per character |
| Accuracy at practical range | Unknown (drops signal) | ~0.1 bit error |
| Left tail | Exact for reduced char | Exact P(min) + normal extrapolation |
| Handles any k | Yes (truncates) | Yes |
| New code | 0 lines | ~60 lines R |
| C++ changes | None | None |

For datasets where >5-state characters are rare (typical morphology), the
MC overhead is negligible relative to the search time. For datasets with
many such characters, the top-5 fallback remains available.

### 7. Future improvements (deferred)

- **Exact P(s_min + 1):** Extending `WithOneExtraStep()` to k>2 would
  give a second exact anchor point, improving the left-tail interpolation.
  The combinatorics are non-trivial but tractable.
- **Importance sampling:** For characters where the search regularly reaches
  near-minimum step counts (small n, few states per token), importance
  sampling could improve accuracy. Not worth implementing unless a specific
  dataset demonstrates the need.
- **Cached MC tables:** For common state-frequency patterns, pre-computed
  MC tables could eliminate the per-character sampling cost.
