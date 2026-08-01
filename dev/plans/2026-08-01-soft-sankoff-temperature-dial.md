# Soft-Sankoff: a temperature dial between parsimony and likelihood

**Date:** 2026-08-01
**Status:** exploration. **Gate B CLOSED and FAILED (2026-08-01): Steps 4 and 5
are dead on cost.** Gate A retired rather than answered, since it existed only to
protect Step 4. Steps 3a and 3b were placed off the gates and are live. **3a is
DONE** over 100 matrices in two independent runs: a genuine interior optimum at
`T = 0.5` (robust, p = 0.002 in both) and soft-Sankoff ranking *within* the MPT
set (quantile 0.33 vs a null of 0.5, p ~ 1e-6, n = 92); homoplasy tracking **not**
supported. **But NEITHER result transfers to 75-tip O'Reilly matrices** — the rank
statistic is chance-level there at every temperature (min p = 0.118, n = 39), with
a U-shaped distribution that looks like the density tilt Gate A predicted. Treat
the Step 3a findings as properties of `congreveLamsdellMatrices`, not of the
criterion. A prototype scorer now exists in
`src/ts_soft_sankoff.{h,cpp}` — new
files only, reachable from no default scoring path, and deliberately not sharing
a struct with `ts_sankoff.{h,cpp}`, which sits on the live x-transformation
pathway with its open T-374/T-385 rooting defects.
**Branch:** `feature/soft-sankoff` (worktree `../worktrees/TS-softsankoff`, from `cpp-search` @ `f13f0c31`)
**Origin:** discussion of Siepel, Hassett & Staklinski, *VINE: Variational inference for scalable Bayesian reconstruction of species and cell-lineage phylogenies*, bioRxiv `10.64898/2025.12.24.696405`. Vine's chain rule terminates at the branch-length gradient, which parsimony does not have; softening the Sankoff `min` is the construction that supplies one.
**Reference implementation:** [`tests/testthat/helper-soft-sankoff.R`](../../tests/testthat/helper-soft-sankoff.R) — pure R, no build required.
**Evidence scripts:** [`dev/soft-sankoff/`](../soft-sankoff/)

---

## The claim, up front

Sankoff's dynamic program with the `min` replaced by a soft-min,

```
softmin_T(x)  =  -T * log( sum_j exp(-x_j / T) )
S~_v(i)       =  sum over children c of  softmin_T over j of [ cost(i,j) + S~_c(j) ]
```

is **exactly** weighted parsimony at `T -> 0` and **exactly** Felsenstein's pruning
algorithm at `T = 1` with `cost(i,j) = -log P_ij(t)`. Substituting
`L_v(i) = exp(-S~_v(i))` into the `T = 1` recursion gives

```
L_v(i) = prod over children c of  sum_j P_ij(t_c) * L_c(j)
```

which is pruning, term for term. Parsimony and likelihood are the same dynamic
program in two semirings — tropical `(min, +)` and probability `(sum, x)` — and
`T` is the dial between them. This is the tropicalisation result of Pachter &
Sturmfels; it is distinct from the Tuffley & Steel (1997) bridge, which reaches
parsimony from likelihood via free per-character-per-branch rates rather than via
temperature.

The interpretive consequence is the payload. Because
`softmin_T(x) ~= min(x) - T * log(number of near-minima)`, **temperature controls
how much the criterion integrates over ancestral reconstructions rather than
optimising them.** Hard Sankoff commits to the single best reconstruction;
Felsenstein sums over all of them; `T` in `(0, 1)` is principled partial
integration.

## Why this is cheaper than it sounds

Three facts about the existing codebase, all verified:

- `src/ts_sankoff.{h,cpp}` already implements a general Sankoff DP over
  double-precision cost matrices, with an up-pass, and already coexists with
  Fitch on a per-character split (`score_tree()` returns
  `fitch_score_ew(...) + sankoff_score(...)`). Adding a temperature is a local
  change to code that exists, not a new kernel.
- `ts_sankoff_test()` is already exposed to R (`R/RcppExports.R:229`), so a
  hard-Sankoff oracle is available from R without new bindings.
- `TreeSearch(..., TreeScorer = )` / `EdgeListSearch()` accept a pluggable
  scorer, with a "custom optimality criteria" vignette. A soft scorer can be
  driven through the existing search skeleton in pure R with no C++ at all.

So Steps 0–3 below need **no compiled change**. Only Steps 4 and 5 do, and both
are gated.

---

## The two gates

Everything expensive downstream is conditional on these. Both are cheap, and
both can fail.

### Gate A — does the tilt point the right way?

`softmin` rewards topologies admitting **more** near-optimal reconstructions, and
reconstruction ambiguity is precisely the property that inflates MPT sets. So
smoothing may preferentially walk the search *toward* the broad plateaux it was
meant to escape. This is a signed prediction, not a generic uncertainty.

`WideSample()`'s own documentation states the underlying distinction: density in
the parsimony landscape is not support. Soft-Sankoff scores something closer to
density. Gate A asks whether it nonetheless correlates with the thing we want.

**Test:** on `congreveLamsdellMatrices` (generating tree known, in
`data/referenceTree.RData`), collect a large MPT set, then rank it by the soft
score at a grid of `T`. Measure Spearman correlation between soft-score rank and
(i) RF/CID distance to the generating tree, (ii) Mk likelihood under `phangorn`.

**Pass:** correlation is materially negative (better soft score => closer to
truth) at some `T`, consistently across matrices.
**Fail:** correlation is null or positive. Steps 4 and 5 die; Steps 2 and 3
survive on their own merits.

**First read (2026-08-01, `02-tilt-direction.R`, 6 matrices): not yet decidable,
and the design needs widening.**

- **Only 1 of 6 matrices produced a usable MPT set** (>= 8 distinct trees);
  the rest returned 1–7. The Congreve & Lamsdell matrices are simply not
  MPT-rich enough at 22 tips for a rank correlation over strictly-optimal trees.
  **Design correction:** score trees within a few steps of optimal rather than
  only the MPTs. That is better anyway — it supplies a real score gradient
  instead of a set that is tied by construction, and it is closer to what an
  annealed search would traverse.
- On the one usable matrix the two measures **disagreed**: `rho(CID) = -0.30`
  (soft score points toward the truth) but `rho(Mk) = +0.14` (soft score points
  *away* from the likelihood). With `n = 1` neither number means anything yet,
  but the disagreement is the shape to watch: it is what a density tilt would
  look like if topological closeness and likelihood happen to diverge.
- The ranking was **invariant in `T` across 0.02–0.5**, changing only at
  `T = 1`. If that holds up, the tilt is real but has almost no tuneable range,
  which weakens the case for temperature as a *search* control specifically.

The script's verdict logic now refuses to declare a pass on fewer than 5 usable
matrices, and reports a SPLIT rather than a PASS when the topology and
likelihood measures disagree.

#### GATE A: RETIRED, NOT ANSWERED — 2026-08-01

Gate A asked whether a soft-guided **search** tilts toward truth or toward MPT
density. That question existed solely to protect Step 4, and Gate B killed Step 4
on cost. Widening Gate A across the remaining 94 matrices would spend the compute
to inform nothing, so it was **not** widened. `02-tilt-direction.R` is left
untouched as the record of what the gate said; its pass/fail/SPLIT verdict logic
encodes Step 4 semantics and should not be carried forward.

Two of its first-read observations change meaning rather than standing:

- The **`T`-invariance across 0.02–0.5** was read as a defect ("almost no
  tuneable range for a search control"). With search out of the picture that is
  simply a *result about the criterion*, not a problem.
- The **`ρ(CID)` / `ρ(Mk)` disagreement** at `n = 1` was a warning sign for a
  search. For the criterion question the honest primary axis is distance to the
  **generating tree**, which is known for these matrices; Mk likelihood drops to
  a secondary descriptive rather than a pass/fail arm.

The criterion question that survives is Step 3a, now implemented separately as
`dev/soft-sankoff/04-dial-study.R`.

### Gate B — what is the per-score cost against SIMD Fitch?

No number is guessed here on purpose; it must be measured. Two distinct costs,
and the second is likely the larger:

1. Loss of bit-parallelism: `O(k^2)` doubles plus `exp`/`log` per character per
   node, against word-packed bitwise Fitch (`src/ts_simd.h`).
2. Loss of incremental rescoring. `src/ts_sankoff.h` states the implementation is
   *full-rescore only, no incremental variant*. TBR draws much of its throughput
   from rescoring only the affected subtree; losing that changes the complexity
   of the inner loop, not merely a constant.

**Pass threshold:** to be set by the measurement, but as an orientation — if a
soft score costs more than ~50x a Fitch score, an annealed search cannot pay for
itself against the ratchet, and Step 4 is dead regardless of Gate A.

**First read (2026-08-01, `03-speed-budget.R`): unfavourable.** The
implementation-independent operation-count ratio is **x350** for binary
characters (`nNode * nChars * (k^2 + k + 1)` against `nNode * ceil(nChars/64)`
bitwise word operations), flat across a 16–128 tip ladder. That is 7x the
orientation threshold before any of the following are counted:

- `exp` is not one operation. The compiled ratio is likely to be *worse* than
  x350, not better, unless a fast vectorised `exp` carries most of it.
- the incremental-rescore loss is still entirely unmeasured, and is the term
  that changes complexity rather than a constant.

Pure-R wall ratio was x253–x1015, which is an upper bound only (pure R against
compiled SIMD) and is not the gate.

#### GATE B: CLOSED, FAILED — compiled, 2026-08-01

The C++ prototype (`src/ts_soft_sankoff.{h,cpp}`, Step 4 of the open items) was
built and Gate B re-measured with a compiled kernel on both sides. **The
constant term is now measured, and it is ~3x worse than the op count predicted
— exactly the direction the `exp` caveat above anticipated.**

| Comparison | Median | Range |
|---|---|---|
| compiled soft / Fitch full rescore, **k = 2** | **x1147** | x886–x1232 |
| compiled soft / Fitch full rescore, k = 4 | x2981 | x2500–x3889 |
| prior op-count estimate | x879 | x438–x1313 |
| pure R / compiled soft (speedup the prototype bought) | x100 | x43–x243 |

All 8 cells clear the x50 orientation threshold, by 20–70x. **The k = 2 number
is the headline** for real morphological data: Fitch's cost is in words and so
near-flat in `k`, while the soft kernel pays `k^2`, so the k = 4 rows say
"multistate is worse", not "this is what a matrix costs".

**The complexity term remains OPEN, but is now bounded.** A full-rescore
prototype cannot measure an incremental soft kernel that does not exist. What it
can measure is how much Fitch gains from incrementality — full-rescore cost
against per-candidate incremental cost, both reported from inside C++ by
`ts_bench_tbr_phases()`. That is a **median x386**, and it is the multiplier a
full-rescore soft kernel forfeits *on top of* the constant above. This bounds
the term's size; it does not measure an incremental soft implementation. It is
also moot: the constant alone already fails.

Measurement corrections worth keeping, because each was an error in the
direction that would have flattered the soft kernel:

- The Fitch baseline must be `time_full_rescore_us` from
  `ts_bench_tbr_phases()`, clocked inside C++ with dataset construction
  excluded. Timing `TreeLength()` — or even `ts_fitch_score()` — from R folds
  per-call marshalling of the contrast and tip-data matrices into the
  denominator, inflating it and so understating the ratio.
- That clock counts **whole microseconds**, and a Fitch rescore of 100 patterns
  is 1–7 of them. One run reported `0 us` and an infinite ratio. Both kernels
  are linear in pattern count, so the comparison moved to 2000 patterns, where
  the denominator is 10–94 ticks.
- The soft side is timed by differencing the binding's own `n_rep`, chosen per
  cell from a pilot pass, so R-side marshalling drops out rather than being
  assumed small. Its residual quantisation (`system.time`, ~15 ms, at the low
  `n_rep` the expensive cells need) is 10–20% — immaterial against a conclusion
  clearing its threshold by 20–70x.

**Consequence: Step 4 (annealed search) is DEAD on cost, and Step 5 inherits
that.** Steps 3a and 3b were placed off the gates and survive. Anyone reviving
Step 4 needs a mitigation — SIMD over the character axis, a `k = 2` special
case, or restricting soft scoring to a coarse outer loop — not just a better
`exp`, and would need to find roughly three orders of magnitude, not one.

---

## Stepped plan

### Step 0 — Reference implementation and identity tests *(no build)*

Establish that the mathematics is right before anything is built on it.

- `tests/testthat/helper-soft-sankoff.R`: pure-R `softmin()`, soft-Sankoff
  down-pass and up-pass, an independent Felsenstein pruning reference, and a
  hard-Sankoff reference.
- `tests/testthat/test-ts-soft-sankoff.R` (Tier 2): the falsifiable claims —
  - `softmin_T -> min` as `T -> 0`, and `softmin_T <= min` always;
  - soft-Sankoff at `T -> 0` equals the hard reference and `ts_sankoff_test()`;
  - **soft-Sankoff at `T = 1` with `cost = -log P` equals Felsenstein pruning**,
    to machine precision, and equals `phangorn::pml()` where available;
  - equal-cost soft-Sankoff at `T -> 0` equals `TreeLength()`'s Fitch score;
  - total score is monotone non-increasing in `T`;
  - the up-pass at `T = 1` reproduces marginal ancestral-state probabilities.

  Exit criterion: all pass. If the `T = 1` identity fails, this whole document
  is wrong and the branch should be abandoned, not debugged.

**Status: DONE, 2026-08-01. 63/63 testthat assertions pass; 8/8 standalone
checks pass.** The load-bearing numbers:

| Claim | Result |
|---|---|
| `soft(T=1, cost = -log P)` vs own pruning reference | agrees to `1.8e-15` |
| `soft(T=1)` vs `phangorn::pml` | agrees to `0.0e+00` |
| `soft(T->0)` vs `ts_sankoff_test()` (compiled kernel) | exact |
| `soft(T->0)`, equal costs, vs `TreeLength()` Fitch | exact |
| `(min - softmin)/T` vs `log(#minima)` | agrees to `9e-10` |

The multiplicity law is worth stating as a result rather than a tolerance:
`softmin` does not converge to `min` at rate `T`, it converges to
`min - T*log(#minima)`, and that gap *is* the entropy of the near-optimal set.
It is the quantitative form of "temperature controls how much the criterion
integrates over reconstructions".

**Trap found, worth recording for anyone else writing an oracle against the
Sankoff bridge:** `ts_sankoff_test()` takes **0-based** tip states — the guard
is `state >= 0 && state < ns_ch` (`src/ts_rcpp.cpp`). The "1-based" comment at
`R/recode_hierarchy.R:188` means present states begin at index 1 of a 0-based
array, state 0 being "absent". An out-of-range index does not error; it leaves
every state at `INF`, so the whole score comes back `Inf`.

### Step 1 — Gate A *(no build)*

`dev/soft-sankoff/02-tilt-direction.R`, as specified above. Report the
correlation table into this document. **Blocks Steps 4 and 5.**

### Step 2 — Gate B *(needs a tarball build for the Fitch baseline)*

`dev/soft-sankoff/03-speed-budget.R`. Time soft-Sankoff (R, then a C++
prototype if R is too slow to be informative) against `TreeLength()` per tree
score, over a size ladder. **Blocks Step 4; informs Step 5.**

### Step 3 — Deliverables that depend on neither gate

These need only the scorer, and are worth doing even if both gates fail.

**3a. The parsimony/likelihood dial study — IMPLEMENTED, running.**
`dev/soft-sankoff/04-dial-study.R`. Search is held constant on purpose: the
candidate pool comes from ordinary hard-parsimony search (MPTs plus TBR
neighbours within 5 steps of optimal), and `T` only chooses among its members.
That is what dissolves the implementation confound. It also bounds the claim —
it measures which tree the criterion at `T` **prefers out of a common pool**,
not what a soft-objective search would find, and a criterion can only be
credited with recovering a tree the pool contains.

Reported per matrix per `T`: normalised `ClusteringInfoDist` to the generating
tree (averaged over trees the criterion cannot separate, since ties at `T = 0`
*are* the MPT set and are real), the hard parsimony score of the selection, a
consistency index as the homoplasy axis, and the Mk log-likelihood of the
selection as a secondary descriptive.

Sized at 3 matrices: ~4.5 s each with the Mk oracle, so ~7.5 min for 100 —
a Hamilton job. Submitted 2026-08-01 as job **18146187**
(`/nobackup/pjjg18/soft-sankoff/`), which installs the branch's TreeSearch into
a project-local library first, because `tsLib`'s TreeSearch 2.0.0 predates the
soft kernel.

#### RESULT — 100 matrices, two independent runs, 2026-08-01

Jobs **18146187** (run 1) and **18146295** (run 2, which adds a random-MPT null).
Run 2 is authoritative; run 1 is retained as
`dev/soft-sankoff/04-dial-study-run1.csv` because **the disagreement between the
two runs is itself part of the evidence.** Run 2 perturbs the RNG stream (one
extra `sample.int` for the null), which reshuffles pool construction — so where
the two runs disagree, the effect was never stable to begin with. 800/800 Mk fits
in both.

**Recovery by temperature** (median normalised `ClusteringInfoDist` to the
generating tree, lower better; run 2):

| `T` | median CID | median Mk logLik |
|---|---|---|
| 0 (hard parsimony) | 0.2481 | −649.38 |
| 0.02 / 0.05 / 0.10 | 0.2463 | −649.38 |
| 0.25 | 0.2474 | −649.38 |
| **0.50** | **0.2363** | **−648.81** |
| 1.00 | worse | −651.1 |
| 2.00 | worse | −651.9 |

**What is robust: `T = 0.5`.** Fixed-`T` per-matrix sign test against the MPT-set
mean, both runs:

| `T` | run 1 | run 2 | verdict |
|---|---|---|---|
| 0.02 | 56/33, p = 0.019 | 50/39, p = 0.289 | **not stable** |
| 0.10 | 56/33, p = 0.019 | 51/38, p = 0.203 | **not stable** |
| 0.25 | 56/36, p = 0.047 | 51/39, p = 0.246 | **not stable** |
| **0.50** | **65/33, p = 0.002** | **63/32, p = 0.002** | **STABLE** |
| 1.00 | 40/60, p = 0.057 | 41/59, p = 0.089 | consistently worse |

The low-`T` sign test moved from `p = 0.019` to `p = 0.289` under nothing but an
RNG perturbation. **Do not quote it.** Only `T = 0.5` survives: median CID gain
≈ 0.010, about 4% relative, with ~32 of 100 matrices still worse.

Also robust: `T >= 1` is worse on **both** axes, and `T = 2` badly so (20/80,
`p = 1e-9`). The dial has a genuine interior optimum. That was read here as
disposing of the Gate-A density worry in its own terms — a density tilt should
have kept helping as `T` rose. **That argument is too weak**; see the
generalisation test below, where the density signature reappears at 75 tips. An
interior optimum in mean CID is compatible with a criterion that chooses
decisively and is uncorrelated with truth.

**The mechanism: tie-breaking IS real, but the earlier evidence for it was
wrong.** Two corrections to what was first written here:

- `suboptimalSelections == 0` at low `T` is near-tautological. The MPT set is by
  definition where the optimal-scoring pool members are, so "the winner is an
  MPT" says nothing about how it was chosen.
- comparing one selected tree against the **mean** of the tied set is a weak test.
  The random-MPT null confirms it directly: a uniformly random MPT beats its own
  set's mean **47/42, `p = 0.67`, win rate 0.53**. So the ~0.56 win rates at low
  `T` were barely distinguishable from drawing at random.

The statistic that does settle it is the winner's **quantile rank among its own
MPT set's CIDs** — uniform at 0.5 under the null, and `NA` when the winner is not
itself an MPT:

| `T` | n | mean rank | `p` (Wilcoxon vs 0.5) |
|---|---|---|---|
| 0.02 | 92 | 0.338 | 2.8e−06 |
| 0.10 | 92 | **0.328** | **8.5e−07** |
| 0.25 | 86 | 0.337 | 4.8e−06 |
| 0.50 | 64 | 0.336 | 4.6e−05 |

**Soft-Sankoff genuinely ranks within the MPT set**: its pick sits around the
33rd percentile of the set's distance-to-truth distribution, not the 50th, at
`p ~ 1e-6` across 92 matrices. Ties at `T = 0` are a median of 11 trees and run to
104; only 8 of 100 matrices have a unique MPT, so this applies to 92% of them.

The two findings reconcile cleanly: **ranking is real but its CID payoff is small**
— small enough that a sign test against a mean cannot reliably see it, which is
exactly why the low-`T` sign test was unstable. `T = 0.5` wins more decisively
because it adds mild tolerance of suboptimality (33/100 selections leave the MPT
set, mean +0.50 steps) on top of the ranking.

| `T` | selections leaving the MPT set | mean extra steps |
|---|---|---|
| 0.02 / 0.05 / 0.10 | 0 / 100 | 0.00 |
| 0.25 | 7 / 100 | 0.08 |
| 0.50 | 33 / 100 | 0.50 |
| 1.00 | 84 / 100 | 2.39 |
| 2.00 | 100 / 100 | 4.03 |

**This is a live application that Gate B does not kill.** Ranking a retained tree
set is `O(pool)` rescores of ~10² trees, paid once — not the `O(candidates)` of a
search, which is what the x1147 penalty was priced against. A post-hoc
soft-Sankoff ranking of an MPT set is affordable today, and needs no annealing, no
incremental kernel, and no change to any search path.

**Homoplasy tracking: NOT supported, and now definitively so.**
Spearman(consistency index, best `T`) was `rho = -0.133` in run 1 and
`rho = +0.106` in run 2 — **the sign is not even stable**. The plan's "does the
answer track homoplasy?" question is answered no on this dataset.

**Caveats.** The pool is hard-parsimony-derived, so the criterion is only ever
credited with trees the pool contains. These are 22-tip simulated binary matrices;
nothing has been checked on empirical data or on `OReillyEtAl2016`. The per-matrix
*best*-`T` figure the script also prints (83/4/13, `p = 1.2e-13`) is a **ceiling,
not a result** — `bestT` is chosen using the answer. And the Gate-A CID/Mk
*disagreement* at `n = 1` does not reproduce: at `n = 100` the two measures agree
throughout.
#### GENERALISATION TEST — it does NOT transfer to 75 tips (jobs 18146497, 18146627)

The Step 3a result above is measured entirely on `congreveLamsdellMatrices`: 22
tips, 54 patterns, binary, low homoplasy. `dev/soft-sankoff/05-poolsize-pilot.R`
retests its headline statistic on O'Reilly et al. (2016) matrices — **75 tips**,
100 characters, 1000 available — before any full sweep was run.

**The rank statistic does not replicate.** Winner's quantile rank among its own
MPT set, `poolMaxSize = 100`, n = 39 per cell:

| `T` | mean rank | median | Wilcoxon vs 0.5 | chosen beats set mean |
|---|---|---|---|---|
| 0.02 | 0.442 | 0.515 | 0.210 | 19/39 |
| 0.10 | 0.445 | 0.535 | 0.215 | 19/39 |
| 0.25 | 0.469 | 0.535 | 0.424 | 18/39 |
| **0.50** | **0.475** | 0.535 | **0.597** | 18/39 |
| 1.00 | 0.386 | 0.273 | 0.118 | 25/39 |

Against the Congreve–Lamsdell reference of **0.328 at p ~ 1e-6 over 92
matrices**. **No cell is significant: the minimum `p` across all ten (cap, `T`)
combinations is 0.118.** `T = 0.5`, the cell where the C-L result was
*strongest*, lands at `p = 0.597` — indistinguishable from chance. The random-MPT
null beats its set mean 16/39, and the criterion manages 18–19/39.

**It is not a pool-truncation artifact.** `MaximizeParsimony()` returns *exactly*
`SearchControl()$poolMaxSize` on these matrices at both 100 and 300, so the MPT
sets genuinely exceed 300 and are search-truncated rather than exhausted. But
ranks are stable across the two caps (mean 0.445 vs 0.470, paired Wilcoxon
`p = 0.31`), so the statistic is measuring the criterion, not retention. The cap
was the obvious confound and it is ruled out.

**The failure has a shape, and it is the Gate-A shape.** The rank distribution is
**U-shaped** at every temperature (quintile counts ~11 / 4 / 6 / 5 / 10):

```
T=0.02  11  4  7  4 10
T=0.10  11  4  6  5 10
T=0.25  10  4  6  5 11
T=0.50   9  6  5  4 12
T=1.00  12  6  6  0 10
```

The criterion is **not** choosing at random — it chooses *decisively*, and is
near-best about as often as near-worst. That is a different failure from noise,
and it is precisely what a criterion tracking **reconstruction density** looks
like on data where density and distance-to-truth have come apart. In other words
**Gate A's original signed prediction was right, and it resurfaces here.** The
argument offered above — that a density tilt should have kept helping as `T` rose,
so the interior optimum disposes of the worry — was weaker than it read: an
interior optimum in mean CID is compatible with a criterion that is decisive and
uncorrelated with truth.

**Consequences.**

- The `T = 0.5` recovery advantage and the rank-0.33 result stand **on
  `congreveLamsdellMatrices` only**. They should not be described as properties
  of the criterion.
- The claim that post-hoc MPT ranking is "a live application Gate B does not
  kill" is **withdrawn as a general claim.** It holds where the ranking holds,
  which so far is one 22-tip low-homoplasy simulated dataset.
- The planned full O'Reilly sweep (~18 core-hours over 3 character counts) was
  **not run**, and should not be on this evidence. A pilot at n = 39 answered it.
- What would change the picture: identifying *what* the criterion is decisively
  tracking. The U-shape says there is signal there, just not truth-correlated
  signal. Testing rank against reconstruction ambiguity (MPR-set size per node,
  or the `T`-gap `(min − softmin)/T` = log #minima, which this kernel already
  computes) would say whether density is the thing — and that is a cheap test on
  data already staged.

**Original framing.** `CongreveLamsdell2016`,
`OReillyEtAl2016` and `Mk-prime-model` are all on disk. The standing
methodological weakness of that literature is that every parsimony-vs-likelihood
comparison confounds criterion with implementation — different programs,
different search intensity, different rooting conventions, different stopping
rules. A temperature parameter dissolves the confound: same DP, same search
kernel, same data, one scalar varied. Ask what `T` maximises recovery of the
generating tree, and whether the answer tracks homoplasy. This is a paper, and
the benchmark data are already here.

**3b. Graded ancestral-state reconstruction.** At `T > 0` the up-pass yields
marginal state probabilities per node per character at the same cost as the
down-pass. `PlotCharacter()` / `PaintCharacters.R` currently show MPR sets
all-or-nothing; a temperature-parameterised version shows graded reconstructions
with an explicit statement of how much integration is being done. Does not touch
search.

### Step 4 — Annealed search *(DEAD: Gate B failed, 2026-08-01)*

**Do not build this.** A compiled soft score costs a median x1147 of a Fitch
score for binary characters against a x50 orientation threshold, plus a bounded
x386 incrementality forfeit on top. The propose/arbitrate architecture below
means scoring twice, so it pays that cost on every candidate. Reviving it needs
roughly three orders of magnitude from a mitigation, not one. The design is kept
below for the record, and because Step 5 references it.

Temper the **objective**, not the acceptance rule. `ts_drift.cpp` and
`ts_temper.cpp` already do the latter; this is a different mechanism (graduated
non-convexity): climb at high `T` where basins are wide, anneal to `T -> 0`.

Architecture keeps exactness: the soft score **proposes**, hard SIMD Fitch
**arbitrates**, and the reported score is always the hard one. The cost is
scoring twice, which is what Gate B prices.

Numerical requirements: shifted log-sum-exp throughout; as `T -> 0` the softmax
concentrates and gradients vanish, the classic continuation-method failure, so
the annealing schedule needs its own small study.

### Step 5 — Gradients, and only then the embedding *(DEAD: inherits Gate B)*

With per-branch cost scalars — which soft-Sankoff has natively, since at `T = 1`
they *are* branch lengths — `d(score)/d(b)` exists and Vine's four-Jacobian chain
closes: embedding -> distances -> NJ -> `(tau, b)` -> soft score -> back. This is
the honest route to embedding-space search for a parsimony-family criterion.

It is last for a reason: it is downstream of both gates, and it inherits the two
open questions recorded against the Vine discussion — reachability of every
topology from a `d`-dimensional embedding, and the "differentiate conditional on
the discrete choice" caveat, which here stacks twice (NJ join order, and the
softmax's dominant term as `T -> 0`).

**Independent of all of the above:** the trustworthiness/continuity measurement
on embedding perturbations (Smith 2022, `syab099`/`syab100`) should be run before
any embedding search code is written, and does not depend on soft-Sankoff at all.
If a small embedding perturbation lands ~15 splits away with no intermediates,
the move is a random restart in costume.

---

## Data and compute

**All of Gate A and Step 3a are Hamilton jobs, not local runs.** Anything over
~45 s of compute goes through the `/hamilton` skill. The 2026-08-01 first read
used 6 of the 100 available Congreve & Lamsdell matrices *because* it was sized
to run locally, and produced exactly one usable data point. Do not repeat that:
size the run to the question and submit it.

| Input | State |
|---|---|
| `congreveLamsdellMatrices` | **100 available**, bundled in `data/`. Only 6 used so far. 22 tips, 55 sites. |
| `referenceTree` | bundled; the generating tree for the above |
| O'Reilly 2016 result trees | on disk, `OReillyEtAl2016/data-raw/Trees/` — **NOT usable as a candidate pool.** Each `.sym` is ONE support-tagged tree per matrix, and `SupportSuboptimal()` in that repo's `GenerateData.Rmd` builds its "suboptimal" set by progressively *collapsing low-support nodes* (`ReduceTreesBySupport` -> `CollapseNode`) — a resolution series, not a score-based near-optimal set. The trees are non-binary, so both `TreeLength()` and the soft kernel reject them. Pools must be built from the matrices instead. |
| O'Reilly 2016 **matrices** | **Fetched by the maintainer 2026-08-01** to `~/Downloads/doi_10_5061_dryad_10qf3__v20160322.zip` (19,061,436 bytes). **Nested archive** — the Dryad wrapper contains a single member, `oreilly2016matrices.zip` (19,061,276 bytes); unzip twice. Not yet unpacked or staged. Suggested home: `OReillyEtAl2016/data-raw/Matrices/`, which that repo already gitignores. |

**Hamilton libraries** are documented in the `/hamilton` skill's
`r-infrastructure.md`: `tsLib = /nobackup/pjjg18/TreeSearch/lib` already carries
TreeSearch 2.0.0, TreeTools, TreeDist and Quartet, and must precede `baseLib`.

**Correction, 2026-08-01: `phangorn` IS present in `tsLib`**
(`/nobackup/pjjg18/TreeSearch/lib/phangorn`), so the Mk oracle needs no
project-local install. The earlier claim that it was absent came from reading
`r-infrastructure.md`'s summary table rather than listing the directory; that
table is a partial listing, not an inventory.

What *does* need a project-local install is **TreeSearch itself**. `tsLib`'s
TreeSearch 2.0.0 predates the soft kernel, so a job using
`ts_soft_sankoff_test()` must install this branch's build into a library ahead of
`tsLib` on `R_LIBS`. The job script at `/nobackup/pjjg18/soft-sankoff/` does
that, and then smoke-tests that the symbol is actually registered before running
anything — a missing entry in `TreeSearch-init.c` makes the binding fail at
*call* time, not at load time, so a job would otherwise get most of the way in
before dying.

## Costs and known limitations

- **Speed** — Gate B. The dominant risk.
- **Surrogate optimum** — the `T > 0` optimum is not the `T = 0` optimum.
  Mitigated by the propose/arbitrate split, at the cost of double scoring.
- **Inapplicables** — Brazeau–Guillerme–Smith is not a Sankoff DP and has no
  obvious soft form. But Goloboff's x-transformation recoding of hierarchical
  characters *is* Sankoff and is already the alternative pathway in-tree
  (`R/recode_hierarchy.R`), so hierarchical data is not excluded; it routes
  differently. Note the live rooting defects on that pathway (T-374, T-385) —
  a soft variant inherits them and must not be built on top of an unsettled
  objective.
- **Topology stays discrete.** Softening `min` does not soften the tree. Steps
  0–4 buy a smoothed objective for the same discrete search; only Step 5
  attempts continuity, and only via the embedding.
- **Interpretation** — at `T > 0` this is no longer parsimony, and reviewers of
  a parsimony package will say so. Step 3a is the answer to that objection
  (the dial is the point), not a defence against it.

## What is deliberately not proposed

- No change to the default scoring path. Nothing here alters what
  `MaximizeParsimony()` does today.
- No `to-do.md` claim or coordination commit has been made; per `AGENTS.md`
  those belong on `cpp-search`, and this exploration has not been assigned a
  task number. Suggest **T-386** if it is taken up.
- No soft-min relaxation of the HSJ DP. HSJ requires genuine rooting-invariance
  (T-374); softening it would convert an unsettled objective into a stably
  unsettled one.
