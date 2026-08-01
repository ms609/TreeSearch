# Soft-Sankoff: a temperature dial between parsimony and likelihood

**Date:** 2026-08-01
**Status:** exploration. No `src/` change proposed yet; two gates must pass first.
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

**Provisional consequence:** Step 4 (annealed search) is in serious doubt, and
Step 5 inherits that doubt. Steps 3a and 3b do not, which is why they were
placed off the gate. Anyone reviving Step 4 needs a mitigation — SIMD over the
character axis, a `k = 2` special case, or restricting soft scoring to a coarse
outer loop — not just a better `exp`.

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

**3a. The parsimony/likelihood dial study.** `CongreveLamsdell2016`,
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

### Step 4 — Annealed search *(gated on A and B)*

Temper the **objective**, not the acceptance rule. `ts_drift.cpp` and
`ts_temper.cpp` already do the latter; this is a different mechanism (graduated
non-convexity): climb at high `T` where basins are wide, anneal to `T -> 0`.

Architecture keeps exactness: the soft score **proposes**, hard SIMD Fitch
**arbitrates**, and the reported score is always the hard one. The cost is
scoring twice, which is what Gate B prices.

Numerical requirements: shifted log-sum-exp throughout; as `T -> 0` the softmax
concentrates and gradients vanish, the classic continuation-method failure, so
the annealing schedule needs its own small study.

### Step 5 — Gradients, and only then the embedding *(gated on B)*

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
