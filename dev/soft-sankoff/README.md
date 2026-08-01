# Soft-Sankoff exploration

Evidence scripts for
[`dev/plans/2026-08-01-soft-sankoff-temperature-dial.md`](../plans/2026-08-01-soft-sankoff-temperature-dial.md).

Sankoff's DP with the `min` softened at temperature `T` is hard weighted
parsimony at `T -> 0` and Felsenstein's pruning algorithm at `T = 1` with
`cost = -log(P)`. Parsimony and likelihood are the same dynamic program in two
semirings; `T` is the dial.

| Script | Purpose | Needs a build? |
|---|---|---|
| `01-identity-check.R` | Step 0 — the semiring identity, numerically | no |
| `02-tilt-direction.R` | **Gate A — RETIRED.** Kept as the record of its first read. Do not widen: it existed to protect Step 4, and Gate B killed Step 4 | yes |
| `03-speed-budget.R` | **Gate B** — cost per score against SIMD Fitch | yes |
| `04-dial-study.R` | **Step 3a** — what `T` best recovers the generating tree, and does it track homoplasy? | yes |
| `build.sh` | install into `.agent-softsankoff` | — |

Reference implementation (pure R): `tests/testthat/helper-soft-sankoff.R`.
Compiled prototype: `src/ts_soft_sankoff.{h,cpp}` plus the
`ts_soft_sankoff_test()` bridge in `src/ts_soft_sankoff_bridge.cpp` — new files
only, reachable from no default scoring path, and deliberately not sharing a
struct with `ts_sankoff.{h,cpp}`, which sits on the live x-transformation pathway
with its open T-374/T-385 rooting defects.

Regression tests, both Tier 2:
`tests/testthat/test-ts-soft-sankoff.R` (R reference, 63 assertions) and
`tests/testthat/test-ts-soft-sankoff-cpp.R` (compiled kernel, 411 assertions,
routed to `FelsensteinLogLik()`, `phangorn::pml()` and the pre-existing
`ts_sankoff_test()` so that a shared conceptual error cannot hide in a passing
test).

## Running

```bash
# Step 0 needs nothing but ape (+ phangorn for the third-party oracle)
Rscript dev/soft-sankoff/01-identity-check.R

# Everything else needs TreeSearch installed to the isolated library.
# CCACHE_DISABLE=1 and --preclean matter: this branch adds a new header, and
# ccache plus a header change is the recipe for a stale .o against an old
# struct layout.
sh dev/soft-sankoff/build.sh

Rscript dev/soft-sankoff/03-speed-budget.R
Rscript dev/soft-sankoff/04-dial-study.R          # all 100 matrices, ~7.5 min
Rscript dev/soft-sankoff/04-dial-study.R 3 200 no # quick local smoke

# The testthat files are Tier 2, so NOT_CRAN must be set or nothing runs.
# Without it testthat reports [ FAIL 0 | PASS 0 ], which reads as success.
NOT_CRAN=true Rscript -e "library(TreeSearch, lib.loc='.agent-softsankoff'); \
  setwd('tests/testthat'); source('helper-soft-sankoff.R'); \
  testthat::test_file('test-ts-soft-sankoff-cpp.R')"
```

`04-dial-study.R` honours `SOFT_SANKOFF_LIB` and `SOFT_SANKOFF_OUT`, so the same
script runs locally and on Hamilton with no cluster-specific fork.

## Status, 2026-08-01

- **Step 0: done.** 63/63 R-reference assertions, 8/8 standalone checks. `T = 1`
  matches `phangorn::pml` exactly; `T -> 0` matches the compiled Sankoff kernel
  and `TreeLength()` exactly.
- **Compiled prototype: done.** 411/411 assertions against independent oracles.
  Bought a median **x100** speedup over the pure-R reference, which is the
  figure the plan had asserted without measuring.
- **Gate B: CLOSED, FAILED.** Compiled soft against an in-C++ Fitch full
  rescore is a median **x1147** for binary characters (x2981 for `k = 4`),
  against a x50 orientation threshold — all 8 cells over, by 20–70x. Roughly 3x
  worse than the x350 op-count estimate, because `exp` is not one operation.
  The *complexity* term (loss of incremental rescoring) stays open but is
  bounded: Fitch's own full-rescore-to-incremental-candidate ratio is a median
  **x386**, and that is what a full-rescore soft kernel forfeits on top. Moot
  anyway — the constant alone fails.
  **Consequence: Steps 4 and 5 are dead on cost.**
- **Gate A: RETIRED, not answered.** It protected Step 4. See the plan doc for
  why two of its first-read observations change meaning rather than standing.
- **Step 3a: DONE**, 100 matrices, **two independent runs** (Hamilton `18146187`
  then `18146295`, which adds a random-MPT null). Results in `04-dial-study.csv`
  (run 2, authoritative), `04-dial-study-run1.csv`, `04-dial-study-per-matrix.csv`.
  Run 2 perturbs the RNG stream, so **where the runs disagree the effect was never
  stable** — that disagreement is part of the evidence.
  - **Robust: `T = 0.5` beats hard parsimony.** Median normalised CID to the
    generating tree 0.2363 vs 0.2481. Fixed `T = 0.5`, per matrix: 65/33
    (p = 0.002) in run 1 and 63/32 (p = 0.002) in run 2. Modest — ~4% relative,
    ~32 of 100 matrices still worse.
  - **NOT robust: the low-`T` sign test.** `T <= 0.25` swung from p = 0.019 to
    p = 0.289 on nothing but an RNG perturbation. Do not quote it.
  - **`T >= 1` is worse on both axes**, `T = 2` badly (20/80, p = 1e-9). The dial
    has a genuine interior optimum — which also disposes of the Gate-A density
    worry in its own terms, since a density tilt should have kept helping as `T`
    rose.
  - **Tie-breaking is real, but the first evidence for it was wrong.** The
    random-MPT null shows a random MPT beats its own set's mean 47/42 (p = 0.67,
    rate 0.53), so the ~0.56 low-`T` win rates were near-indistinguishable from
    chance; and `suboptimalSelections == 0` is near-tautological. What settles it
    is the winner's **quantile rank among its own MPT set's CIDs**: mean 0.328–0.338
    against a null of 0.5, **p ~ 1e-6 across 92 matrices**. Soft-Sankoff's pick
    sits near the 33rd percentile of the set's distance-to-truth distribution.
    Ranking is real; its CID payoff is just small enough that a sign test against
    a mean cannot reliably see it.
  - **Gate B does not kill this application.** Ranking a retained set is
    `O(pool)` rescores of ~10² trees paid once, not the `O(candidates)` the x1147
    penalty was priced against. No annealing, no incremental kernel, no search
    change.
  - **Homoplasy tracking: definitively not supported.** rho = −0.133 in run 1,
    **+0.106** in run 2 — the sign is not even stable.
  - The per-matrix *best*-`T` figure the script prints (83/4/13, p = 1.2e-13) is a
    **ceiling, not a result**: `bestT` is chosen using the answer.
  - The Gate-A CID/Mk **disagreement at `n = 1` does not reproduce** at n = 100.

- **GENERALISATION: the Step 3a findings do NOT transfer.** Pilot
  `05-poolsize-pilot.R` (jobs `18146497`, `18146627`) retested the rank statistic
  on O'Reilly 2016 matrices — 75 tips, 100 characters, n = 39.
  - Chance-level at **every** temperature: mean rank 0.442–0.475, and the
    minimum `p` across all ten (cap, `T`) cells is **0.118**. `T = 0.5`, the
    strongest Congreve–Lamsdell cell, gives `p = 0.597`. Reference: 0.328 at
    p ~ 1e-6 on 22 tips.
  - **Not** a pool-truncation artifact. `MaximizeParsimony()` returns exactly
    `poolMaxSize` at both 100 and 300 (so the sets are truncated), yet ranks are
    stable across the two caps (paired `p = 0.31`).
  - The distribution is **U-shaped** (~11/4/6/5/10) at every `T`: the criterion
    chooses *decisively* and is near-best about as often as near-worst — the
    signature of a density-tracking criterion where density and truth diverge,
    i.e. **Gate A's original prediction resurfacing**.
  - So the MPT-ranking application is **withdrawn as a general claim**; it holds
    on one 22-tip low-homoplasy simulated dataset. The planned ~18 core-hour
    O'Reilly sweep was **not run** and should not be.
  - Cheap next test if anyone picks this up: regress rank against reconstruction
    ambiguity — `(min − softmin)/T` = log #minima, which this kernel already
    computes — to find out what the criterion *is* decisively tracking.

### Gotchas worth not rediscovering

- **A new `// [[Rcpp::export]]` is not enough.** `Rcpp::compileAttributes()`
  writes the wrapper, but `src/TreeSearch-init.c` hand-maintains the
  `R_CallMethodDef` table and needs its own entry. Without it the binding fails
  at *call* time with `object '_TreeSearch_...' not found`, not at load time —
  so every test in a file errors while the file still reports as having run.
- **`phangorn` IS in `tsLib`** on Hamilton
  (`/nobackup/pjjg18/TreeSearch/lib/phangorn`), contrary to the earlier note in
  the plan. What needs a project-local install is TreeSearch itself: `tsLib`'s
  2.0.0 predates the soft kernel.
- **`MaximizeParsimony()` returns POLYTOMIES by default.** `collapse = TRUE`
  contracts zero-length (unsupported) branches before returning, and both
  `TreeLength()` and the soft kernel refuse a non-binary tree. Pass
  `collapse = FALSE`. The first Hamilton submission died at matrix 11 of 100 on
  exactly this — the sweep now also resolves defensively, so a future default
  change cannot fail it late again.
- **`ts_bench_tbr_phases()` clocks whole microseconds.** At 100 patterns a Fitch
  rescore is 1–7 ticks and one run reported `0 us`. Compare at a pattern count
  high enough to resolve the denominator; both kernels are linear in it.
- **Parsimony trees have no branch lengths** and `phangorn::pml()` refuses them
  outright. Seed them before `optim.pml`, and do not wrap the call in a bare
  `try()` — that turned a real failure into an all-NA column on a run that
  otherwise looked fine.
