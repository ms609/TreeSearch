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
- **Step 3a: running.** Hamilton job **18146187**,
  `/nobackup/pjjg18/soft-sankoff/`. A 3-matrix pilot showed an interior optimum
  near `T = 0.25`–`0.5` beating hard parsimony on distance to truth (0.235 vs
  0.262), with the Mk arm agreeing rather than disagreeing. `n = 3`; the full
  run is what counts.

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
