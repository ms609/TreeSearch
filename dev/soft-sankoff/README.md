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
| `02-tilt-direction.R` | **Gate A** — does the soft score point toward the truth, or toward MPT density? | yes |
| `03-speed-budget.R` | **Gate B** — cost per score against SIMD Fitch | yes |

Reference implementation (pure R): `tests/testthat/helper-soft-sankoff.R`.
Regression tests: `tests/testthat/test-ts-soft-sankoff.R` (Tier 2).

## Running

```bash
# Step 0 needs nothing but ape (+ phangorn for the third-party oracle)
Rscript dev/soft-sankoff/01-identity-check.R

# Gates A and B need TreeSearch installed to an isolated library
SRC=$(pwd) && TMPBUILD=$(mktemp -d) && rm -f src/*.o src/*.dll && \
  (cd "$TMPBUILD" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$SRC") && \
  R CMD INSTALL --library=.agent-softsankoff "$TMPBUILD"/TreeSearch_*.tar.gz

Rscript dev/soft-sankoff/02-tilt-direction.R 20 48
Rscript dev/soft-sankoff/03-speed-budget.R

# The testthat file is Tier 2, so NOT_CRAN must be set or nothing runs
NOT_CRAN=true Rscript -e "library(TreeSearch, lib.loc='.agent-softsankoff'); \
  testthat::test_file('tests/testthat/test-ts-soft-sankoff.R')"
```

## Status, 2026-08-01

- **Step 0: done.** 63/63 test assertions, 8/8 standalone checks. `T = 1`
  matches `phangorn::pml` exactly; `T -> 0` matches the compiled Sankoff kernel
  and `TreeLength()` exactly.
- **Gate A: insufficient data.** Only 1 of 6 Congreve & Lamsdell matrices had
  >= 8 distinct MPTs. Design needs widening to near-optimal trees.
- **Gate B: unfavourable.** Operation-count ratio x350 against a x50
  orientation threshold, before counting `exp` cost or the incremental-rescore
  loss.
