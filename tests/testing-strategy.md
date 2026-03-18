# TreeSearch Test Tiering Strategy

Three tiers of tests with distinct run conditions.

---

## Tier 1 — CRAN unit tests (always run)

**Purpose:** Catch breaking changes caused by edits to dependency packages
(ape, TreeTools, TreeDist, etc.). Run on every `R CMD check` including CRAN.

**Criterion:** Fast (< ~2 s per file), test the R-level public API or basic
data-structure invariants, no `skip_on_cran()`.

**Files (ts- engine):**
| File | What it covers |
|------|----------------|
| `test-ts-constraint-small.R` | Constraint logic, small dataset |
| `test-ts-memory-layout.R` | Internal memory layout invariants |
| `test-ts-pool.R` | Tree-pool deduplication |
| `test-ts-simd.R` | SIMD bit-parallel scoring correctness |
| `test-ts-splits.R` | Split hashing and comparison |
| `test-ts-start-tree.R` | Starting-tree API |

**Files (R-level API):** All `test-*.R` files that do NOT carry a `ts-` prefix
(e.g. `test-tree_length.R`, `test-AdditionTree.R`, `test-Morphy.R`, ...). These
verify the public R interface against dependency changes.

---

## Tier 2 — CI coverage tests (`skip_on_cran()`)

**Purpose:** Guarantee code coverage of the C++ engine internals. Run on every
CI platform (`NOT_CRAN=true`) but not on CRAN.

**Guard:** `skip_on_cran()` — either file-level (first executable line in the
test file, before any `test_that()`) or per-test.

**Files:** All remaining `test-ts-*.R` files not in Tier 1 or Tier 3:

```
test-ts-char-ordering.R     test-ts-ratchet-opt.R
test-ts-css.R               test-ts-ratchet-search.R
test-ts-drift-search.R      test-ts-resample.R
test-ts-driven.R            test-ts-sector.R
test-ts-fuse.R              test-ts-simplify.R
test-ts-iw.R                test-ts-spr-nni-opt.R
test-ts-na-incremental.R    test-ts-tabu.R
test-ts-parallel.R          test-ts-tbr-search.R
test-ts-profile.R           test-ts-tbr-symmetry.R
test-ts-progress.R          test-ts-wagner.R
```

---

## Tier 3 — Extended algorithmic / stress tests

**Purpose:** Verify algorithmic correctness and catch performance regressions
under large or adversarial inputs. Run locally during development and in a
dedicated periodic CI workflow, but NOT on every push/PR build.

**Guard:** `skip_extended()` (defined in `tests/testthat/helper-ts.R`).
Enabled by setting `TREESEARCH_EXTENDED_TESTS=true` in the environment.

**Files:**

| File | Nature |
|------|--------|
| `test-ts-timings.R` | Timing measurements (fragile on shared CI runners) |
| `test-ts-tbr-bench.R` | TBR optimization benchmark / correctness |
| `test-ts-ratchet-stress.R` | Ratchet stress test across many datasets |
| `test-ts-resample-stress.R` | Resample + SA stress test |

**Running extended tests locally:**

```bash
TREESEARCH_EXTENDED_TESTS=true Rscript -e \
  "testthat::test_dir('tests/testthat', filter='ts-')"
```

Or to run the full suite with extended tests enabled:
```bash
TREESEARCH_EXTENDED_TESTS=true R CMD check --no-build-vignettes .
```

---

## GHA workflows

| Workflow | `NOT_CRAN` | `TREESEARCH_EXTENDED_TESTS` | Tiers run |
|----------|-----------|----------------------------|-----------|
| `R-CMD-check.yml` (push/PR, 6 platforms) | `true` | unset | 1 + 2 |
| `extended-tests.yml` (scheduled weekly) | `true` | `true` | 1 + 2 + 3 |

---

## Adding a new test — checklist

1. **Tier 1:** No guard. Must complete in < ~2 s. Tests the R-level API or
   a data-structure invariant.
2. **Tier 2:** Add `skip_on_cran()` as the first line inside each `test_that()`
   (or once at file level). Tests internal C++ correctness.
3. **Tier 3:** Add `skip_extended()` as the first line inside each `test_that()`
   (or once at file level). Suitable for stress tests, benchmarks, or anything
   that takes > ~10 s.
