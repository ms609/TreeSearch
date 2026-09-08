# Test File Conventions

Load this when: adding or modifying `tests/testthat/test-ts-*.R` files,
choosing test tiers, or writing test helpers.

---

## Conventions

All `tests/testthat/test-ts-*.R` files must use `TreeSearch:::` to call
internal C++ bridge functions. Define short local wrappers for readability.

Shared helpers are in `tests/testthat/helper-ts.R` (`make_ts_data()`,
`ts_score()`, `validate_result()`, `skip_extended()`).

**Never use `%in%` on Splits objects in test files** — S3 dispatch fails
in the cloned namespace created by `test_check()`. Use `as.logical()`
matrix comparison instead.

---

## Test tiering

Every new `test-ts-*.R` file must be assigned to one of three tiers.
See `tests/testing-strategy.md` for the full rationale.

| Tier | Guard | When it runs | Use for |
|------|-------|-------------|---------|
| 1 — CRAN | none | always (CRAN + CI + local) | Fast (< ~2 s) API and data-structure unit tests |
| 2 — CI | `skip_on_cran()` at **file level** (first executable line) | CI + local | C++ engine correctness, scoring, search algorithms |
| 3 — Extended | `skip_extended()` at **file level** | `TREESEARCH_EXTENDED_TESTS=true` only | Stress tests, benchmarks, timing measurements |

**Default for new `test-ts-*` files: Tier 2.** Add `skip_on_cran()` as the
very first executable line (before any helpers or `test_that()` calls):

```r
# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
```

Use Tier 3 only for tests that take > ~10 s or are sensitive to machine load.
