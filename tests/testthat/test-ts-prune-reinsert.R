# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for taxon pruning-reinsertion perturbation (T-266).
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search with prune-reinsert enabled
ts_driven_pri <- function(ds, pruneReinsertCycles = 2L,
                          pruneReinsertDrop = 0.15,
                          pruneReinsertSelection = 0L,
                          maxReplicates = 3L, targetHits = 2L,
                          ratchetCycles = 2L, driftCycles = 0L,
                          xssRounds = 0L, fuseInterval = 5L,
                          verbosity = 0L, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    driftCycles = driftCycles,
    xssRounds = xssRounds,
    fuseInterval = fuseInterval,
    pruneReinsertCycles = pruneReinsertCycles,
    pruneReinsertDrop = pruneReinsertDrop,
    pruneReinsertSelection = pruneReinsertSelection,
    verbosity = verbosity,
    ...
  )
}

# ---------- Test datasets ----------

small_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
small_dataset <- MatrixToPhyDat(small_mat)
small_ds <- make_ts_data(small_dataset)

set.seed(6817)
med_mat <- matrix(sample(0:1, 20 * 15, replace = TRUE),
                  nrow = 20,
                  dimnames = list(paste0("t", 1:20), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)


# ---------- Tests ----------

test_that("Prune-reinsert produces valid search result", {
  set.seed(3491)
  result <- ts_driven_pri(small_ds, pruneReinsertCycles = 2L)

  expect_true(is.list(result))
  expect_true("trees" %in% names(result))
  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
  expect_true(result$best_score > 0)
  expect_true(result$best_score < Inf)
})

test_that("Prune-reinsert does not degrade search quality", {
  set.seed(8204)
  with_pri <- ts_driven_pri(med_ds, pruneReinsertCycles = 3L,
                            maxReplicates = 5L, ratchetCycles = 3L)
  set.seed(8204)
  without_pri <- ts_driven_pri(med_ds, pruneReinsertCycles = 0L,
                               maxReplicates = 5L, ratchetCycles = 3L)

  # Allow small tolerance — on small datasets the effect may be zero
  expect_true(with_pri$best_score <= without_pri$best_score + 5L)
})

test_that("Prune-reinsert with 0 cycles is effectively disabled", {
  set.seed(5723)
  result <- ts_driven_pri(small_ds, pruneReinsertCycles = 0L)

  expect_true(is.list(result))
  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
  # Timing should be zero when disabled
  expect_equal(result$timings[["prune_reinsert_ms"]], 0)
})

test_that("Prune-reinsert timing is reported", {
  set.seed(7340)
  result <- ts_driven_pri(small_ds, pruneReinsertCycles = 2L)

  expect_true("prune_reinsert_ms" %in% names(result$timings))
  expect_true(result$timings[["prune_reinsert_ms"]] >= 0)
})

test_that("Prune-reinsert with instability selection runs without error", {
  set.seed(9127)
  result <- ts_driven_pri(med_ds,
                          pruneReinsertCycles = 2L,
                          pruneReinsertSelection = 1L,
                          maxReplicates = 3L)

  expect_true(is.list(result))
  expect_true(length(result$trees) > 0)
  validate_result(result, 20L)
})

test_that("Prune-reinsert with high drop fraction runs safely", {
  # Drop 40% of 10 tips = 4 tips, leaving 6 (above minimum of 4)
  set.seed(2158)
  result <- ts_driven_pri(small_ds,
                          pruneReinsertCycles = 1L,
                          pruneReinsertDrop = 0.40)

  expect_true(is.list(result))
  validate_result(result, 10L)
})

test_that("Prune-reinsert with very high drop fraction does not crash", {
  # Drop 90% of 10 tips = 9, but capped at n-4 = 6
  set.seed(4610)
  result <- ts_driven_pri(small_ds,
                          pruneReinsertCycles = 1L,
                          pruneReinsertDrop = 0.90)

  expect_true(is.list(result))
  validate_result(result, 10L)
})

test_that("Prune-reinsert works with medium-sized dataset", {
  set.seed(1856)
  result <- ts_driven_pri(med_ds,
                          pruneReinsertCycles = 3L,
                          pruneReinsertDrop = 0.10,
                          maxReplicates = 3L)

  expect_true(is.list(result))
  validate_result(result, 20L)
  expect_true(result$best_score > 0)
})

test_that("Prune-reinsert interacts correctly with ratchet", {
  set.seed(6089)
  result <- ts_driven_pri(med_ds,
                          pruneReinsertCycles = 2L,
                          ratchetCycles = 4L,
                          maxReplicates = 3L)

  expect_true(is.list(result))
  validate_result(result, 20L)
  expect_true(result$best_score > 0)
  expect_true(result$timings[["ratchet_ms"]] >= 0)
  expect_true(result$timings[["prune_reinsert_ms"]] >= 0)
})

test_that("Prune-reinsert with outer cycles divides evenly", {
  set.seed(3795)
  result <- ts_driven_pri(med_ds,
                          pruneReinsertCycles = 4L,
                          maxReplicates = 2L,
                          outerCycles = 2L)

  expect_true(is.list(result))
  validate_result(result, 20L)
})

test_that("SearchControl accepts prune-reinsert parameters", {
  ctrl <- TreeSearch::SearchControl(
    pruneReinsertCycles = 3L,
    pruneReinsertDrop = 0.15,
    pruneReinsertSelection = 1L
  )
  expect_equal(ctrl$pruneReinsertCycles, 3L)
  expect_equal(ctrl$pruneReinsertDrop, 0.15)
  expect_equal(ctrl$pruneReinsertSelection, 1L)
})

test_that("SearchControl prune-reinsert defaults are correct", {
  ctrl <- TreeSearch::SearchControl()
  expect_equal(ctrl$pruneReinsertCycles, 0L)
  expect_equal(ctrl$pruneReinsertDrop, 0.10)
  expect_equal(ctrl$pruneReinsertSelection, 0L)
})
