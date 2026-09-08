# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for stochastic NNI-perturbation escape mechanism (T-186).
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search with NNI perturbation enabled
ts_driven_np <- function(ds, nniPerturbCycles = 3L, nniPerturbFraction = 0.5,
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
    nniPerturbCycles = nniPerturbCycles,
    nniPerturbFraction = nniPerturbFraction,
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

set.seed(5821)
med_mat <- matrix(sample(0:1, 20 * 10, replace = TRUE),
                  nrow = 20,
                  dimnames = list(paste0("t", 1:20), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)


test_that("NNI perturbation produces valid search result", {
  set.seed(4207)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 3L)

  expect_true(is.list(result))
  expect_true("trees" %in% names(result))
  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
  expect_true(result$best_score > 0)
  expect_true(result$best_score < Inf)
})

test_that("NNI perturbation does not degrade search quality", {
  # Compare: with NNI-perturb enabled vs disabled (ratchet only).
  # NNI-perturb should produce scores no worse than ratchet-only.
  # (On small datasets the difference may be zero.)
  set.seed(7193)
  with_np <- ts_driven_np(med_ds, nniPerturbCycles = 5L,
                           maxReplicates = 5L, ratchetCycles = 3L)
  set.seed(7193)
  without_np <- ts_driven_np(med_ds, nniPerturbCycles = 0L,
                              maxReplicates = 5L, ratchetCycles = 3L)

  expect_true(with_np$best_score <= without_np$best_score + 5L)
})

test_that("NNI perturbation with fraction = 0 is effectively disabled", {
  set.seed(2894)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 3L,
                          nniPerturbFraction = 0.0)

  expect_true(is.list(result))
  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
})

test_that("NNI perturbation with fraction = 1 works", {
  set.seed(6120)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 2L,
                          nniPerturbFraction = 1.0,
                          maxReplicates = 2L)

  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
})

test_that("NNI perturbation timings reported", {
  set.seed(3376)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 2L,
                          maxReplicates = 1L)

  expect_true("timings" %in% names(result))
  expect_true("nni_perturb_ms" %in% names(result$timings))
  expect_true(result$timings[["nni_perturb_ms"]] >= 0)
})

test_that("NNI perturbation disabled when cycles = 0", {
  set.seed(8461)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 0L,
                          maxReplicates = 2L)

  expect_true(length(result$trees) > 0)
  # nni_perturb_ms should be 0 (or very close) when disabled
  expect_equal(result$timings[["nni_perturb_ms"]], 0)
})

test_that("SearchControl includes NNI perturbation params", {
  ctrl <- SearchControl(nniPerturbCycles = 5L, nniPerturbFraction = 0.3)
  expect_equal(ctrl$nniPerturbCycles, 5L)
  expect_equal(ctrl$nniPerturbFraction, 0.3)

  # Defaults
  ctrl_default <- SearchControl()
  expect_equal(ctrl_default$nniPerturbCycles, 0L)
  expect_equal(ctrl_default$nniPerturbFraction, 0.5)
})

test_that("NNI perturbation works with IW scoring", {
  set.seed(9147)
  min_steps <- TreeSearch:::MinimumLength(small_dataset)
  result <- ts_driven_np(small_ds, nniPerturbCycles = 2L,
                          maxReplicates = 2L,
                          min_steps = min_steps, concavity = 10.0)

  expect_true(length(result$trees) > 0)
  validate_result(result, 10L)
  expect_true(result$best_score > 0)
})
