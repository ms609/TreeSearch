# Tier 1: no skip guard — fast API tests only; see tests/testing-strategy.md
# Tests for the maxReplicates adequacy warning in MaximizeParsimony().
#
# Formula: min_reps = max(10, ceiling(nTip * sum(weight) / 5000))
# Warning fires only when user explicitly passes maxReplicates < min_reps,
# and only for datasets with nTip >= 30.

library("TreeTools", quietly = TRUE)

data("inapplicable.phyData", package = "TreeSearch")

# Helper: make a minimal 30-taxon binary phyDat (floor test datasets).
# Use a fixed seed so sum(weight) is deterministic.
.make_30tip_ds <- function() {
  set.seed(8812)
  mat <- matrix(
    sample(c("0", "1"), 30L * 30L, replace = TRUE),
    nrow = 30L, ncol = 30L,
    dimnames = list(paste0("t", seq_len(30L)), NULL)
  )
  MatrixToPhyDat(mat)
}

# ---- Floor-driven cases (nTip=30, formula gives max(10, ...) = 10) ----------

test_that("explicit maxReplicates below floor (10) triggers warning", {
  ds30 <- .make_30tip_ds()
  expect_warning(
    MaximizeParsimony(ds30, maxReplicates = 3L, targetHits = 1L,
                      maxSeconds = 0.5, verbosity = 1L),
    regexp = "replicates are recommended",
    ignore.case = TRUE
  )
})

test_that("explicit maxReplicates at or above floor (10) triggers no warning", {
  ds30 <- .make_30tip_ds()
  expect_no_warning(
    MaximizeParsimony(ds30, maxReplicates = 10L, targetHits = 1L,
                      maxSeconds = 0.5, verbosity = 0L)
  )
})

test_that("default maxReplicates (not user-supplied) triggers no warning", {
  ds30 <- .make_30tip_ds()
  expect_no_warning(
    MaximizeParsimony(ds30, targetHits = 1L,
                      maxSeconds = 0.5, verbosity = 0L)
  )
})

test_that("no warning for nTip < 30 regardless of maxReplicates", {
  ds_small <- inapplicable.phyData[["Vinther2008"]]  # 22 taxa
  expect_no_warning(
    MaximizeParsimony(ds_small, maxReplicates = 1L, targetHits = 1L,
                      maxSeconds = 0.5, verbosity = 0L)
  )
})

# ---- Formula-driven case (min_reps > 10) -------------------------------------
# Construct a phyDat with 50 taxa and sum(weight) = 1200 so that:
#   min_reps = max(10, ceiling(50 * 1200 / 5000)) = max(10, 12) = 12

test_that("explicit maxReplicates below formula threshold triggers warning", {
  set.seed(6241)
  mat <- matrix(
    sample(c("0", "1"), 50L * 1200L, replace = TRUE),
    nrow = 50L, ncol = 1200L,
    dimnames = list(paste0("t", seq_len(50L)), NULL)
  )
  ds_large <- MatrixToPhyDat(mat)
  # Verify our assumption: sum(weight) should equal number of informative
  # patterns, which will be close to (but may be less than) 1200
  n_char <- sum(attr(ds_large, "weight"))
  min_reps <- max(10L, ceiling(50L * n_char / 5000L))

  # Test only makes sense if the formula threshold exceeds 10
  skip_if(min_reps <= 10L, "synthetic dataset too small to test formula threshold")

  expect_warning(
    MaximizeParsimony(ds_large, maxReplicates = min_reps - 1L,
                      targetHits = 1L, maxSeconds = 0.1, verbosity = 1L),
    regexp = "replicates are recommended",
    ignore.case = TRUE
  )
})

test_that("explicit maxReplicates at formula threshold triggers no warning", {
  set.seed(6241)
  mat <- matrix(
    sample(c("0", "1"), 50L * 1200L, replace = TRUE),
    nrow = 50L, ncol = 1200L,
    dimnames = list(paste0("t", seq_len(50L)), NULL)
  )
  ds_large <- MatrixToPhyDat(mat)
  n_char <- sum(attr(ds_large, "weight"))
  min_reps <- max(10L, ceiling(50L * n_char / 5000L))

  skip_if(min_reps <= 10L, "synthetic dataset too small to test formula threshold")

  expect_no_warning(
    MaximizeParsimony(ds_large, maxReplicates = min_reps,
                      targetHits = 1L, maxSeconds = 0.1, verbosity = 0L)
  )
})
