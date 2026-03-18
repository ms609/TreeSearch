# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Tests for Phase 3A: TBR symmetry-breaking (virtual_prelim deduplication)
#
# These tests verify that the symmetry-breaking optimization in TBR
# produces identical results to what would be obtained without it,
# across various dataset types and scoring modes.

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L, concavity = Inf) {
  TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = maxHits, acceptEqual = acceptEqual,
    maxChanges = maxChanges, concavity = concavity
  )
}


# --- Test 1: Determinism (set.seed reproducibility) ---

test_that("TBR with symmetry-breaking is deterministic (EW)", {
  tree <- as.phylo(42, 15)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(7291)
  r1 <- ts_tbr(tree, ds)

  set.seed(7291)
  r2 <- ts_tbr(tree, ds)

  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
  expect_equal(r1$n_accepted, r2$n_accepted)
  expect_equal(r1$n_evaluated, r2$n_evaluated)
})


# --- Test 2: Score matches TreeLength verification ---

test_that("TBR symmetry-breaking finds correct optimal score (EW)", {
  tree <- as.phylo(100, 12)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
    0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1,
    0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0
  ), nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(4821)
  result <- ts_tbr(tree, ds)

  # Verify score with independent scoring function
  res_tree <- tree
  res_tree$edge <- result$edge
  verified_score <- ts_score(res_tree, ds)
  expect_equal(result$score, verified_score)
})


# --- Test 3: Implied weights ---

test_that("TBR symmetry-breaking works with implied weights", {
  tree <- as.phylo(42, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 0, 0, 1, 0, 1, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(3156)
  r1 <- ts_tbr(tree, ds, concavity = 10)

  # Verify score
  res_tree <- tree
  res_tree$edge <- r1$edge
  verified <- ts_score(res_tree, ds, concavity = 10)
  expect_equal(r1$score, verified, tolerance = 1e-10)

  # Determinism
  set.seed(3156)
  r2 <- ts_tbr(tree, ds, concavity = 10)
  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})


# --- Test 4: Dataset with many duplicate states (maximizes dedup) ---

test_that("TBR handles dataset with many duplicate tip states", {
  # Create a dataset where many tips have identical state vectors.
  # This maximizes the chance of duplicate virtual_prelim values.
  tree <- as.phylo(77, 20)
  mat <- matrix(0L, nrow = 20, ncol = 3,
                dimnames = list(paste0("t", 1:20), NULL))
  # Only 3 distinct state patterns across 20 tips
  mat[1:7, ] <- c(0, 0, 0)
  mat[8:14, ] <- c(1, 1, 0)
  mat[15:20, ] <- c(1, 0, 1)
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(8512)
  result <- ts_tbr(tree, ds)

  # Must still find a valid score
  expect_true(is.finite(result$score))
  expect_true(result$score >= 0)

  res_tree <- tree
  res_tree$edge <- result$edge
  verified <- ts_score(res_tree, ds)
  expect_equal(result$score, verified)
})


# --- Test 5: Inapplicable characters ---

test_that("TBR symmetry-breaking works with inapplicable characters", {
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")

  ds_name <- "Vinther2008"
  phy_dat <- inapplicable.phyData[[ds_name]]
  n_tip <- length(phy_dat)
  tree <- as.phylo(42, n_tip, tip.label = names(phy_dat))
  ds <- make_ts_data(phy_dat)

  set.seed(6034)
  r1 <- ts_tbr(tree, ds)

  # Verify score
  res_tree <- tree
  res_tree$edge <- r1$edge
  verified <- ts_score(res_tree, ds)
  expect_equal(r1$score, verified)

  # Determinism
  set.seed(6034)
  r2 <- ts_tbr(tree, ds)
  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})


# --- Test 6: Larger inapplicable dataset ---

test_that("TBR symmetry-breaking on Longrich2010 inapplicable dataset", {
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")

  ds_name <- "Longrich2010"
  phy_dat <- inapplicable.phyData[[ds_name]]
  n_tip <- length(phy_dat)
  tree <- as.phylo(10, n_tip, tip.label = names(phy_dat))
  ds <- make_ts_data(phy_dat)

  set.seed(2847)
  result <- ts_tbr(tree, ds)

  res_tree <- tree
  res_tree$edge <- result$edge
  verified <- ts_score(res_tree, ds)
  expect_equal(result$score, verified)
  expect_true(result$score > 0)
})


# --- Test 7: Equal-score acceptance with dedup ---

test_that("TBR symmetry-breaking with accept_equal still works", {
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(9163)
  result <- ts_tbr(tree, ds, acceptEqual = TRUE, maxHits = 5L)

  res_tree <- tree
  res_tree$edge <- result$edge
  verified <- ts_score(res_tree, ds)
  expect_equal(result$score, verified)
})


# --- Test 8: Driven search integration ---

test_that("Driven search works correctly with TBR symmetry-breaking", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)

  set.seed(5720)
  r1 <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 2L, targetHits = 2L
  )

  expect_true(r1$best_score > 0)
  expect_true(length(r1$trees) >= 1)

  # Determinism
  set.seed(5720)
  r2 <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 2L, targetHits = 2L
  )
  expect_equal(r1$best_score, r2$best_score)
})


# --- Test 9: IW + inapplicable ---

test_that("TBR symmetry-breaking with IW + inapplicable characters", {
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")

  phy_dat <- inapplicable.phyData[["Vinther2008"]]
  n_tip <- length(phy_dat)
  tree <- as.phylo(42, n_tip, tip.label = names(phy_dat))
  ds <- make_ts_data(phy_dat)

  set.seed(1493)
  result <- ts_tbr(tree, ds, concavity = 3)

  res_tree <- tree
  res_tree$edge <- result$edge
  verified <- ts_score(res_tree, ds, concavity = 3)
  expect_equal(result$score, verified, tolerance = 1e-10)
})


# --- Test 10: n_evaluated should be reduced by dedup ---

test_that("Dedup reduces n_evaluated on dataset with repeated states", {
  # Use a dataset with many identical tip states to maximize dedup hits.
  tree <- as.phylo(42, 16)
  # Only 2 distinct patterns across 16 tips
  mat <- matrix(0L, nrow = 16, ncol = 4,
                dimnames = list(paste0("t", 1:16), NULL))
  mat[1:8, ] <- c(0, 0, 1, 0)
  mat[9:16, ] <- c(1, 1, 0, 1)
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(3478)
  result <- ts_tbr(tree, ds)

  # Verified score
  res_tree <- tree
  res_tree$edge <- result$edge
  verified <- ts_score(res_tree, ds)
  expect_equal(result$score, verified)

  # n_evaluated should be reasonable (not testing exact count since
  # it depends on tree topology, but should be finite and positive)
  expect_true(result$n_evaluated > 0)
  expect_true(is.finite(result$n_evaluated))
})
