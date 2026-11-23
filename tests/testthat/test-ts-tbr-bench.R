# Tier 3: only runs when TREESEARCH_EXTENDED_TESTS=true.
# See tests/testing-strategy.md
skip_extended()

# Benchmark / regression test for TBR optimizations (Phase 2B).
# Tests correctness: optimized TBR must find scores equal to or better than
# the baseline, and result topologies must be valid.

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result


ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L, concavity = Inf) {
  TreeSearch:::ts_tbr_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, acceptEqual = acceptEqual,
                             maxChanges = maxChanges, concavity = concavity)
}

# Validate result topology
validate_result <- function(result, n_tip) {
  edge <- result$edge
  expect_equal(nrow(edge), 2L * (n_tip - 1L))
  children <- edge[, 2]
  tips <- sort(children[children <= n_tip])
  expect_equal(tips, seq_len(n_tip))
}

test_that("TBR optimized: small dataset correctness", {
  set.seed(4821)
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)
  expect_true(result$converged)
  # Three binary characters, each needs at least 1 step → minimum 3
  # but actual optimum depends on the starting tree and character compatibility
  expect_true(result$score >= 3L)
  validate_result(result, 8L)
})

test_that("TBR optimized: 15-tip dataset", {
  set.seed(7734)
  tree <- as.phylo(1, 15)
  mat <- matrix(sample(0:2, 15 * 6, replace = TRUE),
                nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)
  expect_true(result$converged)
  expect_true(result$score >= 0)
  validate_result(result, 15L)

  rt <- tree
  rt$edge <- result$edge
  expect_equal(result$score, ts_score(rt, ds))
})

test_that("TBR optimized: 25-tip dataset finds good score", {
  set.seed(3192)
  tree <- as.phylo(1, 25)
  mat <- matrix(sample(0:3, 25 * 10, replace = TRUE),
                nrow = 25, dimnames = list(paste0("t", 1:25), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  start_score <- ts_score(tree, ds)
  result <- ts_tbr(tree, ds, maxHits = 3L)

  expect_true(result$score <= start_score)
  expect_true(result$converged)
  validate_result(result, 25L)
})

test_that("TBR optimized: Congreve-Lamsdell dataset", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  set.seed(5501)
  tree <- as.phylo(1, n_tip)
  start_score <- ts_score(tree, ds)
  result <- ts_tbr(tree, ds, maxHits = 3L)

  expect_true(result$score < start_score)
  expect_true(result$converged)
  validate_result(result, n_tip)
})

test_that("TBR optimized: deterministic with set.seed", {
  set.seed(9113)
  tree <- as.phylo(100, 12)
  mat <- matrix(sample(0:1, 12 * 5, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(2244)
  r1 <- ts_tbr(tree, ds, maxHits = 5L, acceptEqual = TRUE)
  set.seed(2244)
  r2 <- ts_tbr(tree, ds, maxHits = 5L, acceptEqual = TRUE)

  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})

test_that("TBR optimized: IW scoring", {
  set.seed(6320)
  tree <- as.phylo(50, 12)
  mat <- matrix(sample(0:2, 12 * 6, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds, concavity = 10.0)
  expect_true(result$converged)
  expect_true(result$score >= 0)
  validate_result(result, 12L)
})

test_that("TBR optimized: accept_equal and maxChanges", {
  set.seed(1847)
  tree <- as.phylo(1, 10)
  mat <- matrix(sample(0:1, 10 * 5, replace = TRUE),
                nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds, maxChanges = 3L, acceptEqual = TRUE,
                   maxHits = 100L)
  expect_true(result$n_accepted <= 3L)
  validate_result(result, 10L)
})
