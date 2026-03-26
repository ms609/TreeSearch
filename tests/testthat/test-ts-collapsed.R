# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L) {
  TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = maxHits, acceptEqual = acceptEqual,
    maxChanges = maxChanges
  )
}

test_that("TBR result includes n_zero_skipped field", {
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)

  expect_true("n_zero_skipped" %in% names(result))
  expect_true(is.numeric(result$n_zero_skipped))
  expect_gte(result$n_zero_skipped, 0)
})

test_that("Score equivalence: TBR finds same optima with collapsed flags", {
  # Informative 10-tip dataset: enough signal to have a unique optimum
  tree <- as.phylo(100, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 1, 1, 1, 1, 0, 0, 1,
    1, 1, 0, 0, 1, 0, 0, 1, 1, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds, maxHits = 5L)

  result_tree <- tree
  result_tree$edge <- result$edge
  full_score <- ts_score(result_tree, ds)
  expect_equal(result$score, full_score)
})

test_that("Sparse data produces zero-length edges (non-zero skip count)", {
  # Construct a dataset where most tips are identical → many zero-length
  # internal edges at the optimum.
  n_tip <- 15
  mat <- matrix(0L, nrow = n_tip, ncol = 2,
                dimnames = list(paste0("t", seq_len(n_tip)), NULL))
  # Only 2 tips differ from the majority

  mat[1, ] <- c(1, 0)
  mat[2, ] <- c(0, 1)
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Find optimal tree first
  set.seed(3714)
  tree <- as.phylo(1, n_tip)
  result <- ts_tbr(tree, ds, maxHits = 10L)

  # Optimal tree for this data has score = 2 (one step per character)
  expect_lte(result$score, 2)

  # Now run TBR again from the already-optimal tree
  opt_tree <- tree
  opt_tree$edge <- result$edge
  result2 <- ts_tbr(opt_tree, ds, maxHits = 10L)

  # At a converged tree with many identical tips, most edges are zero-length
  # and the collapsed optimization should skip them
  expect_gt(result2$n_zero_skipped, 0)
})

test_that("Collapsed flags work with inapplicable characters", {
  tree <- PectinateTree(8)
  tree$tip.label <- paste0("t", 1:8)

  # Dataset with inapplicable characters but most tips share the same state
  mat <- matrix(c(
    "1", "1", "1", "1", "1", "1", "1", "1",
    "-", "-", "-", "-", "1", "1", "1", "2"
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- TreeTools::MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds, maxHits = 5L)

  result_tree <- tree
  result_tree$edge <- result$edge
  full_score <- ts_score(result_tree, ds)
  expect_equal(result$score, full_score)
})

test_that("Collapsed flags work with implied weighting", {
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  min_steps <- TreeSearch:::ts_char_steps(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels
  )
  # Sanity: min_steps must have correct length
  expect_equal(length(min_steps), length(ds$weight))

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = 5L, min_steps = min_steps, concavity = 3.0
  )

  expect_true("n_zero_skipped" %in% names(result))
  expect_gte(result$n_zero_skipped, 0)
  expect_true(is.numeric(result$score))
})

test_that("Regraft merging: sparse data search succeeds with region skipping", {
  # Dataset with many identical tips → large collapsed regions.
  # Search should find the optimum despite skipping interior regraft positions.
  n_tip <- 20
  mat <- matrix(0L, nrow = n_tip, ncol = 4,
                dimnames = list(paste0("t", seq_len(n_tip)), NULL))
  mat[1, ] <- c(1, 0, 0, 0)
  mat[2, ] <- c(0, 1, 0, 0)
  mat[3, ] <- c(0, 0, 1, 0)
  mat[4, ] <- c(0, 0, 0, 1)
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(6714)
  tree <- as.phylo(1, n_tip)
  result <- ts_tbr(tree, ds, maxHits = 10L)

  # Optimal score = 4 (one step per character)
  expect_equal(result$score, 4)

  # Converged tree: many collapsed edges → non-trivial skip count
  opt_tree <- tree
  opt_tree$edge <- result$edge
  result2 <- ts_tbr(opt_tree, ds, maxHits = 10L)
  expect_gt(result2$n_zero_skipped, 0)
})

test_that("Collapsed pool dedup: driven search works with collapsed dedup", {
  # Small dataset: run driven search and verify it completes without error.
  # The collapsed pool dedup is exercised in the driven pipeline.
  ds <- make_ts_data(TreeSearch::inapplicable.phyData[["Vinther2008"]])

  set.seed(7192)
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 3L, targetHits = 1L,
    ratchetCycles = 1L, driftCycles = 1L,
    xssRounds = 0L, xssPartitions = 3L,
    rssRounds = 0L, cssRounds = 0L,
    fuseInterval = 3L, maxSeconds = 15,
    verbosity = 0L, nThreads = 1L
  )

  expect_true(result$best_score <= 85)
  expect_true(result$pool_size >= 1)
  validate_result(result, length(TreeSearch::inapplicable.phyData[["Vinther2008"]]))
})

test_that("SPR search works with collapsed regraft skipping", {
  ts_spr <- function(tree, ds, maxHits = 1L) {
    TreeSearch:::ts_spr_search(
      tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxHits = maxHits
    )
  }

  n_tip <- 15
  mat <- matrix(0L, nrow = n_tip, ncol = 3,
                dimnames = list(paste0("t", seq_len(n_tip)), NULL))
  mat[1, ] <- c(1, 0, 0)
  mat[2, ] <- c(0, 1, 0)
  mat[3, ] <- c(0, 0, 1)
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(2934)
  tree <- as.phylo(1, n_tip)
  result <- ts_spr(tree, ds, maxHits = 10L)

  # Optimal score = 3
  expect_equal(result$score, 3)
})

test_that("Drift search works with collapsed regraft skipping", {
  ds <- make_ts_data(TreeSearch::inapplicable.phyData[["Vinther2008"]])

  tree <- as.phylo(42, length(TreeSearch::inapplicable.phyData[["Vinther2008"]]))
  tree$tip.label <- names(TreeSearch::inapplicable.phyData[["Vinther2008"]])

  set.seed(8451)
  result <- TreeSearch:::ts_drift_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 20L, afdLimit = 5L, rfdLimit = 0.2
  )

  expect_true(is.numeric(result$score))
  expect_true(result$score < Inf)
})

test_that("Driven search produces valid results with collapsed flags active", {
  ds <- make_ts_data(TreeSearch::inapplicable.phyData[["Vinther2008"]])

  set.seed(5839)
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 2L, targetHits = 1L,
    ratchetCycles = 2L, driftCycles = 1L,
    xssRounds = 1L, xssPartitions = 3L,
    rssRounds = 0L, cssRounds = 0L,
    fuseInterval = 3L, maxSeconds = 15,
    verbosity = 0L, nThreads = 1L
  )

  # Allow small margin: 2 replicates may not reach global optimum (79)
  expect_true(result$best_score <= 85)
  expect_true(result$pool_size >= 1)
  validate_result(result, length(TreeSearch::inapplicable.phyData[["Vinther2008"]]))
})
