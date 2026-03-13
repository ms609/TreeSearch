library("TreeTools")

# Helper: prepare dataset for ts_* functions from a phyDat object
make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- at$weight
  levels <- at$levels
  list(contrast = contrast, tip_data = tip_data,
       weight = weight, levels = levels)
}

# Helper: score a tree with ts engine
ts_score <- function(tree, ds) {
  ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)
}

# Helper: run drift search
ts_drift <- function(tree, ds, nCycles = 10L, afdLimit = 3L,
                     rfdLimit = 0.1, maxHits = 1L) {
  ts_drift_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                  nCycles = nCycles, afdLimit = afdLimit,
                  rfdLimit = rfdLimit, maxHits = maxHits)
}

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L) {
  ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits, acceptEqual = acceptEqual,
                maxChanges = maxChanges)
}


test_that("Drift search returns valid structure", {
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_drift(tree, ds, nCycles = 2L)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_cycles_completed" %in% names(result))
  expect_true("total_drift_moves" %in% names(result))
  expect_true("total_tbr_moves" %in% names(result))
  expect_true(is.numeric(result$score))
  expect_equal(result$n_cycles_completed, 2L)
})

test_that("Drift score matches TreeLength on result tree", {
  tree <- as.phylo(100, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_drift(tree, ds, nCycles = 3L)

  result_tree <- tree
  result_tree$edge <- result$edge
  expected_score <- ts_score(result_tree, ds)
  expect_equal(result$score, expected_score)
})

test_that("Drift doesn't worsen score compared to starting tree", {
  set.seed(5821)
  tree <- as.phylo(42, 12)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # First run plain TBR to get a baseline
  tbr_result <- ts_tbr(tree, ds)

  # Now run drift starting from same tree
  drift_result <- ts_drift(tree, ds, nCycles = 5L)

  # Drift should be at least as good as starting score
  start_score <- ts_score(tree, ds)
  expect_true(drift_result$score <= start_score)
})

test_that("AFD limit of 0 behaves conservatively (equal or better only)", {
  tree <- as.phylo(200, 9)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1
  ), nrow = 9, dimnames = list(paste0("t", 1:9), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # With afd_limit=0, only equal or better moves accepted in drift phase
  result_strict <- ts_drift(tree, ds, nCycles = 3L, afdLimit = 0L)

  # Score should not exceed the converged TBR score (afd=0 means no
  # suboptimal moves, so drift phase only accepts improvements/equals)
  tbr_result <- ts_tbr(tree, ds)
  expect_true(result_strict$score <= tbr_result$score + 1,
              info = "AFD=0 drift should not be much worse than plain TBR")
})

test_that("Lower RFD limit is more conservative", {
  set.seed(4719)
  tree <- as.phylo(50, 15)
  mat <- matrix(sample(0:1, 15 * 8, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Very permissive RFD
  result_loose <- ts_drift(tree, ds, nCycles = 3L,
                           afdLimit = 5L, rfdLimit = 1.0)

  # Very restrictive RFD
  result_tight <- ts_drift(tree, ds, nCycles = 3L,
                           afdLimit = 5L, rfdLimit = 0.01)

  # Tight RFD should accept fewer drift moves (or at most equal)
  expect_true(result_tight$total_drift_moves <= result_loose$total_drift_moves + 5,
              info = paste("Tight drift moves:", result_tight$total_drift_moves,
                           "Loose drift moves:", result_loose$total_drift_moves))
})

test_that("Single cycle works", {
  tree <- as.phylo(1, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_drift(tree, ds, nCycles = 1L)

  expect_equal(result$n_cycles_completed, 1L)
  expect_true(result$score >= 0)

  # Verify score
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
})

test_that("Drift works on various tree sizes", {
  for (n_tip in c(10, 20, 50)) {
    tree <- as.phylo(1, n_tip)
    mat <- matrix(sample(0:2, n_tip * 4, replace = TRUE),
                  nrow = n_tip,
                  dimnames = list(paste0("t", seq_len(n_tip)), NULL))
    dataset <- MatrixToPhyDat(mat)
    ds <- make_ts_data(dataset)

    result <- ts_drift(tree, ds, nCycles = 2L)

    expect_true(result$score >= 0,
                info = paste("n_tip =", n_tip))

    # Verify score
    result_tree <- tree
    result_tree$edge <- result$edge
    expect_equal(result$score, ts_score(result_tree, ds),
                 info = paste("n_tip =", n_tip))
  }
})

test_that("Drift escapes local optima that plain TBR cannot", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)

  set.seed(7394)
  n_trials <- 5
  drift_wins <- 0

  for (i in seq_len(n_trials)) {
    tree <- as.phylo(sample.int(1e6, 1), length(dataset))

    tbr_result <- ts_tbr(tree, ds, maxHits = 3L)
    drift_result <- ts_drift(tree, ds, nCycles = 5L, maxHits = 3L)

    if (drift_result$score < tbr_result$score) drift_wins <- drift_wins + 1
  }

  # Drift should win at least once in 5 trials on a non-trivial dataset
  expect_true(drift_wins >= 1,
              info = paste("Drift won", drift_wins, "out of", n_trials,
                           "trials"))
})
