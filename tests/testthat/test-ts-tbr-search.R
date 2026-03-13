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

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L) {
  ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits, acceptEqual = acceptEqual,
                maxChanges = maxChanges)
}

# Helper: run SPR search
ts_spr <- function(tree, ds, maxHits = 20L) {
  ts_spr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits)
}


test_that("TBR search returns valid structure", {
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_accepted" %in% names(result))
  expect_true("n_evaluated" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true(is.integer(result$score) || is.numeric(result$score))
  expect_true(result$converged)
})

test_that("TBR score matches TreeLength on result tree", {
  tree <- as.phylo(100, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)

  # Verify reported score matches independent calculation
  result_tree <- tree
  result_tree$edge <- result$edge
  expected_score <- ts_score(result_tree, ds)
  expect_equal(result$score, expected_score)
})

test_that("TBR finds competitive scores compared to SPR", {
  set.seed(8347)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tbr_wins <- 0
  spr_wins <- 0
  for (i in 1:10) {
    tree <- as.phylo(sample.int(1e6, 1), 12)
    spr_result <- ts_spr(tree, ds)
    tbr_result <- ts_tbr(tree, ds)
    if (tbr_result$score < spr_result$score) tbr_wins <- tbr_wins + 1
    if (spr_result$score < tbr_result$score) spr_wins <- spr_wins + 1
  }

  # TBR should not be consistently worse than SPR
  expect_true(spr_wins <= tbr_wins + 3,
              info = paste("SPR won", spr_wins, "times, TBR won", tbr_wins))
})

test_that("TBR finds optimal score on small known cases", {
  # 7-tip tree with a single informative character: optimal = 1 step
  tree <- as.phylo(1, 7)
  mat <- matrix(c(0, 0, 0, 0, 1, 1, 1),
                nrow = 7,
                dimnames = list(paste0("t", 1:7), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)
  expect_equal(result$score, 1L)
})

test_that("TBR accept_equal allows lateral moves", {
  tree <- as.phylo(200, 9)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1
  ), nrow = 9, dimnames = list(paste0("t", 1:9), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result_strict <- ts_tbr(tree, ds, acceptEqual = FALSE)
  result_equal  <- ts_tbr(tree, ds, acceptEqual = TRUE, maxHits = 5L)

  # With equal acceptance, should accept at least as many moves
  # (or same number if already at optimum)
  expect_true(result_equal$n_accepted >= result_strict$n_accepted ||
              result_equal$score <= result_strict$score)
})

test_that("TBR max_accepted_changes stops early", {
  tree <- as.phylo(1, 10)
  mat <- matrix(sample(0:1, 10 * 5, replace = TRUE),
                nrow = 10,
                dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds, maxChanges = 2L, acceptEqual = TRUE,
                   maxHits = 100L)

  # Should have accepted at most 2 changes

  expect_true(result$n_accepted <= 2L)
})

test_that("TBR works on various tree sizes", {
  for (n_tip in c(5, 8, 15, 25)) {
    tree <- as.phylo(1, n_tip)
    mat <- matrix(sample(0:2, n_tip * 4, replace = TRUE),
                  nrow = n_tip,
                  dimnames = list(paste0("t", seq_len(n_tip)), NULL))
    dataset <- MatrixToPhyDat(mat)
    ds <- make_ts_data(dataset)

    result <- ts_tbr(tree, ds)

    # Score should be positive and reasonable
    expect_true(result$score >= 0)
    expect_true(result$converged)

    # Verify score
    result_tree <- tree
    result_tree$edge <- result$edge
    expect_equal(result$score, ts_score(result_tree, ds))
  }
})

test_that("TBR result tree has valid topology", {
  tree <- as.phylo(100, 12)
  mat <- matrix(sample(0:1, 12 * 5, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_tbr(tree, ds)

  # Build a phylo object and check it's valid
  result_tree <- tree
  result_tree$edge <- result$edge
  result_tree$Nnode <- tree$Nnode

  # Should have the right number of edges
  expect_equal(nrow(result$edge), 2 * (length(tree$tip.label) - 1))

  # All tips should be present
  tips_in_tree <- sort(result$edge[result$edge[, 2] <= length(tree$tip.label), 2])
  expect_equal(tips_in_tree, seq_len(length(tree$tip.label)))
})

test_that("TBR on Congreve-Lamsdell dataset improves over random start", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)

  tree <- as.phylo(1, length(dataset))
  start_score <- ts_score(tree, ds)
  result <- ts_tbr(tree, ds, maxHits = 3L)

  expect_true(result$score <= start_score)
  expect_true(result$score < start_score,
              info = "TBR should improve a random starting tree")
})
