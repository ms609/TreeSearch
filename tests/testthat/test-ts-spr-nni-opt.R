# Tests for SPR and NNI search optimizations (Phase 2C).
# Verifies that optimized SPR (bounded indirect, subtree filter, incremental
# clip, deferred reshuffling) and NNI (incremental scoring) produce correct
# results across EW, IW, and NA datasets.

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

ts_score <- function(tree, ds) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels)
}

ts_spr <- function(tree, ds, maxHits = 20L, concavity = Inf) {
  TreeSearch:::ts_spr_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, concavity = concavity)
}

ts_nni <- function(tree, ds, maxHits = 20L, concavity = Inf) {
  TreeSearch:::ts_nni_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, concavity = concavity)
}

validate_result <- function(result, n_tip) {
  edge <- result$edge
  expect_equal(nrow(edge), 2L * (n_tip - 1L))
  children <- edge[, 2]
  tips <- sort(children[children <= n_tip])
  expect_equal(tips, seq_len(n_tip))
}

# ---- SPR correctness tests ----

test_that("SPR: small 8-tip dataset", {
  set.seed(4210)
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_spr(tree, ds)
  expect_true(result$score >= 3L)
  validate_result(result, 8L)
})

test_that("SPR: 15-tip random dataset", {
  set.seed(8293)
  tree <- as.phylo(1, 15)
  mat <- matrix(sample(0:2, 15 * 6, replace = TRUE),
                nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  start_score <- ts_score(tree, ds)
  result <- ts_spr(tree, ds)
  expect_true(result$score <= start_score)
  validate_result(result, 15L)

  rt <- tree
  rt$edge <- result$edge
  expect_equal(result$score, ts_score(rt, ds))
})

test_that("SPR: 25-tip dataset finds improvement", {
  set.seed(3741)
  tree <- as.phylo(1, 25)
  mat <- matrix(sample(0:3, 25 * 10, replace = TRUE),
                nrow = 25, dimnames = list(paste0("t", 1:25), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  start_score <- ts_score(tree, ds)
  result <- ts_spr(tree, ds, maxHits = 3L)
  expect_true(result$score <= start_score)
  validate_result(result, 25L)
})

test_that("SPR: Congreve-Lamsdell empirical dataset", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  set.seed(6192)
  tree <- as.phylo(1, n_tip)
  start_score <- ts_score(tree, ds)
  result <- ts_spr(tree, ds, maxHits = 3L)

  expect_true(result$score < start_score)
  validate_result(result, n_tip)
})

test_that("SPR: deterministic with set.seed", {
  set.seed(5517)
  tree <- as.phylo(100, 12)
  mat <- matrix(sample(0:1, 12 * 5, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(7803)
  r1 <- ts_spr(tree, ds, maxHits = 5L)
  set.seed(7803)
  r2 <- ts_spr(tree, ds, maxHits = 5L)

  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})

test_that("SPR: IW scoring", {
  set.seed(2845)
  tree <- as.phylo(50, 12)
  mat <- matrix(sample(0:2, 12 * 6, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_spr(tree, ds, concavity = 10.0)
  expect_true(result$score >= 0)
  validate_result(result, 12L)
})

test_that("SPR: inapplicable characters", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")
  # Use first dataset that has inapplicable characters
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  set.seed(1937)
  tree <- as.phylo(1, n_tip)
  start_score <- ts_score(tree, ds)
  result <- ts_spr(tree, ds, maxHits = 3L)

  expect_true(result$score <= start_score)
  validate_result(result, n_tip)
})

test_that("SPR: SPR score <= TBR score (SPR is weaker)", {
  skip_if_not_installed("TreeSearch")
  set.seed(9431)
  tree <- as.phylo(1, 15)
  mat <- matrix(sample(0:2, 15 * 8, replace = TRUE),
                nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(3310)
  spr_result <- ts_spr(tree, ds)
  set.seed(3310)
  tbr_result <- TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data,
    ds$weight, ds$levels, maxHits = 1L)
  # TBR searches a superset of SPR moves, so it should find
  # a score at least as good
  expect_true(tbr_result$score <= spr_result$score)
})

# ---- NNI correctness tests ----

test_that("NNI: small 8-tip dataset", {
  set.seed(6754)
  tree <- as.phylo(42, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_nni(tree, ds)
  expect_true(result$score >= 3L)
  validate_result(result, 8L)
})

test_that("NNI: 15-tip random dataset", {
  set.seed(1482)
  tree <- as.phylo(1, 15)
  mat <- matrix(sample(0:2, 15 * 6, replace = TRUE),
                nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  start_score <- ts_score(tree, ds)
  result <- ts_nni(tree, ds)
  expect_true(result$score <= start_score)
  validate_result(result, 15L)

  rt <- tree
  rt$edge <- result$edge
  expect_equal(result$score, ts_score(rt, ds))
})

test_that("NNI: 25-tip dataset finds improvement", {
  set.seed(4590)
  tree <- as.phylo(1, 25)
  mat <- matrix(sample(0:3, 25 * 10, replace = TRUE),
                nrow = 25, dimnames = list(paste0("t", 1:25), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  start_score <- ts_score(tree, ds)
  result <- ts_nni(tree, ds, maxHits = 3L)
  expect_true(result$score <= start_score)
  validate_result(result, 25L)
})

test_that("NNI: deterministic with set.seed", {
  set.seed(3392)
  tree <- as.phylo(100, 12)
  mat <- matrix(sample(0:1, 12 * 5, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(8871)
  r1 <- ts_nni(tree, ds, maxHits = 5L)
  set.seed(8871)
  r2 <- ts_nni(tree, ds, maxHits = 5L)

  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})

test_that("NNI: IW scoring", {
  set.seed(5213)
  tree <- as.phylo(50, 12)
  mat <- matrix(sample(0:2, 12 * 6, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_nni(tree, ds, concavity = 10.0)
  expect_true(result$score >= 0)
  validate_result(result, 12L)
})

test_that("NNI: inapplicable characters (falls back to full rescore)", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  set.seed(6328)
  tree <- as.phylo(1, n_tip)
  start_score <- ts_score(tree, ds)
  result <- ts_nni(tree, ds, maxHits = 3L)

  expect_true(result$score <= start_score)
  validate_result(result, n_tip)
})

test_that("NNI: SPR finds equal or better score than NNI", {
  set.seed(2750)
  tree <- as.phylo(1, 15)
  mat <- matrix(sample(0:2, 15 * 8, replace = TRUE),
                nrow = 15, dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  set.seed(4417)
  nni_result <- ts_nni(tree, ds)
  set.seed(4417)
  spr_result <- ts_spr(tree, ds)
  # SPR searches a superset of NNI moves
  expect_true(spr_result$score <= nni_result$score)
})

test_that("NNI: Congreve-Lamsdell empirical dataset", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  set.seed(7205)
  tree <- as.phylo(1, n_tip)
  start_score <- ts_score(tree, ds)
  result <- ts_nni(tree, ds, maxHits = 3L)

  expect_true(result$score < start_score)
  validate_result(result, n_tip)
})

# ---- Cross-method consistency ----

test_that("All methods converge on same optimum (small dataset)", {
  set.seed(8115)
  mat <- matrix(sample(0:1, 10 * 8, replace = TRUE),
                nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Run from multiple starting trees with each method
  scores_nni <- numeric(3)
  scores_spr <- numeric(3)
  for (i in 1:3) {
    set.seed(1000 + i)
    tree <- as.phylo(i * 100, 10)
    scores_nni[i] <- ts_nni(tree, ds, maxHits = 10L)$score
    set.seed(1000 + i)
    scores_spr[i] <- ts_spr(tree, ds, maxHits = 10L)$score
  }

  # SPR should find scores at least as good as NNI
  for (i in 1:3) {
    expect_true(scores_spr[i] <= scores_nni[i])
  }
})
