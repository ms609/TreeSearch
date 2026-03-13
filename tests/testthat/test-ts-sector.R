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
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                               ds$weight, ds$levels)
}

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L) {
  TreeSearch:::ts_tbr_search(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels, maxHits = maxHits)
}

# Helper: run RSS search
ts_rss <- function(tree, ds, minSize = 6L, maxSize = 50L,
                   acceptEqual = FALSE, rssPicks = 0L,
                   ratchetCycles = 6L, maxHits = 1L) {
  TreeSearch:::ts_rss_search(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels,
                              minSectorSize = minSize,
                              maxSectorSize = maxSize,
                              acceptEqual = acceptEqual,
                              rssPicks = rssPicks,
                              ratchetCycles = ratchetCycles,
                              maxHits = maxHits)
}

# Helper: run XSS search
ts_xss <- function(tree, ds, nPartitions = 4L, xssRounds = 3L,
                   acceptEqual = FALSE, ratchetCycles = 6L,
                   maxHits = 1L) {
  TreeSearch:::ts_xss_search(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels,
                              nPartitions = nPartitions,
                              xssRounds = xssRounds,
                              acceptEqual = acceptEqual,
                              ratchetCycles = ratchetCycles,
                              maxHits = maxHits)
}

# ---------- Test datasets ----------

# Small dataset: 10 tips, 4 characters
small_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
small_dataset <- MatrixToPhyDat(small_mat)
small_ds <- make_ts_data(small_dataset)

# Medium dataset: 30 tips with more conflict
set.seed(5471)
med_mat <- matrix(sample(0:1, 30 * 15, replace = TRUE),
                  nrow = 30,
                  dimnames = list(paste0("t", 1:30), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)


test_that("RSS returns valid structure", {
  tree <- as.phylo(42, 10)
  result <- ts_rss(tree, small_ds, minSize = 4L, maxSize = 8L,
                   rssPicks = 3L, ratchetCycles = 0L)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_sectors_searched" %in% names(result))
  expect_true("n_sectors_improved" %in% names(result))
  expect_true("total_steps_saved" %in% names(result))
})

test_that("XSS returns valid structure", {
  tree <- as.phylo(42, 10)
  result <- ts_xss(tree, small_ds, nPartitions = 2L, xssRounds = 1L,
                   ratchetCycles = 0L)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_sectors_searched" %in% names(result))
})

test_that("RSS score matches independent verification", {
  tree <- as.phylo(100, 10)
  result <- ts_rss(tree, small_ds, minSize = 4L, maxSize = 8L,
                   rssPicks = 5L, ratchetCycles = 0L)

  result_tree <- tree
  result_tree$edge <- result$edge
  expected_score <- ts_score(result_tree, small_ds)
  expect_equal(result$score, expected_score)
})

test_that("XSS score matches independent verification", {
  tree <- as.phylo(100, 10)
  result <- ts_xss(tree, small_ds, nPartitions = 2L, xssRounds = 1L,
                   ratchetCycles = 0L)

  result_tree <- tree
  result_tree$edge <- result$edge
  expected_score <- ts_score(result_tree, small_ds)
  expect_equal(result$score, expected_score)
})

test_that("RSS does not worsen score", {
  tree <- as.phylo(42, 10)
  initial_score <- ts_score(tree, small_ds)

  result <- ts_rss(tree, small_ds, minSize = 4L, maxSize = 8L,
                   rssPicks = 5L, ratchetCycles = 0L)

  expect_true(result$score <= initial_score)
})

test_that("XSS does not worsen score", {
  tree <- as.phylo(42, 10)
  initial_score <- ts_score(tree, small_ds)

  result <- ts_xss(tree, small_ds, nPartitions = 2L, xssRounds = 1L,
                   ratchetCycles = 0L)

  expect_true(result$score <= initial_score)
})

test_that("RSS on 30-tip data produces valid result", {
  tree <- as.phylo(1234, 30)
  initial_score <- ts_score(tree, med_ds)

  rss_result <- ts_rss(tree, med_ds, minSize = 6L, maxSize = 20L,
                       rssPicks = 10L, ratchetCycles = 0L)

  # RSS should not worsen compared to start (includes global TBR)
  expect_true(rss_result$score <= initial_score)
  expect_true(rss_result$score > 0)

  # Verify returned tree matches reported score
  result_tree <- tree
  result_tree$edge <- rss_result$edge
  expect_equal(rss_result$score, ts_score(result_tree, med_ds))
})

test_that("XSS on 30-tip data produces valid result", {
  tree <- as.phylo(1234, 30)
  initial_score <- ts_score(tree, med_ds)

  xss_result <- ts_xss(tree, med_ds, nPartitions = 3L, xssRounds = 2L,
                       ratchetCycles = 0L)

  expect_true(xss_result$score <= initial_score)
  expect_true(xss_result$score > 0)

  result_tree <- tree
  result_tree$edge <- xss_result$edge
  expect_equal(xss_result$score, ts_score(result_tree, med_ds))
})

test_that("RSS works with various tree sizes", {
  for (n in c(8, 15, 20)) {
    mat <- matrix(sample(0:1, n * 5, replace = TRUE),
                  nrow = n,
                  dimnames = list(paste0("t", seq_len(n)), NULL))
    dat <- MatrixToPhyDat(mat)
    ds <- make_ts_data(dat)
    tree <- as.phylo(1, n)

    result <- ts_rss(tree, ds, minSize = 4L,
                     maxSize = min(n - 2L, 12L),
                     rssPicks = 3L, ratchetCycles = 0L)

    expect_true(result$score > 0)

    result_tree <- tree
    result_tree$edge <- result$edge
    expect_equal(result$score, ts_score(result_tree, ds))
  }
})

test_that("XSS works with various partition counts", {
  tree <- as.phylo(100, 20)
  mat <- matrix(sample(0:1, 20 * 8, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dat <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dat)

  for (p in c(2, 3, 4)) {
    result <- ts_xss(tree, ds, nPartitions = p, xssRounds = 1L,
                     ratchetCycles = 0L)
    expect_true(result$score > 0)
    expect_true(result$n_sectors_searched >= 1)
  }
})
