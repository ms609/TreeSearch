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

# Helper: run ratchet search
ts_ratchet <- function(tree, ds, nCycles = 10L, perturbProb = 0.04,
                       maxHits = 1L) {
  ts_ratchet_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                    nCycles = nCycles, perturbProb = perturbProb,
                    maxHits = maxHits)
}

# Helper: run TBR search
ts_tbr <- function(tree, ds, maxHits = 1L) {
  ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits)
}


test_that("Ratchet search returns valid structure", {
  tree <- as.phylo(42, 10)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
    0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
    0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
    1, 0, 0, 0, 1, 1, 1, 1, 0, 0
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 2L)

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("n_cycles" %in% names(result))
  expect_true("total_tbr_moves" %in% names(result))
  expect_true(result$score > 0)
  expect_equal(result$n_cycles, 2L)
})

test_that("Ratchet score matches TreeLength on result tree", {
  tree <- as.phylo(100, 12)
  set.seed(7341)
  mat <- matrix(sample(0:1, 12 * 8, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 3L)

  result_tree <- tree
  result_tree$edge <- result$edge
  expected_score <- ts_score(result_tree, ds)
  expect_equal(result$score, expected_score)
})

test_that("Ratchet does not worsen score vs plain TBR", {
  set.seed(5612)
  mat <- matrix(sample(0:1, 15 * 10, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tree <- as.phylo(1, 15)
  tbr_result <- ts_tbr(tree, ds)
  ratchet_result <- ts_ratchet(tree, ds, nCycles = 5L)

  expect_true(ratchet_result$score <= tbr_result$score)
})

test_that("Ratchet escapes local optima on Congreve-Lamsdell data", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)

  # Start from several random trees; compare TBR-only vs ratchet
  set.seed(4821)
  ratchet_wins <- 0
  n_trials <- 5
  for (i in seq_len(n_trials)) {
    tree <- as.phylo(sample.int(1e6, 1), length(dataset))
    tbr_result <- ts_tbr(tree, ds, maxHits = 3L)
    ratchet_result <- ts_ratchet(tree, ds, nCycles = 10L, maxHits = 3L)
    if (ratchet_result$score < tbr_result$score) {
      ratchet_wins <- ratchet_wins + 1
    }
  }
  # Ratchet should beat plain TBR at least once
  expect_true(ratchet_wins >= 1,
              info = paste("Ratchet won", ratchet_wins, "/", n_trials,
                           "trials"))
})

test_that("Single cycle completes without error", {
  tree <- as.phylo(1, 8)
  mat <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1,
    0, 0, 1, 1, 0, 0, 1, 1,
    0, 1, 0, 1, 0, 1, 0, 1
  ), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 1L)

  expect_equal(result$n_cycles, 1L)
  expect_true(result$score > 0)
})

test_that("Higher perturbation probability changes search behavior", {
  set.seed(9183)
  mat <- matrix(sample(0:2, 20 * 15, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tree <- as.phylo(1, 20)

  result_low  <- ts_ratchet(tree, ds, nCycles = 5L, perturbProb = 0.02)
  result_high <- ts_ratchet(tree, ds, nCycles = 5L, perturbProb = 0.20)

  # Both should return valid scores; high perturbation typically
  # produces more TBR moves during the perturbation phase
  expect_true(result_low$score > 0)
  expect_true(result_high$score > 0)
  expect_true(result_high$total_tbr_moves != result_low$total_tbr_moves ||
              result_high$score != result_low$score,
              info = "Different perturbation rates should produce different behavior")
})

test_that("Ratchet result tree has valid topology", {
  tree <- as.phylo(100, 14)
  set.seed(2765)
  mat <- matrix(sample(0:1, 14 * 6, replace = TRUE),
                nrow = 14,
                dimnames = list(paste0("t", 1:14), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 3L)

  # Right number of edges
  expect_equal(nrow(result$edge), 2 * (length(tree$tip.label) - 1))

  # All tips present
  tips_in_tree <- sort(
    result$edge[result$edge[, 2] <= length(tree$tip.label), 2]
  )
  expect_equal(tips_in_tree, seq_len(length(tree$tip.label)))
})
