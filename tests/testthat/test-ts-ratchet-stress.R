library("TreeTools")

# Stress tests for the C++ ratchet implementation (Agent A).
# These go beyond the basic tests to probe edge cases, correctness
# invariants, and robustness under adversarial inputs.

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
  ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)
}

ts_ratchet <- function(tree, ds, nCycles = 10L, perturbProb = 0.04,
                       maxHits = 1L) {
  ts_ratchet_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                    nCycles = nCycles, perturbProb = perturbProb,
                    maxHits = maxHits)
}

ts_tbr <- function(tree, ds, maxHits = 1L, acceptEqual = FALSE,
                   maxChanges = 0L) {
  ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                maxHits = maxHits, acceptEqual = acceptEqual,
                maxChanges = maxChanges)
}


# --- 1. Score integrity: reported score always matches independent rescore ---

test_that("Ratchet score exactly matches independent rescore (20 random starts)", {
  set.seed(6482)
  mat <- matrix(sample(0:1, 15 * 8, replace = TRUE),
                nrow = 15,
                dimnames = list(paste0("t", 1:15), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (i in 1:20) {
    tree <- as.phylo(sample.int(1e6, 1), 15)
    result <- ts_ratchet(tree, ds, nCycles = 3L)
    result_tree <- tree
    result_tree$edge <- result$edge
    actual <- ts_score(result_tree, ds)
    expect_equal(result$score, actual,
                 info = paste("Trial", i, ": reported", result$score,
                              "vs actual", actual))
  }
})


# --- 2. Ratchet output is a local optimum ---
# The ratchet's final search phase should converge, so a fresh TBR on the
# output tree should not improve the score. (We don't compare against a
# separate TBR call, since different RNG seeds find different optima.)

test_that("Ratchet output is a TBR local optimum", {
  set.seed(3917)
  mat <- matrix(sample(0:2, 20 * 12, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (i in 1:10) {
    tree <- as.phylo(sample.int(1e6, 1), 20)
    ratchet_result <- ts_ratchet(tree, ds, nCycles = 5L)

    # Run TBR on the ratchet's output — should not improve
    tbr_on_output <- ts_tbr_search(
      ratchet_result$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)

    expect_equal(tbr_on_output$score, ratchet_result$score,
                 info = paste("Trial", i, ": ratchet", ratchet_result$score,
                              "but TBR improved to", tbr_on_output$score))
  }
})

# --- 2b. Ratchet should not worsen vs starting tree score ---

test_that("Ratchet final score <= starting tree score", {
  set.seed(3917)
  mat <- matrix(sample(0:2, 20 * 12, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (i in 1:10) {
    tree <- as.phylo(sample.int(1e6, 1), 20)
    start_score <- ts_score(tree, ds)
    ratchet_result <- ts_ratchet(tree, ds, nCycles = 5L)
    expect_true(ratchet_result$score <= start_score,
                info = paste("Trial", i, ": start", start_score,
                             "ratchet", ratchet_result$score))
  }
})


# --- 3. Edge case: minimum tree size (5 tips) ---

test_that("Ratchet works on minimum-size trees (5 tips)", {
  tree <- as.phylo(1, 5)
  mat <- matrix(c(0, 0, 0, 1, 1,
                  0, 1, 1, 0, 0),
                nrow = 5,
                dimnames = list(paste0("t", 1:5), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 5L)
  expect_true(result$score >= 0)
  expect_equal(result$n_cycles, 5L)

  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
})


# --- 4. Edge case: single character ---

test_that("Ratchet works with a single binary character", {
  tree <- as.phylo(1, 10)
  mat <- matrix(c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
                nrow = 10,
                dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 3L)
  # Optimal score for a single binary character is 1
  expect_equal(result$score, 1)
})


# --- 5. Edge case: perturbProb = 0 (no perturbation) ---
# Should behave identically to plain TBR.

test_that("Ratchet with perturbProb=0 matches plain TBR", {
  tree <- as.phylo(42, 12)
  mat <- matrix(c(
    0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1,
    0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1
  ), nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tbr_result <- ts_tbr(tree, ds)
  ratchet_result <- ts_ratchet(tree, ds, nCycles = 3L, perturbProb = 0)

  # With zero perturbation, the perturbation phase is just TBR with
  # accept_equal=true. The search phase then re-converges. Final
  # score should be <= TBR score.
  expect_true(ratchet_result$score <= tbr_result$score)
})


# --- 6. Edge case: perturbProb = 1.0 (all characters zeroed) ---
# All active_mask bits zeroed => score becomes 0 during perturbation.
# The perturbation TBR should accept anything. Search phase should
# recover a meaningful score.

test_that("Ratchet with perturbProb=1 doesn't crash and recovers", {
  tree <- as.phylo(1, 10)
  set.seed(8150)
  mat <- matrix(sample(0:1, 10 * 5, replace = TRUE),
                nrow = 10,
                dimnames = list(paste0("t", 1:10), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 3L, perturbProb = 1.0)

  # Should still return a valid score (not 0, since search phase
  # uses original weights)
  expect_true(result$score > 0)
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
})


# --- 7. Weights not corrupted: dataset usable after ratchet ---
# The ratchet modifies active_mask internally but should restore it.
# Verify by scoring the same tree before and after ratchet call.

test_that("Dataset active_masks are restored after ratchet", {
  tree <- as.phylo(100, 12)
  set.seed(2243)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  # Score a known tree before ratchet
  score_before <- ts_score(tree, ds)

  # Run ratchet (which internally perturbs and restores active_mask)
  result <- ts_ratchet(tree, ds, nCycles = 5L, perturbProb = 0.2)

  # Score the SAME original tree again — should get the same score
  # This only works because ts_fitch_score creates a fresh DataSet
  # from the R arguments each time. But if the ratchet were to corrupt
  # the R-side ds object (it shouldn't), we'd catch it here.
  score_after <- ts_score(tree, ds)
  expect_equal(score_before, score_after,
               info = "Dataset should not be corrupted by ratchet")
})


# --- 8. Many cycles: 50 cycles on a moderate dataset ---

test_that("Ratchet survives 50 cycles without crash or corruption", {
  set.seed(1598)
  tree <- as.phylo(1, 20)
  mat <- matrix(sample(0:2, 20 * 10, replace = TRUE),
                nrow = 20,
                dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 50L)

  expect_equal(result$n_cycles, 50L)
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))

  # Topology sanity
  expect_equal(nrow(result$edge), 2 * 19)
  tips <- sort(result$edge[result$edge[, 2] <= 20, 2])
  expect_equal(tips, 1:20)
})


# --- 9. Large tree (75 tips) ---

test_that("Ratchet handles 75-tip tree", {
  set.seed(7023)
  n <- 75
  tree <- as.phylo(1, n)
  mat <- matrix(sample(0:3, n * 20, replace = TRUE),
                nrow = n,
                dimnames = list(paste0("t", seq_len(n)), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 5L)

  expect_true(result$score > 0)
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
  expect_equal(nrow(result$edge), 2 * (n - 1))
})


# --- 10. Multi-state characters (4+ states) ---

test_that("Ratchet works with multi-state characters", {
  set.seed(5501)
  n <- 15
  # Characters with 5 states
  mat <- matrix(sample(0:4, n * 6, replace = TRUE),
                nrow = n,
                dimnames = list(paste0("t", seq_len(n)), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tree <- as.phylo(1, n)
  result <- ts_ratchet(tree, ds, nCycles = 5L)

  expect_true(result$score > 0)
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
})


# --- 11. Repeated ratchet calls on same tree object ---
# Verify no accumulated state corruption.

test_that("Repeated ratchet calls don't accumulate corruption", {
  set.seed(4437)
  tree <- as.phylo(1, 12)
  mat <- matrix(sample(0:1, 12 * 6, replace = TRUE),
                nrow = 12,
                dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  scores <- integer(5)
  for (i in 1:5) {
    result <- ts_ratchet(tree, ds, nCycles = 3L)
    result_tree <- tree
    result_tree$edge <- result$edge
    scores[i] <- ts_score(result_tree, ds)
    expect_equal(result$score, scores[i],
                 info = paste("Iteration", i))
  }
  # All runs should find valid (positive) scores
  expect_true(all(scores > 0))
})


# --- 12. Congreve-Lamsdell with many cycles ---

test_that("Ratchet on Congreve-Lamsdell: 20 cycles, score verified", {
  skip_if_not_installed("TreeSearch")
  data(congreveLamsdellMatrices, package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)

  set.seed(2811)
  tree <- as.phylo(sample.int(1e6, 1), length(dataset))
  start_score <- ts_score(tree, ds)

  result <- ts_ratchet(tree, ds, nCycles = 20L, maxHits = 3L)

  # Must improve over random start
  expect_true(result$score < start_score)

  # Score must be verified
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))

  expect_equal(result$n_cycles, 20L)
})


# --- 13. Uniform data (all tips identical) ---
# Optimal score should be 0 for any topology.

test_that("Ratchet on uniform data returns score 0", {
  tree <- as.phylo(1, 8)
  mat <- matrix(rep(0, 8 * 3),
                nrow = 8,
                dimnames = list(paste0("t", 1:8), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 2L)
  expect_equal(result$score, 0)
})


# --- 14. All-different data (each tip unique) ---

test_that("Ratchet on all-unique-tip data", {
  n <- 8
  # Identity-like matrix: each tip has a unique state pattern
  mat <- matrix(0, nrow = n, ncol = n,
                dimnames = list(paste0("t", 1:n), NULL))
  for (i in 1:n) mat[i, i] <- 1
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  tree <- as.phylo(1, n)
  result <- ts_ratchet(tree, ds, nCycles = 3L)

  # Each character needs exactly 1 step (one tip differs from rest)
  # Score = n characters * 1 step each = n
  expect_equal(result$score, n)
})


# --- 15. Topology validity under heavy perturbation ---

test_that("Topology valid after heavy perturbation (prob=0.5, 20 cycles)", {
  set.seed(9361)
  n <- 25
  tree <- as.phylo(1, n)
  mat <- matrix(sample(0:1, n * 10, replace = TRUE),
                nrow = n,
                dimnames = list(paste0("t", seq_len(n)), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  result <- ts_ratchet(tree, ds, nCycles = 20L, perturbProb = 0.5)

  # Topology checks
  edge <- result$edge
  expect_equal(nrow(edge), 2 * (n - 1))
  tips <- sort(edge[edge[, 2] <= n, 2])
  expect_equal(tips, seq_len(n))

  # All internal nodes should appear as parents
  internal <- sort(unique(edge[, 1]))
  expect_equal(internal, (n + 1):(2 * n - 1))

  # Score verified
  result_tree <- tree
  result_tree$edge <- result$edge
  expect_equal(result$score, ts_score(result_tree, ds))
})
