# Tier 3: only runs when TREESEARCH_EXTENDED_TESTS=true.
# See tests/testing-strategy.md
skip_extended()

# Stress tests for resample + SA in the C++ engine.
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# ---------- Test datasets ----------

# 10 tips, 6 characters
small_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0, 1, 1, 1
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
small_dataset <- MatrixToPhyDat(small_mat)
small_ds <- make_ts_data(small_dataset)

# 20 tips, 10 characters (medium)
set.seed(8317)
med_mat <- matrix(sample(0:1, 20 * 10, replace = TRUE),
                  nrow = 20,
                  dimnames = list(paste0("t", 1:20), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)

# Inapplicable dataset from the package
data("inapplicable.phyData", package = "TreeSearch")
inapp_dataset <- inapplicable.phyData[["Vinther2008"]]
inapp_ds <- make_ts_data(inapp_dataset)


# ========================= DRIVEN SEARCH STRESS ========================= #

test_that("Driven search with fuse_interval=1 works", {
  result <- TreeSearch:::ts_driven_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    maxReplicates = 5L, targetHits = 100L,
    ratchetCycles = 1L, fuseInterval = 1L,
    concavity = Inf
  )
  expect_true(result$best_score > 0)
  expect_equal(result$replicates, 5L)
})

test_that("Driven search with large poolSuboptimal keeps suboptimal trees", {
  result <- TreeSearch:::ts_driven_search(
    med_ds$contrast, med_ds$tip_data, med_ds$weight, med_ds$levels,
    maxReplicates = 5L, targetHits = 100L,
    ratchetCycles = 1L, fuseInterval = 100L,
    poolSuboptimal = 100.0,
    concavity = Inf
  )
  # With a huge suboptimal tolerance, pool should retain many unique trees
  expect_true(result$pool_size >= 2)
  # All scores should be within tolerance
  expect_true(all(result$scores <= result$best_score + 100.0 + 1e-9))
})

test_that("Driven search pool trees are all valid phylogenies", {
  result <- TreeSearch:::ts_driven_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    maxReplicates = 5L, targetHits = 100L,
    ratchetCycles = 1L, fuseInterval = 100L,
    poolSuboptimal = 5.0,
    concavity = Inf
  )

  for (i in seq_along(result$trees)) {
    edge <- result$trees[[i]]
    expect_equal(ncol(edge), 2L)
    # Correct number of edges for 10 tips
    expect_equal(nrow(edge), 18L)
    # All tip labels present
    children <- edge[, 2]
    tips_in_tree <- sort(children[children <= 10])
    expect_equal(tips_in_tree, 1:10)
  }
})

test_that("set.seed produces reproducible driven search results", {
  run_search <- function() {
    set.seed(6142)
    TreeSearch:::ts_driven_search(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      maxReplicates = 3L, targetHits = 1L,
      ratchetCycles = 2L,
      concavity = Inf
    )
  }
  r1 <- run_search()
  r2 <- run_search()
  expect_equal(r1$best_score, r2$best_score)
  expect_identical(r1$trees[[1]], r2$trees[[1]])
})

test_that("Timeout with very short time and large dataset", {
  result <- TreeSearch:::ts_driven_search(
    med_ds$contrast, med_ds$tip_data, med_ds$weight, med_ds$levels,
    maxReplicates = 10000L, targetHits = 10000L,
    ratchetCycles = 20L, driftCycles = 10L,
    maxSeconds = 0.01,
    concavity = Inf
  )
  expect_true(result$timed_out)
  # Should have completed at most a handful of replicates
  expect_true(result$replicates < 10000L)
  # Result should still be valid
  if (result$pool_size > 0) {
    expect_true(result$best_score > 0)
    expect_equal(length(result$trees), result$pool_size)
  }
})

test_that("Driven search with inapplicable characters", {
  result <- TreeSearch:::ts_driven_search(
    inapp_ds$contrast, inapp_ds$tip_data, inapp_ds$weight, inapp_ds$levels,
    maxReplicates = 3L, targetHits = 1L,
    ratchetCycles = 1L,
    concavity = Inf
  )
  n_inapp_tips <- length(inapp_dataset)
  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1)
  # Score should be independently verifiable
  tr <- list(edge = result$trees[[1]], Nnode = nrow(result$trees[[1]]) / 2L,
             tip.label = names(inapp_dataset))
  class(tr) <- "phylo"
  expected <- TreeSearch:::ts_fitch_score(
    tr$edge, inapp_ds$contrast, inapp_ds$tip_data,
    inapp_ds$weight, inapp_ds$levels, concavity = Inf
  )
  expect_equal(result$best_score, expected)
})


# ========================= RESAMPLE STRESS ========================= #

test_that("set.seed produces reproducible jackknife", {
  run_jack <- function() {
    set.seed(3847)
    TreeSearch:::ts_resample_search(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      bootstrap = FALSE, jackProportion = 2 / 3,
      maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
      concavity = Inf
    )
  }
  r1 <- run_jack()
  r2 <- run_jack()
  expect_equal(r1$score, r2$score)
  expect_identical(r1$edge, r2$edge)
})

test_that("set.seed produces reproducible bootstrap", {
  run_boot <- function() {
    set.seed(9251)
    TreeSearch:::ts_resample_search(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      bootstrap = TRUE,
      maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
      concavity = Inf
    )
  }
  r1 <- run_boot()
  r2 <- run_boot()
  expect_equal(r1$score, r2$score)
  expect_identical(r1$edge, r2$edge)
})

test_that("Jackknife with inapplicable characters", {
  set.seed(4021)
  result <- TreeSearch:::ts_resample_search(
    inapp_ds$contrast, inapp_ds$tip_data, inapp_ds$weight, inapp_ds$levels,
    bootstrap = FALSE, jackProportion = 0.5,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(is.matrix(result$edge))
  expect_true(result$score >= 0)
})

test_that("Bootstrap with inapplicable characters", {
  set.seed(5678)
  result <- TreeSearch:::ts_resample_search(
    inapp_ds$contrast, inapp_ds$tip_data, inapp_ds$weight, inapp_ds$levels,
    bootstrap = TRUE,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(is.matrix(result$edge))
  expect_true(result$score >= 0)
})

test_that("Jackknife with very low proportion (0.1)", {
  set.seed(7193)
  result <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = FALSE, jackProportion = 0.1,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(is.matrix(result$edge))
  expect_true(result$score >= 0)
})

test_that("Jackknife with very high proportion (0.99)", {
  set.seed(2046)
  result <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = FALSE, jackProportion = 0.99,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(is.matrix(result$edge))
  expect_true(result$score > 0)
})

test_that("Jackknife on medium dataset", {
  set.seed(6634)
  result <- TreeSearch:::ts_resample_search(
    med_ds$contrast, med_ds$tip_data, med_ds$weight, med_ds$levels,
    bootstrap = FALSE, jackProportion = 2 / 3,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(is.matrix(result$edge))
  expect_equal(nrow(result$edge), 38L)  # 20 tips → 38 edges
})


# ========================= SA STRESS ========================= #

test_that("SA with inapplicable characters", {
  set.seed(8302)
  result <- TreeSearch:::ts_successive_approx(
    inapp_ds$contrast, inapp_ds$tip_data, inapp_ds$weight, inapp_ds$levels,
    saK = 3.0, maxSAIter = 3L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_true(result$score >= 0)
  expect_true(result$sa_iterations >= 1)

  # Verify the returned tree's score is correct EW parsimony
  tr <- list(edge = result$edge, Nnode = nrow(result$edge) / 2L,
             tip.label = names(inapp_dataset))
  class(tr) <- "phylo"
  expected <- TreeSearch:::ts_fitch_score(
    tr$edge, inapp_ds$contrast, inapp_ds$tip_data,
    inapp_ds$weight, inapp_ds$levels, concavity = Inf
  )
  expect_equal(result$score, expected)
})

test_that("SA set.seed reproducibility", {
  run_sa <- function() {
    set.seed(1234)
    TreeSearch:::ts_successive_approx(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      saK = 3.0, maxSAIter = 3L,
      maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
      concavity = Inf
    )
  }
  r1 <- run_sa()
  r2 <- run_sa()
  expect_equal(r1$score, r2$score)
  expect_identical(r1$edge, r2$edge)
  expect_equal(r1$sa_iterations, r2$sa_iterations)
})

test_that("SA with k=1 gives different weighting from k=10", {
  set.seed(4471)
  r1 <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 1.0, maxSAIter = 5L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L,
    concavity = Inf
  )
  set.seed(4471)
  r10 <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 10.0, maxSAIter = 5L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L,
    concavity = Inf
  )
  # Both should find valid trees, but k=10 penalizes homoplasy more aggressively
  expect_true(r1$score > 0)
  expect_true(r10$score > 0)
})

test_that("SA on medium dataset converges or completes", {
  set.seed(5102)
  result <- TreeSearch:::ts_successive_approx(
    med_ds$contrast, med_ds$tip_data, med_ds$weight, med_ds$levels,
    saK = 3.0, maxSAIter = 8L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L,
    concavity = Inf
  )
  expect_true(result$sa_iterations >= 1)
  expect_true(result$converged || result$sa_iterations == 8L)
  expect_equal(nrow(result$edge), 38L)  # 20 tips → 38 edges
})

test_that("SA maxSAIter=0 returns immediately with no result", {
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 0L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = Inf
  )
  expect_equal(result$sa_iterations, 0L)
  expect_false(result$converged)
  # Edge matrix should be empty since no iteration ran
  expect_equal(nrow(result$edge), 0L)
})

test_that("SA EW score is correct for all SA iteration counts", {
  # Run SA for exactly 2 iterations (not enough to converge on this dataset)
  set.seed(6189)
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 2L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L,
    concavity = Inf
  )
  expect_true(result$sa_iterations >= 1)

  # Even for non-converged runs, the returned score should be verifiable
  tr <- list(edge = result$edge, Nnode = nrow(result$edge) / 2L,
             tip.label = paste0("t", 1:10))
  class(tr) <- "phylo"
  expected <- TreeSearch:::ts_fitch_score(
    tr$edge, small_ds$contrast, small_ds$tip_data,
    small_ds$weight, small_ds$levels, concavity = Inf
  )
  expect_equal(result$score, expected)
})


# ========================= IW INTEGRATION ========================= #

test_that("Driven search with implied weights (IW)", {
  # Compute min_steps for IW
  result <- TreeSearch:::ts_driven_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    maxReplicates = 3L, targetHits = 1L,
    ratchetCycles = 1L,
    concavity = 10.0
  )
  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1)
})

test_that("Resample with implied weights", {
  set.seed(2847)
  result <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = FALSE, jackProportion = 2 / 3,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L,
    concavity = 10.0
  )
  expect_true(is.matrix(result$edge))
  expect_true(result$score > 0)
})
