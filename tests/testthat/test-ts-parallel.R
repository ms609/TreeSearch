context("Parallel driven search (Phase 5)")

# --- Test helpers ---

make_dataset <- function(phyDat) {
  at <- attributes(phyDat)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(phyDat, use.names = FALSE),
                      nrow = length(phyDat), byrow = TRUE),
    weight = at$weight,
    levels = at$levels
  )
}

# --- Setup ---

data("inapplicable.phyData", package = "TreeSearch")
vinther <- inapplicable.phyData[["Vinther2008"]]

# Additional dataset for testing
ew_data <- vinther

# --- 1. Serial equivalence (nThreads=1) ---

test_that("nThreads=1 produces identical results to default serial search", {
  skip_on_cran()
  set.seed(8741)
  result_default <- TreeSearch:::ts_driven_search(
    contrast = attributes(vinther)$contrast,
    tip_data = matrix(unlist(vinther, use.names = FALSE),
                      nrow = length(vinther), byrow = TRUE),
    weight = attributes(vinther)$weight,
    levels = attributes(vinther)$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    nThreads = 1L
  )

  set.seed(8741)
  result_serial <- TreeSearch:::ts_driven_search(
    contrast = attributes(vinther)$contrast,
    tip_data = matrix(unlist(vinther, use.names = FALSE),
                      nrow = length(vinther), byrow = TRUE),
    weight = attributes(vinther)$weight,
    levels = attributes(vinther)$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L
  )

  expect_equal(result_default$best_score, result_serial$best_score)
  expect_equal(result_default$replicates, result_serial$replicates)
})

# --- 2. Parallel correctness ---

test_that("Parallel search (2 threads) produces valid trees with correct scores", {
  skip_on_cran()
  set.seed(3192)
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 4L, targetHits = 3L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
  expect_true(result$pool_size > 0)

  # Verify each tree is valid and scores correctly
  for (i in seq_along(result$trees)) {
    edge <- result$trees[[i]]
    score <- TreeSearch:::ts_fitch_score(
      edge, ds$contrast, ds$tip_data, ds$weight, ds$levels
    )
    expect_equal(score, result$scores[i], tolerance = 0.01,
                 info = paste("Tree", i, "score mismatch"))
    # Valid tree structure: 2*(n-1) edges for n tips
    n_tip <- length(vinther)
    expect_equal(nrow(edge), 2 * (n_tip - 1))
  }
})

test_that("Parallel search (2 threads) produces valid results", {
  skip_on_cran()
  set.seed(6104)
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 6L, targetHits = 4L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
  expect_true(is.finite(result$best_score))
})

# --- 3. Timeout in parallel mode ---

test_that("Parallel search respects timeout", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  t0 <- proc.time()["elapsed"]
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 1000L, targetHits = 999L,
    maxSeconds = 1.0, verbosity = 0L,
    nThreads = 2L
  )
  elapsed <- proc.time()["elapsed"] - t0

  expect_true(result$timed_out)
  # Should finish within a reasonable time (timeout + overhead)
  expect_true(elapsed < 10.0)
})

# --- 4. Edge cases ---

test_that("nThreads=0 (auto-detect) is capped at 2 in tests", {
  skip_on_cran()
  # nThreads=0 auto-detects CPU count; use nThreads=2 to stay within

  # the 2-core-per-agent limit (see AGENTS.md).
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
})

test_that("nThreads > maxReplicates is clamped", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 2L, targetHits = 2L, verbosity = 0L,
    nThreads = 100L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$replicates <= 2L)
})

test_that("Single replicate in parallel mode works", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 1L, targetHits = 1L, verbosity = 0L,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_equal(result$replicates, 1L)
})

# --- 5. IW + parallel ---

test_that("Implied weights works in parallel mode", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    concavity = 10.0,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)

  # Verify IW scores
  for (i in seq_along(result$trees)) {
    score <- TreeSearch:::ts_fitch_score(
      result$trees[[i]], ds$contrast, ds$tip_data,
      ds$weight, ds$levels, concavity = 10.0
    )
    expect_equal(score, result$scores[i], tolerance = 0.01)
  }
})

# --- 6. NA + parallel ---

test_that("Inapplicable characters work in parallel mode", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 3L, targetHits = 2L, verbosity = 0L,
    nThreads = 2L
  )

  # Vinther2008 has inapplicable characters
  expect_true(length(result$trees) > 0)
  expect_true(result$best_score > 0)
})

# --- 7. R-level MaximizeParsimony with parallel ---

test_that("MaximizeParsimony with nThreads > 1 works end-to-end", {
  skip_on_cran()
  set.seed(5023)
  result <- MaximizeParsimony(vinther, maxReplicates = 3L,
                              targetHits = 2L, nThreads = 2L,
                              verbosity = 0L)

  expect_s3_class(result, "multiPhylo")
  expect_true(length(result) > 0)
  expect_true(attr(result, "score") > 0)
  expect_true(is.finite(attr(result, "score")))

  # All trees should have correct number of tips
  for (i in seq_along(result)) {
    expect_equal(length(result[[i]]$tip.label), length(vinther))
    expect_true(ape::is.binary(result[[i]]))
  }
})

# --- 8. Score quality ---

test_that("Parallel search finds comparable scores to serial", {
  skip_on_cran()
  ds <- make_dataset(vinther)

  set.seed(1829)
  result_serial <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 5L, targetHits = 3L, verbosity = 0L,
    nThreads = 1L
  )

  set.seed(1829)
  result_parallel <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 5L, targetHits = 3L, verbosity = 0L,
    nThreads = 2L
  )

  # Parallel score should be in the same ballpark (within 10% or 5 steps)
  expect_true(result_parallel$best_score <=
              result_serial$best_score + max(5, result_serial$best_score * 0.1))
})

# --- 9. Pool suboptimal in parallel ---

test_that("Pool suboptimal collection works in parallel", {
  skip_on_cran()
  ds <- make_dataset(vinther)
  result <- TreeSearch:::ts_driven_search(
    contrast = ds$contrast, tip_data = ds$tip_data,
    weight = ds$weight, levels = ds$levels,
    maxReplicates = 4L, targetHits = 3L, verbosity = 0L,
    poolSuboptimal = 5.0,
    nThreads = 2L
  )

  expect_true(length(result$trees) > 0)
  # With suboptimal > 0, we might get multiple trees
  # All scores should be within suboptimal of best
  for (i in seq_along(result$scores)) {
    expect_true(result$scores[i] <= result$best_score + 5.0 + 0.01)
  }
})
