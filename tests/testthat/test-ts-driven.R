# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Tests for the C++ driven search engine.
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search
ts_driven <- function(ds, maxReplicates = 5L, targetHits = 2L,
                      ratchetCycles = 3L, xssRounds = 1L,
                      xssPartitions = 2L, fuseInterval = 2L,
                      maxSeconds = 0, verbosity = 0L, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    xssRounds = xssRounds,
    xssPartitions = xssPartitions,
    fuseInterval = fuseInterval,
    maxSeconds = maxSeconds,
    verbosity = verbosity,
    ...
  )
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

# Medium dataset: 20 tips, 10 characters
set.seed(8317)
med_mat <- matrix(sample(0:1, 20 * 10, replace = TRUE),
                  nrow = 20,
                  dimnames = list(paste0("t", 1:20), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)

# Tiny dataset: 8 tips, 3 characters
tiny_mat <- matrix(c(
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1,
  0, 1, 0, 1, 0, 1, 0, 1
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
tiny_dataset <- MatrixToPhyDat(tiny_mat)
tiny_ds <- make_ts_data(tiny_dataset)


test_that("Driven search returns valid structure", {
  result <- ts_driven(small_ds, maxReplicates = 2L, targetHits = 1L,
                      ratchetCycles = 1L)

  expect_true(is.list(result))
  expect_true("trees" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true("best_score" %in% names(result))
  expect_true("replicates" %in% names(result))
  expect_true("hits_to_best" %in% names(result))
  expect_true("pool_size" %in% names(result))
  expect_true("timed_out" %in% names(result))

  # trees is a list of edge matrices
 expect_true(is.list(result$trees))
  expect_true(length(result$trees) >= 1)
  expect_true(is.matrix(result$trees[[1]]))
  expect_equal(ncol(result$trees[[1]]), 2L)

  # scores vector matches trees length
  expect_equal(length(result$scores), length(result$trees))
})

test_that("Driven search score matches independent verification", {
  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 1L)

  # Verify best tree score
  best_edge <- result$trees[[1]]
  result_tree <- list(edge = best_edge, Nnode = nrow(best_edge) / 2L,
                      tip.label = paste0("t", 1:10))
  class(result_tree) <- "phylo"
  expected_score <- ts_score(result_tree, small_ds)
  expect_equal(result$best_score, expected_score)

  # All reported scores should be verifiable
  for (i in seq_along(result$trees)) {
    tr <- list(edge = result$trees[[i]], Nnode = nrow(result$trees[[i]]) / 2L,
               tip.label = paste0("t", 1:10))
    class(tr) <- "phylo"
    expect_equal(result$scores[i], ts_score(tr, small_ds))
  }
})

test_that("Driven search converges with targetHits=1", {
  result <- ts_driven(small_ds, maxReplicates = 20L, targetHits = 1L,
                      ratchetCycles = 1L)

  expect_true(result$hits_to_best >= 1)
  expect_true(result$pool_size >= 1)
})

test_that("Driven search improves over random Wagner tree", {
  set.seed(4291)
  wagner <- TreeSearch:::ts_random_wagner_tree(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    concavity = Inf
  )

  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L)

  expect_true(result$best_score <= wagner$score)
})

test_that("Driven search works on tiny dataset", {
  result <- ts_driven(tiny_ds, maxReplicates = 2L, targetHits = 1L,
                      ratchetCycles = 1L, xssRounds = 0L)

  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1)
})

test_that("Driven search works on medium dataset", {
  result <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L)

  expect_true(result$best_score > 0)
  expect_true(result$replicates >= 1)

  # Verify best tree score
  best_edge <- result$trees[[1]]
  result_tree <- list(edge = best_edge, Nnode = nrow(best_edge) / 2L,
                      tip.label = paste0("t", 1:20))
  class(result_tree) <- "phylo"
  expect_equal(result$best_score, ts_score(result_tree, med_ds))
})

test_that("Multiple replicates improve search quality", {
  r1 <- ts_driven(med_ds, maxReplicates = 1L, targetHits = 1L,
                   ratchetCycles = 1L, fuseInterval = 100L)

  scores <- numeric(2)
  for (i in seq_along(scores)) {
    r <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 2L,
                   ratchetCycles = 2L, fuseInterval = 2L)
    scores[i] <- r$best_score
  }

  expect_true(min(scores) <= r1$best_score)
})

test_that("Pool accumulates trees", {
  result <- ts_driven(small_ds, maxReplicates = 5L, targetHits = 10L,
                      ratchetCycles = 1L, fuseInterval = 100L)

  expect_true(result$pool_size >= 1)
  expect_true(result$replicates == 5L)
  # Pool returns all trees
  expect_equal(length(result$trees), result$pool_size)
})

test_that("Driven search handles edge case parameters", {
  # Zero ratchet cycles
  r <- ts_driven(small_ds, maxReplicates = 2L, targetHits = 1L,
                 ratchetCycles = 0L)
  expect_true(r$best_score > 0)

  # Large targetHits forces all replicates
  r2 <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 100L,
                  ratchetCycles = 1L)
  expect_equal(r2$replicates, 3L)
})

# ---------- New feature tests (Agent C) ----------

test_that("All pool trees are returned", {
  # Force enough replicates to accumulate pool entries
  result <- ts_driven(small_ds, maxReplicates = 5L, targetHits = 100L,
                      ratchetCycles = 1L, fuseInterval = 100L,
                      poolSuboptimal = 0.0)

  expect_equal(length(result$trees), result$pool_size)
  expect_equal(length(result$scores), result$pool_size)
  # All scores should be the best (suboptimal = 0)
  expect_true(all(result$scores == result$best_score))
})

test_that("Suboptimal tree collection works", {
  # Allow suboptimal trees within 2 steps
  result <- ts_driven(med_ds, maxReplicates = 10L, targetHits = 100L,
                      ratchetCycles = 2L, fuseInterval = 100L,
                      poolSuboptimal = 2.0)

  # With suboptimal > 0, we may have trees at different scores
  expect_true(length(result$trees) >= 1)
  # All scores should be within tolerance of best
  expect_true(all(result$scores <= result$best_score + 2.0 + 1e-9))
})

test_that("Timeout stops search early", {
  # Set a very short timeout
  result <- ts_driven(med_ds, maxReplicates = 1000L, targetHits = 1000L,
                      ratchetCycles = 5L, maxSeconds = 0.5,
                      perturbStopFactor = 0L)

  # Should not have completed all 1000 replicates
  expect_true(result$replicates < 1000L)
  expect_true(result$timed_out)

  # Should still have valid results
  if (result$pool_size > 0) {
    expect_true(result$best_score > 0)
    expect_equal(length(result$trees), result$pool_size)
  }
})

test_that("Timeout of 0 means no timeout", {
  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 1L, maxSeconds = 0)

  expect_false(result$timed_out)
})

test_that("Verbosity does not break search", {
  # verbosity=1 and verbosity=2 should work without error
  expect_no_error({
    r1 <- ts_driven(small_ds, maxReplicates = 2L, targetHits = 1L,
                    ratchetCycles = 1L, verbosity = 1L)
  })
  expect_true(r1$best_score > 0)

  expect_no_error({
    r2 <- ts_driven(tiny_ds, maxReplicates = 2L, targetHits = 1L,
                    ratchetCycles = 1L, xssRounds = 0L, verbosity = 2L)
  })
  expect_true(r2$best_score > 0)
})

test_that("Zero replicates returns empty result", {
  result <- ts_driven(small_ds, maxReplicates = 0L)

  expect_equal(length(result$trees), 0)
  expect_equal(length(result$scores), 0)
  expect_equal(result$pool_size, 0)
  expect_false(result$timed_out)
})

test_that("MaximizeParsimony() uses C++ engine", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  result <- MaximizeParsimony(dataset, maxReplicates = 2L, targetHits = 1L,
                              verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(!is.null(attr(result, "score")))
  expect_true(attr(result, "score") > 0)
})

test_that("MaximizeParsimony() supports IW natively", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  result <- MaximizeParsimony(dataset, concavity = 10,
                              maxReplicates = 2L, targetHits = 1L,
                              verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(attr(result, "score") > 0)
})

test_that("Ratchet taper runs without error and finds valid score", {
  result <- ts_driven(med_ds, maxReplicates = 6L, targetHits = 4L,
                      ratchetCycles = 3L, ratchetTaper = TRUE)
  expect_true(is.list(result))
  expect_true(result$best_score > 0)
  expect_true(result$replicates >= 1L)
})

test_that("Ratchet taper produces comparable scores to non-taper", {
  set.seed(6183)
  no_taper <- ts_driven(small_ds, maxReplicates = 8L, targetHits = 4L,
                        ratchetCycles = 3L, ratchetTaper = FALSE)
  set.seed(6183)
  with_taper <- ts_driven(small_ds, maxReplicates = 8L, targetHits = 4L,
                          ratchetCycles = 3L, ratchetTaper = TRUE)
  # Taper should not make things dramatically worse (within 2 steps)
  expect_lte(with_taper$best_score, no_taper$best_score + 2)
})

test_that("Ratchet taper works with adaptive level", {
  result <- ts_driven(med_ds, maxReplicates = 6L, targetHits = 4L,
                      ratchetCycles = 3L, ratchetTaper = TRUE,
                      adaptiveLevel = TRUE)
  expect_true(result$best_score > 0)
})

test_that("SearchControl includes ratchetTaper", {
  ctrl <- SearchControl(ratchetTaper = TRUE)
  expect_true(ctrl$ratchetTaper)
  ctrl2 <- SearchControl()
  expect_false(ctrl2$ratchetTaper)
})

test_that("perturbStopFactor stops search after unsuccessful replicates", {
  # Small dataset with 10 tips: perturbStopFactor=1 -> limit = 10 replicates.
  result <- ts_driven(small_ds, maxReplicates = 100L, targetHits = 100L,
                      ratchetCycles = 1L, xssRounds = 0L,
                      perturbStopFactor = 1L)
  expect_lt(result$replicates, 100L)
  expect_true(result$pool_size >= 1)
  expect_true(result$best_score > 0)
  expect_true(result$perturb_stop)
})

test_that("perturbStopFactor=0 disables the rule", {
  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 1L, perturbStopFactor = 0L)
  expect_true(result$pool_size >= 1)
})

test_that("MaximizeParsimony2() is deprecated alias", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  expect_warning(
    MaximizeParsimony2(dataset, maxReplicates = 2L, targetHits = 1L,
                       verbosity = 0L),
    "deprecated"
  )
})
