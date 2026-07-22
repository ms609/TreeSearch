# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# ---------- Test datasets ----------

# 10 tips, informative characters
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

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

# 8 tips, 3 characters (tiny)
tiny_mat <- matrix(c(
  0, 0, 0, 0, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1,
  0, 1, 0, 1, 0, 1, 0, 1
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
tiny_dataset <- MatrixToPhyDat(tiny_mat)
tiny_ds <- make_ts_data(tiny_dataset)

# ---------- Jackknife tests ----------

test_that("Jackknife returns valid tree", {
  set.seed(7284)
  result <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = FALSE, jackProportion = 2 / 3,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true(is.matrix(result$edge))
  expect_equal(ncol(result$edge), 2L)
  # 10 tips → 18 edges
  expect_equal(nrow(result$edge), 18L)
  expect_true(result$score > 0)
})

test_that("Jackknife produces different trees across runs", {
  edges <- list()
  for (i in 1:5) {
    r <- TreeSearch:::ts_resample_search(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      bootstrap = FALSE, jackProportion = 2 / 3,
      maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
    )
    edges[[i]] <- r$edge
  }
  # At least 2 distinct topologies across 5 runs (almost certain)
  n_unique <- length(unique(lapply(edges, function(e) sort(paste(e[, 1], e[, 2])))))
  expect_true(n_unique >= 2)
})

test_that("Jackknife proportion parameter works", {
  # Very low proportion should still return a valid tree
  set.seed(3916)
  r <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = FALSE, jackProportion = 0.3,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(is.matrix(r$edge))
  expect_true(r$score > 0)
})

# ---------- Bootstrap tests ----------

test_that("Bootstrap returns valid tree", {
  set.seed(5193)
  result <- TreeSearch:::ts_resample_search(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    bootstrap = TRUE,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )

  expect_true(is.matrix(result$edge))
  expect_equal(ncol(result$edge), 2L)
  expect_equal(nrow(result$edge), 18L)
  expect_true(result$score > 0)
})

test_that("Bootstrap produces different trees across runs", {
  edges <- list()
  for (i in 1:5) {
    r <- TreeSearch:::ts_resample_search(
      small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
      bootstrap = TRUE,
      maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
    )
    edges[[i]] <- r$edge
  }
  n_unique <- length(unique(lapply(edges, function(e) sort(paste(e[, 1], e[, 2])))))
  expect_true(n_unique >= 2)
})

test_that("Resample works on tiny dataset", {
  set.seed(6102)
  r_jack <- TreeSearch:::ts_resample_search(
    tiny_ds$contrast, tiny_ds$tip_data, tiny_ds$weight, tiny_ds$levels,
    bootstrap = FALSE, jackProportion = 0.5,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(r_jack$score > 0)

  r_boot <- TreeSearch:::ts_resample_search(
    tiny_ds$contrast, tiny_ds$tip_data, tiny_ds$weight, tiny_ds$levels,
    bootstrap = TRUE,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(r_boot$score > 0)
})

# ---------- Successive Approximations tests ----------

test_that("Successive approximations returns valid structure", {
  set.seed(4517)
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 5L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )

  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("sa_iterations" %in% names(result))
  expect_true("converged" %in% names(result))
  expect_true(is.matrix(result$edge))
  expect_equal(nrow(result$edge), 18L)
  expect_true(result$score > 0)
  expect_true(result$sa_iterations >= 1)
  expect_true(is.logical(result$converged))
})

test_that("SA converges or completes max iterations", {
  set.seed(2758)
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 10L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L
  )

  # Either converged before max, or completed max iterations
  expect_true(result$converged || result$sa_iterations == 10L)
})

test_that("SA with k=1 and k=5 both work", {
  set.seed(8403)
  r1 <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 1.0, maxSAIter = 3L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(r1$score > 0)

  r5 <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 5.0, maxSAIter = 3L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(r5$score > 0)
})

test_that("SA works on tiny dataset", {
  set.seed(1647)
  result <- TreeSearch:::ts_successive_approx(
    tiny_ds$contrast, tiny_ds$tip_data, tiny_ds$weight, tiny_ds$levels,
    saK = 3.0, maxSAIter = 5L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(result$score > 0)
  expect_true(result$sa_iterations >= 1)
})

test_that("SA with maxSAIter=1 completes one iteration", {
  set.seed(9371)
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 1L,
    maxReplicates = 2L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_equal(result$sa_iterations, 1L)
  expect_true(result$score > 0)
})

test_that("SA score is a valid EW parsimony score", {
  set.seed(5539)
  result <- TreeSearch:::ts_successive_approx(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels,
    saK = 3.0, maxSAIter = 5L,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 2L
  )

  # Verify by scoring the returned tree independently
  tr <- list(edge = result$edge, Nnode = nrow(result$edge) / 2L,
             tip.label = paste0("t", 1:10))
  class(tr) <- "phylo"
  expected <- TreeSearch:::ts_fitch_score(
    tr$edge, small_ds$contrast, small_ds$tip_data,
    small_ds$weight, small_ds$levels, concavity = Inf
  )
  expect_equal(result$score, expected)
})
