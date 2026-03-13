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

# Helper: run driven search
ts_driven <- function(ds, maxReplicates = 5L, targetHits = 2L,
                      ratchetCycles = 3L, xssRounds = 1L,
                      xssPartitions = 2L, fuseInterval = 2L, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    xssRounds = xssRounds,
    xssPartitions = xssPartitions,
    fuseInterval = fuseInterval,
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
  expect_true("edge" %in% names(result))
  expect_true("score" %in% names(result))
  expect_true("replicates" %in% names(result))
  expect_true("hits_to_best" %in% names(result))
  expect_true("pool_size" %in% names(result))
})

test_that("Driven search score matches independent verification", {
  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 1L)

  result_tree <- list(edge = result$edge, Nnode = nrow(result$edge) / 2L,
                      tip.label = paste0("t", 1:10))
  class(result_tree) <- "phylo"
  expected_score <- ts_score(result_tree, small_ds)
  expect_equal(result$score, expected_score)
})

test_that("Driven search converges with targetHits=1", {
  result <- ts_driven(small_ds, maxReplicates = 20L, targetHits = 1L,
                      ratchetCycles = 1L)

  expect_true(result$hits_to_best >= 1)
  expect_true(result$pool_size >= 1)
})

test_that("Driven search improves over random Wagner tree", {
  # A random Wagner tree should be suboptimal; driven search should improve
  set.seed(4291)
  wagner <- TreeSearch:::ts_random_wagner_tree(
    small_ds$contrast, small_ds$tip_data, small_ds$weight, small_ds$levels
  )

  result <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L)

  expect_true(result$score <= wagner$score)
})

test_that("Driven search works on tiny dataset", {
  result <- ts_driven(tiny_ds, maxReplicates = 2L, targetHits = 1L,
                      ratchetCycles = 1L, xssRounds = 0L)

  expect_true(result$score > 0)
  expect_true(result$pool_size >= 1)
})

test_that("Driven search works on medium dataset", {
  result <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L)

  expect_true(result$score > 0)
  expect_true(result$replicates >= 1)

  # Verify score
  result_tree <- list(edge = result$edge, Nnode = nrow(result$edge) / 2L,
                      tip.label = paste0("t", 1:20))
  class(result_tree) <- "phylo"
  expect_equal(result$score, ts_score(result_tree, med_ds))
})

test_that("Multiple replicates improve search quality", {
  # Single replicate
  r1 <- ts_driven(med_ds, maxReplicates = 1L, targetHits = 1L,
                   ratchetCycles = 1L, fuseInterval = 100L)

  # Multiple replicates with fusing
  scores <- numeric(5)
  for (i in seq_along(scores)) {
    r <- ts_driven(med_ds, maxReplicates = 5L, targetHits = 3L,
                   ratchetCycles = 3L, fuseInterval = 2L)
    scores[i] <- r$score
  }

  # Best of multiple runs should be <= single replicate (usually)
  expect_true(min(scores) <= r1$score)
})

test_that("Pool accumulates trees", {
  result <- ts_driven(small_ds, maxReplicates = 5L, targetHits = 10L,
                      ratchetCycles = 1L, fuseInterval = 100L)

  # With targetHits=10 and maxReplicates=5, pool should have some trees
  expect_true(result$pool_size >= 1)
  expect_true(result$replicates == 5L)
})

test_that("Driven search handles edge case parameters", {
  # Zero ratchet cycles
  r <- ts_driven(small_ds, maxReplicates = 2L, targetHits = 1L,
                 ratchetCycles = 0L)
  expect_true(r$score > 0)

  # Large targetHits forces all replicates
  r2 <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 100L,
                  ratchetCycles = 1L)
  expect_equal(r2$replicates, 3L)
})
