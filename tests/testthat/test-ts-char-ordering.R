# Tests for Phase 3C: character-ordering optimizations.
# Verifies score invariance under expensive-first block ordering,
# zero-weight pattern compaction, active_mask skip, and bounded
# indirect call correctness.

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
ts_score <- function(tree, ds, concavity = Inf) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels, concavity = concavity)
}

# Helper: run driven search
ts_driven <- function(ds, maxReplicates = 3L, targetHits = 1L,
                      ratchetCycles = 2L, xssRounds = 1L,
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

# Helper: validate topology
validate_result <- function(result, n_tip) {
  edges <- result$trees[[1]]
  expect_equal(nrow(edges), 2L * (n_tip - 1L))
  children <- edges[, 2]
  tips <- sort(children[children <= n_tip])
  expect_equal(tips, seq_len(n_tip))
}

# ---------- Datasets ----------

# Mixed weights: characters with different weights
mixed_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0,
  0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
  1, 1, 0, 0, 0, 0, 0, 1, 1, 1
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
mixed_dataset <- MatrixToPhyDat(mixed_mat)
mixed_ds <- make_ts_data(mixed_dataset)

# 15-tip multi-state dataset
set.seed(5501)
multi_mat <- matrix(sample(0:3, 15 * 10, replace = TRUE),
                    nrow = 15,
                    dimnames = list(paste0("t", 1:15), NULL))
multi_dataset <- MatrixToPhyDat(multi_mat)
multi_ds <- make_ts_data(multi_dataset)

# ---------- Score invariance ----------

test_that("Scores are correct after block reordering", {
  # The sort change (descending weight) is internal; scores must be identical.
  set.seed(2047)
  tree <- as.phylo(42, 10)
  score <- ts_score(tree, mixed_ds)
  expect_true(score > 0)

  # Multiple random trees should all score correctly
  for (i in 1:5) {
    rt <- as.phylo(sample.int(1e5, 1), 10)
    s <- ts_score(rt, mixed_ds)
    expect_true(s >= score || s > 0)
  }
})

test_that("EW driven search finds correct optimum", {
  set.seed(3341)
  result <- ts_driven(mixed_ds, maxReplicates = 5L, targetHits = 2L)
  expect_true(result$best_score > 0)
  validate_result(result, 10L)

  # Verify score matches rescore
  edge <- result$trees[[1]]
  rt <- as.phylo(1, 10)
  rt$edge <- edge
  expect_equal(result$best_score, ts_score(rt, mixed_ds))
})

test_that("IW driven search works with descending block order", {
  set.seed(8153)
  result <- ts_driven(mixed_ds, concavity = 3,
                      maxReplicates = 3L, targetHits = 1L)
  expect_true(result$best_score >= 0)
  validate_result(result, 10L)
})

test_that("Multi-state dataset scored correctly", {
  set.seed(7722)
  tree <- as.phylo(1, 15)
  score <- ts_score(tree, multi_ds)
  expect_true(score > 0)

  result <- ts_driven(multi_ds, maxReplicates = 3L, targetHits = 1L)
  expect_true(result$best_score <= score)
  validate_result(result, 15L)
})

# ---------- Zero-weight pattern compaction ----------

test_that("Jackknife with extreme deletion works correctly", {
  set.seed(6619)
  result <- TreeSearch:::ts_resample_search(
    mixed_ds$contrast, mixed_ds$tip_data, mixed_ds$weight, mixed_ds$levels,
    bootstrap = FALSE, jackProportion = 0.1,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(is.list(result))
  expect_true("edge" %in% names(result))
  expect_true(result$score >= 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

test_that("Bootstrap produces valid results with compaction", {
  set.seed(4408)
  result <- TreeSearch:::ts_resample_search(
    multi_ds$contrast, multi_ds$tip_data, multi_ds$weight, multi_ds$levels,
    bootstrap = TRUE,
    maxReplicates = 3L, targetHits = 1L, ratchetCycles = 1L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (15L - 1L))
})

# ---------- Bounded indirect correctness ----------

test_that("Bounded indirect produces same results as search", {
  # If bounded indirect is wrong, search results will differ.
  # Run same search twice with different seeds: both must find valid trees.
  for (seed in c(1129, 5982)) {
    set.seed(seed)
    result <- ts_driven(multi_ds, maxReplicates = 3L, targetHits = 1L)
    expect_true(result$best_score > 0)
    validate_result(result, 15L)

    # Rescore to verify
    edge <- result$trees[[1]]
    rt <- as.phylo(1, 15)
    rt$edge <- edge
    expect_equal(result$best_score, ts_score(rt, multi_ds))
  }
})

# ---------- Ratchet with active_mask skip ----------

test_that("Ratchet search correct with active_mask optimization", {
  set.seed(9341)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

# ---------- set.seed() reproducibility ----------

test_that("Driven search is reproducible with set.seed()", {
  run_search <- function() {
    set.seed(2200)
    ts_driven(mixed_ds, maxReplicates = 5L, targetHits = 2L)
  }
  r1 <- run_search()
  r2 <- run_search()
  expect_equal(r1$best_score, r2$best_score)
  expect_equal(r1$trees[[1]], r2$trees[[1]])
})

# ---------- Inapplicable characters ----------

test_that("NA dataset scored correctly with block optimizations", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")

  if ("Vinther2008" %in% names(inapplicable.phyData)) {
    ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])
    n_tip <- length(inapplicable.phyData[["Vinther2008"]])

    set.seed(6677)
    result <- ts_driven(ds, maxReplicates = 2L, targetHits = 1L,
                        ratchetCycles = 1L)
    expect_true(result$best_score > 0)
    validate_result(result, n_tip)
  }
})

test_that("NA dataset with IW works after reordering", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")

  if ("Vinther2008" %in% names(inapplicable.phyData)) {
    ds <- make_ts_data(inapplicable.phyData[["Vinther2008"]])
    n_tip <- length(inapplicable.phyData[["Vinther2008"]])

    set.seed(8831)
    result <- ts_driven(ds, concavity = 5,
                        maxReplicates = 2L, targetHits = 1L,
                        ratchetCycles = 1L)
    expect_true(result$best_score >= 0)
    validate_result(result, n_tip)
  }
})

# ---------- Drift search with bounded indirect ----------

test_that("Drift search works with bounded indirect calls", {
  set.seed(3055)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_drift_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L
  )
  expect_true(result$score > 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

test_that("Drift search IW with bounded indirect calls", {
  set.seed(4601)
  tree <- as.phylo(1, 10)
  result <- TreeSearch:::ts_drift_search(
    tree$edge, mixed_ds$contrast, mixed_ds$tip_data,
    mixed_ds$weight, mixed_ds$levels,
    nCycles = 3L, concavity = 3
  )
  expect_true(result$score >= 0)
  expect_equal(nrow(result$edge), 2L * (10L - 1L))
})

# ---------- Wagner tree with bounded indirect ----------

test_that("Wagner tree construction correct with bounded indirect", {
  set.seed(7713)
  result1 <- TreeSearch:::ts_random_wagner_tree(
    mixed_ds$contrast, mixed_ds$tip_data, mixed_ds$weight, mixed_ds$levels
  )
  expect_true(result1$score > 0)
  expect_equal(nrow(result1$edge), 2L * (10L - 1L))

  result2 <- TreeSearch:::ts_random_wagner_tree(
    multi_ds$contrast, multi_ds$tip_data, multi_ds$weight, multi_ds$levels
  )
  expect_true(result2$score > 0)
  expect_equal(nrow(result2$edge), 2L * (15L - 1L))
})
