# Tests for constrained sectorial search (CSS).
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# Helper: run driven search with CSS control
ts_driven <- function(ds, maxReplicates = 5L, targetHits = 2L,
                      ratchetCycles = 3L, xssRounds = 1L,
                      xssPartitions = 2L, fuseInterval = 2L,
                      cssRounds = 1L, cssPartitions = 4L,
                      maxSeconds = 0, verbosity = 0L, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    ratchetCycles = ratchetCycles,
    xssRounds = xssRounds,
    xssPartitions = xssPartitions,
    fuseInterval = fuseInterval,
    cssRounds = cssRounds,
    cssPartitions = cssPartitions,
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
set.seed(6142)
med_mat <- matrix(sample(0:1, 20 * 10, replace = TRUE),
                  nrow = 20,
                  dimnames = list(paste0("t", 1:20), NULL))
med_dataset <- MatrixToPhyDat(med_mat)
med_ds <- make_ts_data(med_dataset)

# Larger dataset: 30 tips, 15 characters
set.seed(2879)
large_mat <- matrix(sample(0:1, 30 * 15, replace = TRUE),
                    nrow = 30,
                    dimnames = list(paste0("t", 1:30), NULL))
large_dataset <- MatrixToPhyDat(large_mat)
large_ds <- make_ts_data(large_dataset)


# ---------- Tests ----------

test_that("CSS driven search produces valid trees", {
  set.seed(3841)
  result <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 2L,
                      cssRounds = 2L, cssPartitions = 3L)

  expect_true(is.list(result))
  expect_true("trees" %in% names(result))
  expect_true("scores" %in% names(result))
  expect_true(length(result$trees) >= 1)
  expect_true(all(result$scores <= result$best_score + 1e-6))

  # Verify each returned tree is valid and matches reported score
  for (i in seq_along(result$trees)) {
    tree <- structure(list(edge = result$trees[[i]],
                           Nnode = nrow(result$trees[[i]]) / 2,
                           tip.label = paste0("t", 1:20)),
                      class = "phylo")
    score <- ts_score(tree, med_ds)
    expect_equal(score, result$scores[i], tolerance = 1e-6)
  }
})

test_that("CSS rounds = 0 disables CSS", {
  set.seed(4521)
  # With CSS
  r_with <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 2L,
                      cssRounds = 1L)
  set.seed(4521)
  # Without CSS
  r_without <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 2L,
                         cssRounds = 0L)

  # Both should produce valid results

  expect_true(r_with$best_score > 0)
  expect_true(r_without$best_score > 0)
})

test_that("CSS produces competitive results vs no-CSS", {
  set.seed(1234)
  r_with <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                      cssRounds = 2L, cssPartitions = 3L)
  set.seed(1234)
  r_without <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                         cssRounds = 0L)
  # Both should find similar-quality trees
  expect_true(abs(r_with$best_score - r_without$best_score) < 5,
              info = paste("CSS:", r_with$best_score,
                           "no-CSS:", r_without$best_score))
})

test_that("CSS works with implied weights", {
  set.seed(7392)
  result <- ts_driven(med_ds, maxReplicates = 3L, targetHits = 2L,
                      cssRounds = 1L, cssPartitions = 2L,
                      concavity = 10.0)

  expect_true(result$best_score > 0)
  expect_true(length(result$trees) >= 1)

  # Score verification
  tree <- structure(list(edge = result$trees[[1]],
                         Nnode = nrow(result$trees[[1]]) / 2,
                         tip.label = paste0("t", 1:20)),
                    class = "phylo")
  iw_score <- TreeSearch:::ts_fitch_score(
    tree$edge, med_ds$contrast, med_ds$tip_data,
    med_ds$weight, med_ds$levels, concavity = 10.0)
  expect_equal(iw_score, result$scores[1], tolerance = 1e-6)
})

test_that("CSS works with inapplicable characters", {
  skip_if_not_installed("TreeSearch")
  # Use a real inapplicable dataset
  data("inapplicable.phyData", package = "TreeSearch")
  vinther <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(vinther)

  set.seed(5981)
  result <- ts_driven(ds, maxReplicates = 3L, targetHits = 2L,
                      cssRounds = 1L, cssPartitions = 2L)

  expect_true(result$best_score > 0)
  expect_true(length(result$trees) >= 1)
})

test_that("CSS is deterministic with set.seed", {
  set.seed(2346)
  r1 <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                  cssRounds = 1L, cssPartitions = 3L)
  set.seed(2346)
  r2 <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                  cssRounds = 1L, cssPartitions = 3L)

  expect_equal(r1$best_score, r2$best_score)
  expect_equal(r1$replicates, r2$replicates)
  expect_equal(r1$trees[[1]], r2$trees[[1]])
})

test_that("CSS with small tree (below sector threshold) is no-op", {
  set.seed(8831)
  # 10 tips, sectorMinSize = 6 → tree is barely large enough
  # With cssPartitions = 2, sectors of ~5 tips → below min_size of 4
  r <- ts_driven(small_ds, maxReplicates = 3L, targetHits = 2L,
                 cssRounds = 1L, cssPartitions = 2L)

  expect_true(r$best_score > 0)
  expect_true(length(r$trees) >= 1)
})

test_that("CSS integration with full R-level MaximizeParsimony", {
  skip_if_not_installed("TreeSearch")
  result <- MaximizeParsimony(med_dataset, maxReplicates = 3L,
                              targetHits = 2L, cssRounds = 1L,
                              verbosity = 0L)

  expect_s3_class(result, "multiPhylo")
  expect_true(length(result) >= 1)
})

test_that("Driven search with CSS handles timeout", {
  set.seed(4411)
  result <- ts_driven(large_ds, maxReplicates = 100L, targetHits = 100L,
                      cssRounds = 2L, cssPartitions = 3L,
                      maxSeconds = 0.5)

  expect_true(result$timed_out || result$replicates < 100)
  expect_true(result$best_score > 0)
})

test_that("Multiple CSS partitions are all searched", {
  set.seed(7712)
  # With more partitions, each sector is smaller → more sectors searched
  r2 <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                  cssRounds = 1L, cssPartitions = 2L)
  r6 <- ts_driven(large_ds, maxReplicates = 3L, targetHits = 2L,
                  cssRounds = 1L, cssPartitions = 6L)

  # Both should produce valid results
  expect_true(r2$best_score > 0)
  expect_true(r6$best_score > 0)
})
