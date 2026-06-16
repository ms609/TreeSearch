# Phase 3D: Memory layout regression tests
#
# Verifies that TBR search produces identical results after postorder
# save/restore optimization (no redundant build_postorder calls).

test_that("ts_bench_tbr_phases returns correct structure", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)
  tree <- TreeTools::RandomTree(names(ds), root = TRUE)

  result <- TreeSearch:::ts_bench_tbr_phases(
    tree$edge, at$contrast, tip_data, at$weight, at$levels
  )

  expect_type(result, "list")
  expect_equal(result$n_tips, 23L)
  expect_true(result$n_node > 0)
  expect_true(result$n_blocks > 0)
  expect_true(result$total_words > 0)
  expect_true(result$score > 0)
  expect_true(result$n_clips > 0)
  expect_true(result$n_candidates > 0)

  # Timing fields exist and are non-negative
  expect_true(result$time_full_rescore_us >= 0)
  expect_true(result$time_clip_incr_us >= 0)
  expect_true(result$time_indirect_us >= 0)
  expect_true(result$time_unclip_us >= 0)
  expect_true(result$time_snapshot_save_us >= 0)
  expect_true(result$time_snapshot_restore_us >= 0)
  expect_true(result$snapshot_bytes > 0)
})

test_that("TBR search with postorder optimization gives correct scores", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  set.seed(6283)
  tree <- TreeTools::RandomTree(names(ds), root = TRUE)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, at$contrast, tip_data, at$weight, at$levels,
    maxHits = 1L, acceptEqual = FALSE
  )

  expect_true(result$score > 0)
  expect_true(result$n_evaluated > 0)

  # Score should match independent full rescore
  score_check <- TreeSearch:::ts_fitch_score(
    result$edge, at$contrast, tip_data, at$weight, at$levels
  )
  expect_equal(result$score, score_check)
})

test_that("TBR search deterministic with set.seed", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  set.seed(4419)
  tree <- TreeTools::RandomTree(names(ds), root = TRUE)
  edge <- tree$edge

  set.seed(7701)
  r1 <- TreeSearch:::ts_tbr_search(
    edge, at$contrast, tip_data, at$weight, at$levels,
    maxHits = 1L, acceptEqual = TRUE
  )

  set.seed(7701)
  r2 <- TreeSearch:::ts_tbr_search(
    edge, at$contrast, tip_data, at$weight, at$levels,
    maxHits = 1L, acceptEqual = TRUE
  )

  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
  expect_equal(r1$n_accepted, r2$n_accepted)
})

test_that("TBR search correct with NA dataset", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  set.seed(3392)
  tree <- TreeTools::RandomTree(names(ds), root = TRUE)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, at$contrast, tip_data, at$weight, at$levels,
    maxHits = 5L, acceptEqual = TRUE
  )

  # Verify score is consistent
  score_check <- TreeSearch:::ts_fitch_score(
    result$edge, at$contrast, tip_data, at$weight, at$levels
  )
  expect_equal(result$score, score_check)

  # Score should improve from random tree
  initial_score <- TreeSearch:::ts_fitch_score(
    tree$edge, at$contrast, tip_data, at$weight, at$levels
  )
  expect_true(result$score <= initial_score)
})

test_that("TBR search correct with IW scoring", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  # Per-character minimum step counts (Farris's m), one value per pattern.
  # `apply(at$contrast, 2, ...)` would iterate over contrast *columns* (states),
  # yielding a length-nStates vector that the scorer rejects (and that formerly
  # read out of bounds). Mirror MaximizeParsimony()'s own IW setup instead.
  min_steps <- as.integer(MinimumLength(ds, compress = TRUE))

  set.seed(8816)
  tree <- TreeTools::RandomTree(names(ds), root = TRUE)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, at$contrast, tip_data, at$weight, at$levels,
    maxHits = 1L, acceptEqual = FALSE,
    min_steps = min_steps, concavity = 10.0
  )

  expect_true(result$score > 0)

  # Score should match independent check
  score_check <- TreeSearch:::ts_fitch_score(
    result$edge, at$contrast, tip_data, at$weight, at$levels,
    min_steps = min_steps, concavity = 10.0
  )
  expect_equal(result$score, score_check, tolerance = 1e-8)
})

test_that("Driven search works correctly after TBR postorder optimization", {
  data("inapplicable.phyData", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  set.seed(1134)
  result <- TreeSearch:::ts_driven_search(
    at$contrast, tip_data, at$weight, at$levels,
    maxReplicates = 2L, targetHits = 2L,
    tbrMaxHits = 1L, ratchetCycles = 2L,
    xssRounds = 1L, fuseInterval = 10L,
    verbosity = 0L
  )

  expect_true(result$best_score > 0)
  expect_true(result$pool_size >= 1)

  # Verify best tree's score
  score_check <- TreeSearch:::ts_fitch_score(
    result$trees[[1]], at$contrast, tip_data, at$weight, at$levels
  )
  expect_equal(result$best_score, score_check)
})

test_that("Bench function works with synthetic binary data", {
  library(TreeTools)

  set.seed(2261)
  n <- 30
  tree <- RandomTree(n, root = TRUE)
  mat <- matrix(sample(c("0", "1"), n * 50, replace = TRUE), n, 50,
                dimnames = list(tree$tip.label, NULL))
  ds <- MatrixToPhyDat(mat)
  at <- attributes(ds)
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  result <- TreeSearch:::ts_bench_tbr_phases(
    tree$edge, at$contrast, tip_data, at$weight, at$levels
  )

  expect_equal(result$n_tips, 30L)
  expect_false(result$has_na)
  expect_true(result$n_blocks > 0)
  expect_true(result$n_candidates > 0)
})
