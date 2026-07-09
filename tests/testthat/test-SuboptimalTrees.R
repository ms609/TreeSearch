library("TreeTools", quietly = TRUE)

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

# --- Scores surfaced only when collapse = FALSE ---

test_that("MaximizeParsimony(collapse = TRUE) attaches no per-tree scores", {
  set.seed(3418)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                              verbosity = 0L)  # collapse = TRUE by default
  expect_null(attr(result, "scores"))
  expect_true(is.finite(attr(result, "score")))
})

test_that("MaximizeParsimony(collapse = FALSE) surfaces aligned pool scores", {
  set.seed(3418)
  result <- MaximizeParsimony(ds, maxReplicates = 3L, targetHits = 1L,
                              verbosity = 0L, collapse = FALSE,
                              control = SearchControl(poolSuboptimal = 5,
                                                      poolMaxSize = 500L))
  scores <- attr(result, "scores")
  expect_type(scores, "double")
  expect_length(scores, length(result))

  # Each tree carries a matching `score` attribute (so Suboptimality() works).
  perTree <- vapply(result, attr, double(1), "score")
  expect_equal(perTree, scores)

  # Scores match an independent re-scoring of the trees.
  expect_equal(scores, TreeLength(result, ds), tolerance = 1e-8)

  # The best of the pool equals the reported best score, and every retained
  # tree is within the requested suboptimality window.
  expect_equal(min(scores), attr(result, "score"))
  expect_true(all(scores <= attr(result, "score") + 5 + 1e-8))
})

# --- SuboptimalTrees() accessor ---

test_that("SuboptimalTrees() returns a scored, bounded pool", {
  set.seed(3418)
  trees <- SuboptimalTrees(ds, maxSuboptimal = 3, maxPool = 200L,
                           maxReplicates = 3L, targetHits = 1L,
                           verbosity = 0L)
  expect_s3_class(trees, "multiPhylo")
  scores <- attr(trees, "scores")
  expect_false(is.null(scores))
  expect_length(scores, length(trees))
  expect_true(length(trees) <= 200L)
  expect_equal(min(scores), attr(trees, "score"))
  expect_true(all(scores <= attr(trees, "score") + 3 + 1e-8))

  # Suboptimality() reads the per-tree `score` attribute we attached.
  subopt <- Suboptimality(trees)
  expect_equal(subopt, scores - min(scores))
  expect_true(all(subopt >= 0))
})

test_that("SuboptimalTrees() validates its arguments", {
  expect_error(SuboptimalTrees(ds, maxSuboptimal = -1),
               "must be a single non-negative number")
  expect_error(SuboptimalTrees(ds, maxPool = 0L),
               "must be a single positive integer")
})
