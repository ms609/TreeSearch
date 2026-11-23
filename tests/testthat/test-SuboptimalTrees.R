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

test_that("SuboptimalTrees() keeps scores aligned under pool eviction (TA-5)", {
  # A tiny maxPool forces the diversity-eviction path; the scores attribute -- and
  # each tree's own score attribute -- must stay aligned with the evicted-down set.
  set.seed(3418)
  trees <- SuboptimalTrees(ds, maxSuboptimal = 8, maxPool = 5L,
                           maxReplicates = 4L, targetHits = 1L, verbosity = 0L)
  expect_true(length(trees) <= 5L)
  scores <- attr(trees, "scores")
  expect_length(scores, length(trees))
  expect_equal(vapply(trees, attr, double(1), "score"), scores)   # alignment held
  expect_equal(scores, TreeLength(trees, ds), tolerance = 1e-8)    # independent re-score
})

test_that("SuboptimalTrees() warns and ignores raw poolSuboptimal/poolMaxSize (B-10)", {
  # These raw SearchControl fields would silently override the values set from
  # maxSuboptimal/maxPool via MaximizeParsimony()'s dots-override-control merge.
  set.seed(3418)
  expect_warning(
    SuboptimalTrees(ds, maxSuboptimal = 2, poolSuboptimal = 9,
                    maxReplicates = 2L, targetHits = 1L, verbosity = 0L),
    "maxSuboptimal")
  expect_warning(
    SuboptimalTrees(ds, maxPool = 100L, poolMaxSize = 5L,
                    maxReplicates = 2L, targetHits = 1L, verbosity = 0L),
    "maxPool")
})
