context("Starting tree (warm-start)")

ts_score <- function(tree, ds) {
  TreeSearch:::ts_fitch_score(
    tree$edge,
    attr(ds, "contrast"),
    matrix(unlist(ds, use.names = FALSE), nrow = length(ds), byrow = TRUE),
    attr(ds, "weight"),
    attr(ds, "levels")
  )
}

data("inapplicable.phyData", package = "TreeSearch")
dataset <- inapplicable.phyData[["Vinther2008"]]

test_that("User-supplied tree is used as starting topology", {
  # Get a good starting tree
  set.seed(6714)
  baseline <- MaximizeParsimony(
    dataset, maxReplicates = 3L, targetHits = 1L, verbosity = 0L
  )
  best_score <- attr(baseline, "score")
  good_tree <- baseline[[1L]]

  # Warm-start from the good tree
  set.seed(6714)
  warm <- MaximizeParsimony(
    dataset, tree = good_tree,
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  warm_score <- attr(warm, "score")

  # Warm-start should find score at least as good

  expect_true(warm_score <= best_score)
})

test_that("multiPhylo input extracts first tree", {
  set.seed(2987)
  res <- MaximizeParsimony(
    dataset, maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  # Pass multiPhylo directly — should extract [[1]]
  set.seed(2987)
  warm <- MaximizeParsimony(
    dataset, tree = res,
    maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  expect_true(attr(warm, "score") <= attr(res, "score"))
})

test_that("Verbosity shows 'Starting tree' instead of 'Wagner'", {
  set.seed(3491)
  good <- MaximizeParsimony(
    dataset, maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  # Capture verbose output
  out <- capture.output({
    warm <- MaximizeParsimony(
      dataset, tree = good[[1L]],
      maxReplicates = 1L, targetHits = 1L, verbosity = 2L
    )
  }, type = "message")
  # Also capture stdout (Rprintf goes to stdout)
  out2 <- capture.output({
    warm2 <- MaximizeParsimony(
      dataset, tree = good[[1L]],
      maxReplicates = 1L, targetHits = 1L, verbosity = 2L
    )
  })
  all_out <- paste(c(out, out2), collapse = "\n")
  expect_true(grepl("Starting tree score", all_out))
})

test_that("Without starting tree, default Wagner path works", {
  set.seed(5172)
  res <- MaximizeParsimony(
    dataset, maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  expect_true(attr(res, "score") > 0)
  expect_s3_class(res, "multiPhylo")
})

test_that("Starting tree with IW mode works", {
  set.seed(8456)
  good <- MaximizeParsimony(
    dataset, concavity = 10,
    maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  set.seed(8456)
  warm <- MaximizeParsimony(
    dataset, tree = good[[1L]], concavity = 10,
    maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  expect_true(attr(warm, "score") <= attr(good, "score"))
})
