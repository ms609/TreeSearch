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

test_that("Unrooted, non-TreeTools start tree is accepted", {
  # ape::rtree() + ape::unroot() yields a structurally valid unrooted binary
  # tree (nrow(edge) == 2 * NTip - 3), distinct from the malformed trees
  # ape::unroot() can produce from a TreeTools `order = "preorder"` tree
  # (covered separately in test-MaximizeParsimony-features.R). Previously
  # this shape reached MakeTreeBinary() before being rooted, which misread
  # the unrooted root's legitimate degree-3 trifurcation as a polytomy,
  # corrupting the tree and surfacing as "argument is of length zero".
  set.seed(9)
  tr <- ape::unroot(ape::rtree(NTip(dataset), tip.label = names(dataset)))
  expect_false(TreeTools::TreeIsRooted(tr))
  res <- MaximizeParsimony(
    dataset, tree = tr, maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  expect_s3_class(res, "multiPhylo")
  expect_true(attr(res, "score") > 0)
})

test_that("multiPhylo input warm-starts from the whole pool", {
  set.seed(2987)
  res <- MaximizeParsimony(
    dataset, maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  # Pass a multiPhylo directly: the search resumes from a previous result.
  # One replicate cannot consume a 50-odd tree pool, and the shortfall is
  # reported rather than swallowed.
  set.seed(2987)
  expect_warning(
    warm <- MaximizeParsimony(
      dataset, tree = res,
      maxReplicates = 1L, targetHits = 1L, verbosity = 0L
    ),
    paste0("Used 1 of the ", length(res), " trees supplied")
  )
  expect_true(attr(warm, "score") <= attr(res, "score"))

  # With replicates to spare, each pool tree seeds one of them
  set.seed(2987)
  pooled <- MaximizeParsimony(
    dataset, tree = res[seq_len(min(3L, length(res)))],
    maxReplicates = 3L, targetHits = 99L, verbosity = 0L
  )
  expect_true(attr(pooled, "score") <= attr(res, "score"))
  expect_equal(attr(pooled, "replicates"), 3L)
})

test_that("Verbosity shows 'Starting tree' instead of 'Wagner'", {
  set.seed(3491)
  good <- MaximizeParsimony(
    dataset, maxReplicates = 1L, targetHits = 1L, verbosity = 0L
  )
  # MaximizeParsimony emits two streams at verbosity = 2: cli messages
  # via message() ("Strategy: ...", "Search complete: ...") and C++
  # Rprintf progress via stdout ("Starting tree score", "Converged: ...").
  # Capture both (nested, single run) so neither leaks to console.
  msg_lines <- character()
  stdout_lines <- capture.output(
    msg_lines <- capture.output(
      warm <- MaximizeParsimony(
        dataset, tree = good[[1L]],
        maxReplicates = 1L, targetHits = 1L, verbosity = 2L
      ),
      type = "message"
    )
  )
  all_out <- paste(c(msg_lines, stdout_lines), collapse = "\n")
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
