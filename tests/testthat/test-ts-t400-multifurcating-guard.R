# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# T-400: `TreeState::init_from_edge()` derives every node count from the edge
# count alone, which holds only for a binary tree.  A multifurcating tree wrote
# past the end of the topology arrays when the edge count was odd, and rooted
# the tree on a real tip when it was even -- leaving a one-element postorder
# whose downpass read the words immediately *before* the state buffer.  Either
# way the caller got a plausible number rather than an error, so both parities
# are exercised below.

# `(a,(b,((e,f),(g,h),(c,d))));`  8 tips, 6 internal nodes: 13 edges (odd)
.OddPolytomy <- function() {
  ape::read.tree(text = "(a,(b,((e,f),(g,h),(c,d))));")
}

# `(a,(b,(c,(d,(e,f,g,h)))));`  8 tips, 5 internal nodes: 12 edges (even)
.EvenPolytomy <- function() {
  ape::read.tree(text = "(a,(b,(c,(d,(e,f,g,h)))));")
}

.EightTaxonData <- function() {
  MatrixToPhyDat(matrix(
    c("0", "0", "1", "1", "0", "1", "0", "1",
      "0", "1", "0", "1", "1", "1", "0", "0",
      "1", "1", "1", "0", "0", "0", "1", "0",
      "0", "0", "0", "1", "1", "0", "1", "1"),
    nrow = 8, dimnames = list(letters[1:8], NULL)))
}

test_that("The test polytomies have the edge counts the guard must handle", {
  # Both parities must be covered: they failed by different mechanisms.
  expect_equal(dim(.OddPolytomy()[["edge"]])[1], 13L)
  expect_equal(dim(.EvenPolytomy()[["edge"]])[1], 12L)
})

test_that("TreeLength() rejects a multifurcating tree in a list", {
  dat <- .EightTaxonData()
  for (tr in list(.OddPolytomy(), .EvenPolytomy())) {
    # Length-1 sets are the dangerous case: a heterogeneous set was already
    # caught, accidentally, by the differing-edge-count check.
    expect_error(TreeLength(structure(list(tr), class = "multiPhylo"), dat),
                 "must be binary")
    expect_error(TreeLength(list(tr, tr), dat), "must be binary")
  }
})

test_that("TreeLength() rejects a multifurcating tree through `[` and `[[`", {
  dat <- .EightTaxonData()
  tr <- .OddPolytomy()
  trees <- structure(list(tr, tr), class = "multiPhylo")
  # `[[` dispatches to the phylo method, which was already guarded; `[` keeps
  # the multiPhylo class and reached the kernel.
  expect_error(TreeLength(trees[[1]], dat), "must be binary")
  expect_error(TreeLength(trees[1], dat), "must be binary")
})

test_that("Repeated scoring of one tree gives one answer", {
  # The score was read from memory before the state buffer, so identical calls
  # could disagree.  Whatever the answer is, it must not vary between calls.
  dat <- .EightTaxonData()
  trees <- structure(list(.OddPolytomy()), class = "multiPhylo")
  outcomes <- vapply(seq_len(5), function(i) {
    tryCatch(paste(TreeLength(trees, dat), collapse = ","),
             error = function(e) conditionMessage(e))
  }, character(1))
  expect_length(unique(outcomes), 1L)
})

test_that("CharacterLength() rejects a multifurcating tree", {
  dat <- .EightTaxonData()
  expect_error(CharacterLength(.OddPolytomy(), dat), "must be binary")
  expect_error(CharacterLength(.EvenPolytomy(), dat), "must be binary")
})

test_that("TreeScore() and EdgeListScore() reject a multifurcating tree", {
  dat <- PrepareData(.EightTaxonData())
  tr <- .OddPolytomy()
  expect_error(TreeScore(tr, dat), "must be binary")
  expect_error(EdgeListScore(tr[["edge"]][, 1], tr[["edge"]][, 2], dat),
               "must be binary")
  ev <- .EvenPolytomy()
  expect_error(TreeScore(ev, dat), "must be binary")
  expect_error(EdgeListScore(ev[["edge"]][, 1], ev[["edge"]][, 2], dat),
               "must be binary")
})

test_that("The scoring kernel itself refuses a multifurcating edge matrix", {
  # The R-level guards above are convenience; this is the boundary that every
  # other kernel entry point sits behind.
  dat <- .EightTaxonData()
  tr <- RenumberTips(Renumber(.OddPolytomy()), names(dat))
  at <- attributes(dat)
  tipData <- matrix(unlist(dat, use.names = FALSE), nrow = length(dat),
                    byrow = TRUE)
  expect_error(
    TreeSearch:::ts_fitch_score(tr[["edge"]], at[["contrast"]], tipData,
                                TreeSearch:::.ScaleWeight(at[["weight"]]),
                                at[["levels"]]),
    "must be binary")
})

test_that("The kernel refuses malformed edge lists of binary length", {
  # A polytomy is not the only edge list that would send the kernel out of
  # bounds; these have the edge count of a four-tip binary tree (root = 5) but
  # a shape it cannot index.
  dat <- MatrixToPhyDat(matrix(c("0", "0", "1", "1",
                                 "0", "1", "0", "1",
                                 "1", "1", "0", "0"),
                               nrow = 4, dimnames = list(letters[1:4], NULL)))
  at <- attributes(dat)
  tipData <- matrix(unlist(dat, use.names = FALSE), nrow = 4, byrow = TRUE)
  Score <- function(edge) {
    TreeSearch:::ts_fitch_score(edge, at[["contrast"]], tipData,
                                TreeSearch:::.ScaleWeight(at[["weight"]]),
                                at[["levels"]])
  }
  Edge <- function(...) matrix(c(...), ncol = 2, byrow = TRUE)

  # The same topology, accepted and scored as `(a,(b,(c,d)));`
  expect_equal(Score(Edge(5, 1, 5, 6, 6, 2, 6, 7, 7, 3, 7, 4)),
               TreeLength(ape::read.tree(text = "(a,(b,(c,d)));"), dat))
  # One node claimed as a child twice, leaving another with no parent
  expect_error(Score(Edge(5, 1, 5, 6, 6, 2, 6, 7, 7, 3, 7, 3)), "must be binary")
  # The root claimed as a child
  expect_error(Score(Edge(5, 1, 5, 6, 6, 2, 6, 5, 7, 3, 7, 4)), "must be binary")
  # A tip used as a parent
  expect_error(Score(Edge(5, 1, 5, 6, 6, 2, 6, 7, 1, 3, 1, 4)), "must be binary")
})

test_that("A non-binary `startEdge` is refused on the main thread", {
  # init_from_edge also runs on a search worker, where a throw would terminate
  # the session rather than raise an R error, so the driven search screens
  # start trees before dispatching.  nThreads = 2 covers the threaded path.
  dat <- MatrixToPhyDat(matrix(
    c("0", "0", "0", "0", "0", "1", "1", "1",
      "0", "0", "1", "1", "1", "0", "0", "1",
      "0", "1", "0", "1", "1", "0", "1", "0"),
    nrow = 8, dimnames = list(letters[1:8], NULL)))
  ds <- make_ts_data(dat)
  Driven <- function(edge, nThreads = 1L) {
    TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = 2L, ratchetCycles = 1L, verbosity = 0L,
      nThreads = nThreads, startEdge = edge)
  }
  binary <- RenumberTips(Preorder(BalancedTree(letters[1:8])),
                         names(dat))[["edge"]]
  expect_error(Driven(.OddPolytomy()[["edge"]]), "binary")
  expect_error(Driven(.EvenPolytomy()[["edge"]]), "binary")
  expect_error(Driven(.EvenPolytomy()[["edge"]], nThreads = 2L), "binary")
  # A binary start must still be accepted, on both paths: a search has to have
  # run and scored something, which an empty `scores` would not show.
  for (nThreads in c(1L, 2L)) {
    scores <- Driven(binary, nThreads = nThreads)[["scores"]]
    expect_gt(length(scores), 0L)
    expect_true(all(is.finite(scores)))
  }
})

test_that("Binary trees are unaffected by the guard", {
  dat <- .EightTaxonData()
  pd <- PrepareData(dat)
  for (tr in lapply(list(MakeTreeBinary(.OddPolytomy()),
                         MakeTreeBinary(.EvenPolytomy()),
                         BalancedTree(letters[1:8]),
                         PectinateTree(letters[1:8])),
                    Preorder)) {
    score <- TreeLength(tr, dat)
    expect_true(is.finite(score))
    expect_equal(unname(TreeLength(structure(list(tr), class = "multiPhylo"),
                                   dat)),
                 score)
    expect_equal(sum(CharacterLength(tr, dat, compress = TRUE) *
                       attr(dat, "weight")),
                 score)
    expect_equal(TreeScore(RenumberTips(tr, names(dat)), pd), score)
  }
})
