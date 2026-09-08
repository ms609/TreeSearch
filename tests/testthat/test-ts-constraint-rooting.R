# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

## T-384: constraint mapping must not depend on the tree's rooting.
##
## A constraint split is an *unrooted* bipartition, so a tree displays it
## whenever EITHER side is a rooted clade.  build_constraint() canonicalizes
## every split mask to exclude tip 0, and map_constraint_nodes() used to accept
## only that side as a clade — true of a tip-0-rooted tree and of no other.  In
## any other rooting the mapping returned -1, which regraft_violates_constraint()
## reads as "the tree already violates", answering by rejecting EVERY regraft:
## a perfectly valid tree the search cannot move away from.  No wrong answer is
## returned — the constraint stays honoured either way — so what this file
## asserts is that the search MOVES, not that its output is legal.
##
## MaximizeParsimony() is shielded by an accident of its input normalisation:
## `R/MaximizeParsimony.R` re-roots any start tree whose root's first child is
## not a tip, and Preorder() orders children by smallest descendant tip, so a
## tip child in that slot can only be tip 1 — a user tree therefore always
## arrives rooted on tip 0.  The engine's own producers (fuse, sector, drift,
## nni-perturb, the parallel driver) map constraints on trees they did not root,
## which is why the contract is exercised here through ts_driven_search(), the
## entry point that takes the rooting as given.

library("TreeTools")

# tbrOnlyRun() lives in helper-ts.R: everything that could rescue a
# rooting-dependent TBR is switched off there, so the scores below are TBR's.

# Phases whose timing must be zero, so a passing test cannot be one that quietly
# searched its way around the mapping.  Each is guarded by an explicit `> 0`
# test in driven_search(), so zero ms means it did not run.  The sectorial
# counters are deliberately absent: xss_search() is *called* unconditionally and
# its round loop is empty at `xssRounds = 0`, so `xss_ms` records phase-lap
# overhead (~0.005 ms) rather than any searching.
mutedPhases <- c("ratchet_ms", "nni_perturb_ms", "drift_ms", "anneal_ms",
                 "prune_reinsert_ms", "fuse_ms")

# Does `tree` display the bipartition `inGroup` | rest?
DisplaysGroup <- function(edge, labels, inGroup) {
  tree <- structure(
    list(edge = edge, Nnode = length(labels) - 1L, tip.label = labels),
    class = "phylo")
  splits <- as.logical(as.Splits(tree, tipLabels = labels))
  if (!is.matrix(splits)) splits <- matrix(splits, nrow = 1)
  target <- labels %in% inGroup
  any(apply(splits, 1, function(row) {
    all(row == target) || all(row == !target)
  }))
}

test_that("T-384: constrained TBR searches from every rooting of one tree", {
  dataset <- congreveLamsdellMatrices[[1]]
  labels <- names(dataset)
  nTip <- length(labels)
  ds <- make_ts_data(dataset)

  # One constraint split; which side is stored is immaterial, as
  # build_constraint() canonicalizes it to exclude tip 0.
  inGroup <- labels[1:6]
  splitMatrix <- matrix(as.integer(labels %in% inGroup), nrow = 1)

  # A constraint-compliant but deliberately poor start tree: an arbitrary
  # resolution of the constraint itself.  MakeTreeBinary() resolves at random,
  # so the seed is set here, not just before the searches.
  constraint <- ape::read.tree(text = paste0(
    "((", paste(inGroup, collapse = ","), "),(",
    paste(labels[7:nTip], collapse = ","), "));"))
  set.seed(384)
  base <- MakeTreeBinary(constraint)
  startScore <- TreeLength(base, dataset)
  expect_true(DisplaysGroup(Preorder(RenumberTips(base, labels))[["edge"]],
                            labels, inGroup))

  scores <- numeric(nTip)
  compliant <- logical(nTip)
  leaked <- character(0)
  for (k in seq_len(nTip)) {
    rooted <- Preorder(RenumberTips(RootTree(base, labels[k]), labels))
    # Guard the premise: every rooting is the same tree unrooted, so all nTip
    # searches start from the same score and the same displayed splits.
    expect_equal(TreeLength(rooted, dataset), startScore)

    set.seed(384)
    result <- tbrOnlyRun(ds, rooted[["edge"]], splitMatrix)
    scores[[k]] <- result$best_score
    timings <- unlist(result$timings)
    leaked <- c(leaked, mutedPhases[timings[mutedPhases] > 0])
    compliant[[k]] <- all(vapply(result$trees, DisplaysGroup, logical(1),
                                 labels = labels, inGroup = inGroup))
  }

  # TBR alone is responsible for the scores below.
  expect_equal(unique(leaked), character(0))

  # The constraint is honoured in every rooting — before the fix as well as
  # after.  The defect cost search progress, never legality.
  expect_true(all(compliant))

  # The defect: from a rooting whose canonical (tip-0-excluded) side is not a
  # clade, every regraft was rejected and the start tree came back unimproved.
  # Pre-fix that held for the 16 of 22 rootings that place the root outside the
  # constrained group; post-fix, 0 of 22.
  expect_equal(sum(scores >= startScore), 0L,
               info = paste("start", startScore, "finals",
                            paste(scores, collapse = " ")))
})
