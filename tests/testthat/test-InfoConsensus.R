# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

library("TreeTools")
library("TreeDist")

# Helpers -------------------------------------------------------------------

.CIDScorer <- TreeSearch:::.CIDScorer
.MakeCIDData <- TreeSearch:::.MakeCIDData
.CIDBootstrap <- TreeSearch:::.CIDBootstrap
.EdgeListToPhylo <- TreeSearch:::.EdgeListToPhylo
.CollapseEdge <- TreeSearch:::.CollapseEdge
.ResolveNode <- TreeSearch:::.ResolveNode
.RenumberNodes <- TreeSearch:::.RenumberNodes
.CollapseSwap <- TreeSearch:::.CollapseSwap
.ResolveSwap <- TreeSearch:::.ResolveSwap
.CollapseRefine <- TreeSearch:::.CollapseRefine
.CollapseSpecificEdge <- TreeSearch:::.CollapseSpecificEdge
.ResolveSpecificPair <- TreeSearch:::.ResolveSpecificPair
.ScoreTree <- TreeSearch:::.ScoreTree
.TopologySearch <- TreeSearch:::.TopologySearch
.RogueRefine <- TreeSearch:::.RogueRefine
.PrescreenMarginalNID <- TreeSearch:::.PrescreenMarginalNID
.PruneTrees <- TreeSearch:::.PruneTrees
.BestInsertion <- TreeSearch:::.BestInsertion
.InsertTipAtEdge <- TreeSearch:::.InsertTipAtEdge

# Small reproducible tree set
set.seed(4817)
smallTrees <- as.phylo(sample.int(100, 20), nTip = 12)


# .EdgeListToPhylo ----------------------------------------------------------

test_that(".EdgeListToPhylo returns valid phylo", {
  tr <- as.phylo(1, 10)
  edge <- tr$edge
  result <- .EdgeListToPhylo(edge[, 1], edge[, 2], tr$tip.label)
  expect_s3_class(result, "phylo")
  expect_equal(NTip(result), 10L)
  expect_equal(result$Nnode, tr$Nnode)
})


# .MakeCIDData ---------------------------------------------------------------

test_that(".MakeCIDData creates correct environment", {
  cidData <- .MakeCIDData(smallTrees, ClusteringInfoDistance,
                          smallTrees[[1]]$tip.label)
  expect_true(is.environment(cidData))
  expect_equal(length(cidData$trees), 20L)
  expect_equal(cidData$nTip, 12L)
  expect_identical(cidData$metric, ClusteringInfoDistance)
})


# .CIDScorer ----------------------------------------------------------------

test_that(".CIDScorer returns negated mean MCI", {
  cidData <- .MakeCIDData(smallTrees, ClusteringInfoDistance,
                          smallTrees[[1]]$tip.label)
  tr <- smallTrees[[1]]
  edge <- tr$edge

  score <- .CIDScorer(edge[, 1], edge[, 2], cidData)
  expected <- -mean(MutualClusteringInfo(tr, smallTrees))

  expect_equal(score, expected, tolerance = 1e-6)
})

test_that(".CIDScorer works with alternative metric", {
  cidData <- .MakeCIDData(smallTrees, MutualClusteringInfo,
                          smallTrees[[1]]$tip.label)
  tr <- smallTrees[[1]]
  edge <- tr$edge

  score <- .CIDScorer(edge[, 1], edge[, 2], cidData)
  expected <- mean(MutualClusteringInfo(tr, smallTrees))
  expect_equal(score, expected, tolerance = 1e-10)
})


# .CIDBootstrap -------------------------------------------------------------

test_that(".CIDBootstrap returns valid edgeList", {
  cidData <- .MakeCIDData(smallTrees, ClusteringInfoDistance,
                          smallTrees[[1]]$tip.label)
  tr <- smallTrees[[1]]
  edge <- tr$edge
  edgeList <- TreeTools::RenumberEdges(edge[, 1], edge[, 2])

  set.seed(7321)
  result <- .CIDBootstrap(edgeList, cidData,
                          EdgeSwapper = RootedNNISwap,
                          maxIter = 5, maxHits = 3,
                          verbosity = 0L)

  expect_length(result, 2L)
  expect_true(is.integer(result[[1]]))
  expect_true(is.integer(result[[2]]))
  expect_equal(length(result[[1]]), length(result[[2]]))
})

test_that(".CIDBootstrap restores original trees", {
  cidData <- .MakeCIDData(smallTrees, ClusteringInfoDistance,
                          smallTrees[[1]]$tip.label)
  origTrees <- cidData$trees
  tr <- smallTrees[[1]]
  edgeList <- TreeTools::RenumberEdges(tr$edge[, 1], tr$edge[, 2])

  set.seed(2199)
  .CIDBootstrap(edgeList, cidData,
                EdgeSwapper = RootedNNISwap,
                maxIter = 3, maxHits = 2,
                verbosity = 0L)

  expect_identical(cidData$trees, origTrees)
})


# InfoConsensus — core tests ------------------------------------------------

test_that("InfoConsensus rejects non-multiPhylo input", {
  expect_error(InfoConsensus(as.phylo(1, 10)),
               "multiPhylo")
})

test_that("InfoConsensus rejects single tree", {
  expect_error(InfoConsensus(c(as.phylo(1, 10))),
               "at least 2")
})

test_that("InfoConsensus runs and returns valid result", {
  set.seed(5103)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 3L, targetHits = 2L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)

  expect_s3_class(result, "phylo")
  expect_equal(NTip(result), 12L)
  resultScore <- attr(result, "score")
  expect_true(is.numeric(resultScore))
  expect_true(resultScore >= 0)
})

test_that("InfoConsensus score attribute is set", {
  set.seed(9241)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)

  expect_true(!is.null(attr(result, "score")))
  expect_true(is.numeric(attr(result, "score")))
  expect_true(attr(result, "score") >= 0)
})

test_that("CIDConsensus deprecated alias works", {
  lifecycle::expect_deprecated(
    result <- CIDConsensus(smallTrees,
                           maxReplicates = 2L, targetHits = 1L,
                           neverDrop = TRUE, collapse = FALSE,
                           verbosity = 0L)
  )
  expect_s3_class(result, "phylo")
})


# .CollapseEdge --------------------------------------------------------------

test_that(".CollapseEdge produces valid non-binary tree", {
  tr <- as.phylo(1, 10)
  edge <- tr$edge
  nTip <- NTip(tr)

  set.seed(1234)
  result <- .CollapseEdge(edge[, 1], edge[, 2], nTip)

  expect_length(result, 2L)
  newParent <- result[[1]]
  newChild <- result[[2]]

  # One fewer edge
  expect_equal(length(newParent), nrow(edge) - 1L)

  # Still a valid tree: tips 1..nTip all present as children
  expect_true(all(seq_len(nTip) %in% newChild))

  # At least one polytomy (degree > 2)
  expect_true(any(tabulate(newParent) > 2L))
})

test_that(".CollapseEdge handles star tree gracefully", {
  nTip <- 5L
  parent <- rep(nTip + 1L, nTip)
  child <- seq_len(nTip)

  result <- .CollapseEdge(parent, child, nTip)
  expect_equal(result[[1]], parent)
  expect_equal(result[[2]], child)
})


# .ResolveNode ----------------------------------------------------------------

test_that(".ResolveNode resolves a polytomy", {
  nTip <- 5L
  parent <- c(6L, 6L, 7L, 7L, 7L, 6L)
  child <-  c(1L, 7L, 2L, 3L, 4L, 5L)

  set.seed(5555)
  result <- .ResolveNode(parent, child, nTip)
  newParent <- result[[1]]
  newChild <- result[[2]]

  # One more edge than original
  expect_equal(length(newParent), length(parent) + 1L)

  # All tips still present
  expect_true(all(seq_len(nTip) %in% newChild))
})

test_that(".ResolveNode returns unchanged on binary tree", {
  tr <- as.phylo(1, 8)
  edge <- tr$edge
  nTip <- NTip(tr)

  result <- .ResolveNode(edge[, 1], edge[, 2], nTip)
  expect_equal(result[[1]], edge[, 1])
  expect_equal(result[[2]], edge[, 2])
})


# .RenumberNodes ---------------------------------------------------------------

test_that(".RenumberNodes fills gaps in node numbering", {
  nTip <- 4L
  parent <- c(5L, 5L, 7L, 7L, 5L)
  child <-  c(1L, 7L, 2L, 3L, 4L)

  result <- .RenumberNodes(parent, child, nTip)
  newParent <- result[[1]]
  newChild <- result[[2]]

  internalNew <- sort(unique(c(newParent, newChild)))
  internalNew <- internalNew[internalNew > nTip]

  # Should be contiguous from nTip+1
  expect_equal(internalNew, seq(nTip + 1L, nTip + length(internalNew)))
})


# Collapse then resolve roundtrip -------------------------------------------

test_that("Collapse then resolve preserves tip set", {
  tr <- as.phylo(1, 10)
  edge <- tr$edge
  nTip <- NTip(tr)
  origEdges <- nrow(edge)

  set.seed(7777)
  collapsed <- .CollapseEdge(edge[, 1], edge[, 2], nTip)
  expect_equal(length(collapsed[[1]]), origEdges - 1L)

  resolved <- .ResolveNode(collapsed[[1]], collapsed[[2]], nTip)
  expect_equal(length(resolved[[1]]), origEdges)

  # All tips present
  expect_true(all(seq_len(nTip) %in% resolved[[2]]))
})


# .CollapseSwap and .ResolveSwap interface -----------------------------------

test_that(".CollapseSwap follows EdgeSwapper interface", {
  tr <- as.phylo(1, 10)
  edge <- tr$edge
  set.seed(4444)
  result <- .CollapseSwap(edge[, 1], edge[, 2])
  expect_length(result, 2L)
  expect_true(is.integer(result[[1]]) || is.numeric(result[[1]]))
})

test_that(".ResolveSwap follows EdgeSwapper interface", {
  nTip <- 5L
  parent <- c(6L, 6L, 7L, 7L, 7L, 6L)
  child <-  c(1L, 7L, 2L, 3L, 4L, 5L)
  set.seed(6666)
  result <- .ResolveSwap(parent, child)
  expect_length(result, 2L)
})


# .CollapseSpecificEdge -------------------------------------------------------

test_that(".CollapseSpecificEdge collapses the targeted edge", {
  tr <- as.phylo(1, 8)
  edge <- tr$edge
  nTip <- NTip(tr)
  internalEdges <- which(edge[, 2] > nTip)

  result <- .CollapseSpecificEdge(edge[, 1], edge[, 2], internalEdges[1], nTip)
  expect_equal(length(result[[1]]), nrow(edge) - 1L)
  expect_true(all(seq_len(nTip) %in% result[[2]]))
})


# .ResolveSpecificPair -------------------------------------------------------

test_that(".ResolveSpecificPair creates a new binary split", {
  nTip <- 5L
  parent <- c(6L, 6L, 6L, 6L, 6L)
  child <-  c(1L, 2L, 3L, 4L, 5L)
  moveEdges <- c(1L, 2L)

  result <- .ResolveSpecificPair(parent, child, 6L, moveEdges, nTip)
  newParent <- result[[1]]
  newChild <- result[[2]]

  # One more edge
  expect_equal(length(newParent), length(parent) + 1L)
  # Tips 1 and 2 now share a new parent
  expect_equal(newParent[1], newParent[2])
  expect_true(newParent[1] != 6L)
})


# .CollapseRefine ---------------------------------------------------------------

test_that(".CollapseRefine can collapse edges to improve score", {
  set.seed(3210)
  goodTree <- as.phylo(1, 10)
  inputTrees <- c(rep(list(goodTree), 19),
                  list(as.phylo(99, 10)))
  class(inputTrees) <- "multiPhylo"

  cidData <- .MakeCIDData(inputTrees, ClusteringInfoDistance,
                          goodTree$tip.label)

  badTree <- as.phylo(50, 10)
  badScore <- .ScoreTree(badTree, cidData)  # internal negated MCI

  result <- .CollapseRefine(badTree, cidData, verbosity = 0L)
  resultScore <- attr(result, "score")  # internal negated MCI

  # Lower negated MCI = higher MCI = better; should be no worse
  expect_true(resultScore <= badScore + sqrt(.Machine$double.eps))
})

test_that(".CollapseRefine returns valid phylo", {
  cidData <- .MakeCIDData(smallTrees, ClusteringInfoDistance,
                          smallTrees[[1]]$tip.label)
  startTree <- smallTrees[[1]]

  result <- .CollapseRefine(startTree, cidData, verbosity = 0L)
  expect_s3_class(result, "phylo")
  expect_true(!is.null(attr(result, "score")))
  expect_equal(length(result$tip.label), 12L)
})


# InfoConsensus with collapse -------------------------------------------------

test_that("InfoConsensus collapse=TRUE produces equal-or-better score", {
  set.seed(7799)
  resultNoCollapse <- InfoConsensus(smallTrees,
                                    maxReplicates = 2L, targetHits = 1L,
                                    neverDrop = TRUE, collapse = FALSE,
                                    verbosity = 0L)
  set.seed(7799)
  resultCollapse <- InfoConsensus(smallTrees,
                                  maxReplicates = 2L, targetHits = 1L,
                                  neverDrop = TRUE, collapse = TRUE,
                                  verbosity = 0L)

  # Higher MCI = better; collapse should improve or equal
  expect_true(attr(resultCollapse, "score") >=
                attr(resultNoCollapse, "score") - sqrt(.Machine$double.eps))
})


# --- MCI scoring (.CIDScoreFast) -------------------------------------------

test_that("Internal MCI score is negative (negated mean MCI)", {
  cidData <- .MakeCIDData(smallTrees, smallTrees[[1]]$tip.label)
  tr <- smallTrees[[1]]
  edge <- tr$edge
  score <- .CIDScorer(edge[, 1], edge[, 2], cidData)
  expect_true(score <= 0)
  # Negated MCI should match -mean(MCI)
  r_mci <- mean(MutualClusteringInfo(tr, smallTrees))
  expect_equal(-score, r_mci, tolerance = 1e-6)
})


# --- .InsertTipAtEdge -------------------------------------------------------

test_that(".InsertTipAtEdge produces valid phylo", {
  tr <- as.phylo(1, 8)
  for (i in 1:nrow(tr$edge)) {
    inserted <- .InsertTipAtEdge(tr, "new_tip", i)
    expect_s3_class(inserted, "phylo")
    expect_equal(NTip(inserted), 9L)
    expect_true("new_tip" %in% inserted$tip.label)
    expect_equal(inserted$Nnode, tr$Nnode + 1L)
  }
})


# --- .PruneTrees ------------------------------------------------------------

test_that(".PruneTrees drops tips from all trees", {
  pruned <- .PruneTrees(smallTrees, "t1")
  expect_equal(length(pruned), length(smallTrees))
  expect_equal(NTip(pruned[[1]]), NTip(smallTrees[[1]]) - 1L)
  expect_false("t1" %in% pruned[[1]]$tip.label)
})

test_that(".PruneTrees returns unchanged on empty drop set", {
  result <- .PruneTrees(smallTrees, character(0))
  expect_identical(result, smallTrees)
})


# --- .ScoreTree -------------------------------------------------------------

test_that(".ScoreTree matches .CIDScorer on edge vectors", {
  cidData <- .MakeCIDData(smallTrees, smallTrees[[1]]$tip.label)
  tr <- smallTrees[[1]]
  edge <- tr$edge
  expect_equal(.ScoreTree(tr, cidData),
               .CIDScorer(edge[, 1], edge[, 2], cidData))
})


# --- Rogue dropping (toy example) ------------------------------------------

test_that("Rogue taxon correctly identified in toy example", {
  base <- ape::read.tree(text = "((a,b),(c,(d,e)));")
  set.seed(2891)
  rogue_trees <- lapply(sample.int(nrow(base$edge), 5), function(i)
    .InsertTipAtEdge(base, "r", i))
  class(rogue_trees) <- "multiPhylo"

  set.seed(2891)
  result <- InfoConsensus(rogue_trees,
                          maxReplicates = 3L, targetHits = 2L,
                          neverDrop = FALSE,
                          verbosity = 0L)
  expect_equal(attr(result, "droppedTips"), "r")
  expect_false("r" %in% result$tip.label)
})


# --- neverDrop parameter ----------------------------------------------------

test_that("neverDrop = TRUE prevents rogue dropping", {
  set.seed(4817)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = TRUE,
                          verbosity = 0L)
  expect_null(attr(result, "droppedTips"))
  expect_equal(NTip(result), 12L)
})

test_that("neverDrop protects specified tips", {
  base <- ape::read.tree(text = "((a,b),(c,(d,e)));")
  set.seed(2891)
  rogue_trees <- lapply(sample.int(nrow(base$edge), 5), function(i)
    .InsertTipAtEdge(base, "r", i))
  class(rogue_trees) <- "multiPhylo"

  set.seed(2891)
  result <- InfoConsensus(rogue_trees,
                          maxReplicates = 3L, targetHits = 2L,
                          neverDrop = c("r", "a"),
                          verbosity = 0L)
  expect_true("r" %in% result$tip.label)
  expect_true("a" %in% result$tip.label)
})


# --- maxDrop parameter -------------------------------------------------------

test_that("maxDrop limits the number of dropped tips", {
  set.seed(4817)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = FALSE,
                          maxDrop = 1L,
                          verbosity = 0L)
  nDropped <- length(attr(result, "droppedTips"))
  expect_true(nDropped <= 1L)
})



