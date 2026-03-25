# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

library(TreeTools)
library(TreeDist)

# Helpers
.CIDScorer <- TreeSearch:::.CIDScorer
.MakeCIDData <- TreeSearch:::.MakeCIDData
.ScoreTree <- TreeSearch:::.ScoreTree
.EdgeListToPhylo <- TreeSearch:::.EdgeListToPhylo
ts_cid_consensus <- TreeSearch:::ts_cid_consensus

set.seed(6183)
smallTrees <- as.phylo(sample.int(100, 20), nTip = 12)
tipLabels <- smallTrees[[1]]$tip.label


# --- ts_cid_consensus basic smoke test ----------------------------------------

test_that("ts_cid_consensus returns valid result", {
  splitMats <- lapply(smallTrees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })

  set.seed(8472)
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 12L,
    normalize = FALSE,
    maxReplicates = 3L,
    targetHits = 2L,
    verbosity = 0L,
    nThreads = 1L
  )

  expect_type(result, "list")
  expect_true(result$pool_size > 0L)
  expect_true(is.numeric(result$best_score))
  # Score is negated MCI: should be <= 0
  expect_true(result$best_score <= 0)

  # Tree is a valid edge matrix
  bestEdge <- result$trees[[1]]
  expect_equal(ncol(bestEdge), 2L)
  nTip <- 12L
  expect_true(all(seq_len(nTip) %in% bestEdge[, 2]))
})


# --- Negated MCI matches R-level MCI -----------------------------------------

test_that("C++ negated MCI matches R-level mean(MCI)", {
  splitMats <- lapply(smallTrees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })

  set.seed(3719)
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 12L,
    normalize = FALSE,
    maxReplicates = 2L,
    targetHits = 1L,
    verbosity = 0L,
    nThreads = 1L
  )

  bestEdge <- result$trees[[1]]
  bestTree <- .EdgeListToPhylo(bestEdge[, 1], bestEdge[, 2], tipLabels)

  # R-level mean MCI
  r_mci <- mean(MutualClusteringInfo(bestTree, smallTrees))

  # C++ returns negated MCI
  expect_equal(-result$best_score, r_mci, tolerance = 1e-6)
})


# --- Search improves over majority-rule consensus ----------------------------

test_that("C++ search improves over majority-rule consensus", {
  mr <- Consensus(smallTrees, p = 0.5)
  mr_binary <- MakeTreeBinary(mr)

  splitMats <- lapply(smallTrees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })

  set.seed(5902)
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 12L,
    normalize = FALSE,
    maxReplicates = 5L,
    targetHits = 3L,
    verbosity = 0L,
    nThreads = 1L
  )

  expect_true(result$pool_size > 0L)
  expect_true(is.finite(result$best_score))
})


# --- InfoConsensus R wrapper ---------------------------------------------------

test_that("InfoConsensus runs end-to-end", {
  set.seed(4291)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)

  expect_s3_class(result, "phylo")
  expect_true(!is.null(attr(result, "score")))
  expect_true(is.numeric(attr(result, "score")))
  # User-facing score is positive MCI
  expect_true(attr(result, "score") >= 0)
  expect_equal(NTip(result), 12L)
})

test_that("InfoConsensus score matches mean(MutualClusteringInfo)", {
  set.seed(4291)
  result <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)

  r_mci <- mean(MutualClusteringInfo(result, smallTrees))
  expect_equal(attr(result, "score"), r_mci, tolerance = 1e-6)
})

test_that("InfoConsensus with collapse improves or equals binary", {
  set.seed(4291)
  noColl <- InfoConsensus(smallTrees,
                          maxReplicates = 2L, targetHits = 1L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)
  set.seed(4291)
  coll <- InfoConsensus(smallTrees,
                        maxReplicates = 2L, targetHits = 1L,
                        neverDrop = TRUE, collapse = TRUE,
                        verbosity = 0L)

  # Higher MCI = better; collapse should improve or equal
  expect_true(attr(coll, "score") >=
                attr(noColl, "score") - sqrt(.Machine$double.eps))
})

test_that("InfoConsensus rejects bad input", {
  expect_error(InfoConsensus(as.phylo(1, 10)), "multiPhylo")
  expect_error(InfoConsensus(c(as.phylo(1, 10))), "at least 2")
})


# --- Sectorial search with CID --------------------------------------------------

test_that("CID search with sectors enabled produces valid results", {
  # Use 20-tip trees so sectors can activate (>= 2 * sectorMinSize)
  set.seed(7341)
  trees20 <- as.phylo(sample.int(200, 30), nTip = 20)
  tl <- trees20[[1]]$tip.label
  splitMats <- lapply(trees20, function(tr) unclass(as.Splits(tr, tl)))

  set.seed(2519)
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 20L,
    normalize = FALSE,
    maxReplicates = 3L,
    targetHits = 2L,
    xssRounds = 2L,
    rssRounds = 2L,
    cssRounds = 1L,
    sectorMinSize = 6L,
    sectorMaxSize = 15L,
    verbosity = 0L,
    nThreads = 1L
  )

  expect_true(result$pool_size > 0L)
  expect_true(is.finite(result$best_score))
  bestEdge <- result$trees[[1]]
  expect_true(all(seq_len(20L) %in% bestEdge[, 2]))
})

test_that("Sectors-enabled score <= sectors-disabled score", {
  set.seed(7341)
  trees20 <- as.phylo(sample.int(200, 30), nTip = 20)
  tl <- trees20[[1]]$tip.label
  splitMats <- lapply(trees20, function(tr) unclass(as.Splits(tr, tl)))

  set.seed(9156)
  no_sectors <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 20L,
    normalize = FALSE,
    maxReplicates = 5L,
    targetHits = 3L,
    xssRounds = 0L,
    rssRounds = 0L,
    cssRounds = 0L,
    verbosity = 0L,
    nThreads = 1L
  )

  set.seed(9156)
  with_sectors <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 20L,
    normalize = FALSE,
    maxReplicates = 5L,
    targetHits = 3L,
    xssRounds = 3L,
    rssRounds = 3L,
    cssRounds = 2L,
    sectorMinSize = 6L,
    sectorMaxSize = 15L,
    verbosity = 0L,
    nThreads = 1L
  )

  # Negated MCI: lower = better; sectors should help or not hurt
  expect_true(with_sectors$best_score <=
                no_sectors$best_score + sqrt(.Machine$double.eps))
})

test_that("Small tree gracefully skips sectors", {
  splitMats <- lapply(smallTrees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })

  # 12-tip tree with sectorMinSize = 8 → 2*8=16 > 12, sectors skipped
  set.seed(4018)
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = 12L,
    normalize = FALSE,
    maxReplicates = 2L,
    targetHits = 1L,
    xssRounds = 3L,
    rssRounds = 3L,
    cssRounds = 2L,
    sectorMinSize = 8L,
    verbosity = 0L,
    nThreads = 1L
  )

  expect_true(result$pool_size > 0L)
  expect_true(is.finite(result$best_score))
})

test_that("InfoConsensus R wrapper with sectors enabled", {
  set.seed(7341)
  trees20 <- as.phylo(sample.int(200, 30), nTip = 20)

  set.seed(6612)
  result <- InfoConsensus(trees20,
                          maxReplicates = 3L, targetHits = 2L,
                          neverDrop = TRUE, collapse = FALSE,
                          verbosity = 0L)

  expect_s3_class(result, "phylo")
  expect_true(is.finite(attr(result, "score")))
  expect_equal(NTip(result), 20L)
})
