library("TreeTools", quietly = TRUE)

# Issue #136 (A15-04): Ratchet()'s early-exit paths (stopAtScore= met either
# before or during search) skipped the bookkeeping that its return value
# depends on, so the returned tree and its "score" attribute could disagree.
# The property that must hold on every exit path: a tree's own recomputed
# TreeLength() must equal its "score" attribute.

trueTree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
startTree <- TreeTools::RenumberTips(ape::read.tree(
  text = "(((1, 6), 3), (2, (4, 5)));"), trueTree$tip.label)
startScore <- TreeLength(startTree, dataset)
trueScore <- TreeLength(trueTree, dataset)
preparedData <- PrepareData(dataset)

test_that("Ratchet(stopAtScore=) already met on entry returns a consistent tree", {
  result <- Ratchet(startTree, preparedData, stopAtScore = startScore,
                    verbosity = 0)
  expect_false(is.null(attr(result, "score")))
  expect_equal(TreeLength(result, dataset), attr(result, "score"))
  expect_equal(attr(result, "score"), startScore)

  resultAll <- Ratchet(startTree, preparedData, stopAtScore = startScore,
                       returnAll = TRUE, verbosity = 0)
  expect_s3_class(resultAll, "multiPhylo")
  expect_length(resultAll, 1)
  expect_equal(TreeLength(resultAll[[1]], dataset), attr(resultAll[[1]], "score"))
})

test_that("Ratchet(stopAtScore=) met mid-search returns a consistent tree", {
  oldSeed <- if (exists(".Random.seed", .GlobalEnv)) .GlobalEnv[[".Random.seed"]] else NULL
  on.exit(if (is.null(oldSeed)) rm(".Random.seed", envir = .GlobalEnv) else
    assign(".Random.seed", oldSeed, envir = .GlobalEnv))
  set.seed(1)

  result <- Ratchet(startTree, preparedData, stopAtScore = trueScore,
                    swappers = list(TBRSwap, SPRSwap, NNISwap),
                    ratchIter = 3, searchHits = 5, verbosity = 0)
  expect_equal(attr(result, "score"), trueScore)
  # This is the assertion that failed pre-fix: the tree returned carried the
  # improved score but was, independently, the untouched input tree.
  expect_equal(TreeLength(result, dataset), attr(result, "score"))

  resultAll <- Ratchet(startTree, preparedData, stopAtScore = trueScore,
                       swappers = list(TBRSwap, SPRSwap, NNISwap),
                       ratchIter = 3, searchHits = 5, returnAll = TRUE,
                       verbosity = 0)
  expect_s3_class(resultAll, "multiPhylo")
  # A stopAtScore hit mid-search can only ever bank the single hitting
  # candidate: every earlier iteration scored above stopAtScore + suboptimal,
  # so only the BREAK-path forest slot survives the keepers filter.
  expect_length(resultAll, 1)
  expect_equal(TreeLength(resultAll[[1]], dataset), attr(resultAll[[1]], "score"))
})

test_that("MultiRatchet() survives an already-met stopAtScore", {
  expect_silent(result <- MultiRatchet(startTree, preparedData,
                                       stopAtScore = startScore,
                                       nSearch = 2, verbosity = 0)
  )
  expect_s3_class(result, "multiPhylo")
  for (phy in result) {
    expect_equal(TreeLength(phy, dataset), attr(phy, "score"))
  }
})
