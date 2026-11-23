test_that("TaxonInfluence() works", {
  library("TreeTools") # for phyDat manipulation
  library("TreeDist")  # for ClusteringEntropy
  data("congreveLamsdellMatrices", package = "TreeSearch")
  set.seed(0)
  dataset <- congreveLamsdellMatrices[[42]][1:6, ]
  expect_error(TaxonInfluence(dataset, list(list(StarTree(dataset)))), 
               " class \"phylo\"")
  
  inf <- TaxonInfluence(
    dataset, maxReplicates = 2L, targetHits = 1L, verbosity = 0L
  )
  expect_equal(colnames(inf), names(dataset))
  expect_true(all(inf >= 0))
  expect_true(all(inf <= ClusteringEntropy(BalancedTree(dataset)) * 2))
  expect_true(all(inf["min", ] <= inf["dwMean", ]))
  expect_true(all(inf["max", ] >= inf["dwMean", ]))
  
  # Check distance can be specified
  rf <- TaxonInfluence(dataset, tree = StarTree(dataset),
                       Distance = TreeDist::RobinsonFoulds,
                       calcWeighted = FALSE,
                       maxReplicates = 2L, targetHits = 1L,
                       verbosity = 0L)[c("min", "max"), ]
  expect_true(all(rf == as.integer(rf)))
})

test_that("TaxonInfluence() saves intermediate trees", {
  library("TreeTools") # for phyDat manipulation
  data("congreveLamsdellMatrices", package = "TreeSearch")
  set.seed(0)
  dataset <- congreveLamsdellMatrices[[42]][1:5, ]
  tree <- BalancedTree(dataset)
  
  testDir <- withr::local_tempdir()
  inf <- TaxonInfluence(
    dataset, tree, maxReplicates = 2L, targetHits = 1L, verbosity = 0L,
    savePath = paste0(testDir, "/tmp-")
  )
  expect_false(file.exists(basename(testDir)))
  expect_true(file.exists(paste0(testDir, "/tmp-5.nex")))
  expect_error(TaxonInfluence(dataset, useCache = TRUE),
               "Specify cache path using `savePath` parameter")
  expect_equal(
    expect_silent(
      TaxonInfluence(dataset, tree, savePath = paste0(testDir, "/tmp-"),
                     useCache = TRUE, verbosity = 1L)),
    inf)
})

test_that("TaxonInfluence() normalizes Distance() matrix orientation", {
  library("TreeTools", quietly = TRUE)
  tree <- as.phylo(1:2, nTip = 4)         # stands in for the reference trees
  resultTrees <- as.phylo(1:3, nTip = 4)  # stands in for a leave-one-out re-search

  # A `Distance` whose two-argument form returns dim(x) x dim(y) -- the
  # "matched labels" convention -- rather than TaxonInfluence's real
  # (mismatched-label) dim(y) x dim(x). Exercises the shape-normalization
  # rather than assuming either orientation.
  mockDistance <- function(x, y = NULL) {
    if (is.null(y)) {
      n <- length(x)
      if (n == 2) return(matrix(c(1, 2, 0, 0), 2, 2))  # rowSums = c(1, 2)
      if (n == 3) return(diag(c(1, 10, 100)))          # colSums = c(1, 10, 100)
      stop("unexpected call")
    }
    matrix(seq_len(length(x) * length(y)), nrow = length(x), ncol = length(y))
  }
  testthat::local_mocked_bindings(
    MaximizeParsimony = function(...) resultTrees, .package = "TreeSearch"
  )

  dataset <- list(a = 1, b = 2)
  inf <- TaxonInfluence(dataset, tree = tree, Distance = mockDistance,
                        calcWeighted = TRUE, verbosity = 0L)

  expect_equal(unname(inf["dwMean", "a"]), 1815 / 333)
})
