test_that("NativeLength matches TreeLength", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  dataset <- Lobo.phy
  tree <- NJTree(dataset)
  edge <- tree[["edge"]]

  nd <- PrepareNativeData(dataset)
  expect_equal(NativeLength(edge[, 1], edge[, 2], nd),
               TreeLength(tree, dataset))

  nd_iw <- PrepareNativeData(dataset, concavity = 10)
  expect_equal(NativeLength(edge[, 1], edge[, 2], nd_iw),
               TreeLength(tree, dataset, concavity = 10, extended_iw = FALSE))
})

test_that("PrepareNativeData validates concavity", {
  data("Lobo", package = "TreeTools")
  expect_error(PrepareNativeData(Lobo.phy, concavity = 0), "must be positive")
  expect_error(PrepareNativeData(Lobo.phy, concavity = -5), "must be positive")
})

test_that("CleanNativeData is a no-op", {
  nd <- list(weight = 1:3)
  expect_identical(CleanNativeData(nd), nd)
})

test_that("NativeLength works in EdgeListSearch", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  tree <- NJTree(Lobo.phy)
  nd <- PrepareNativeData(Lobo.phy)
  tree2 <- RenumberTips(tree, names(Lobo.phy))
  edgeList <- RenumberEdges(tree2[["edge"]][, 1], tree2[["edge"]][, 2])
  startScore <- NativeLength(edgeList[[1]], edgeList[[2]], nd)
  result <- EdgeListSearch(edgeList, nd, TreeScorer = NativeLength,
                           EdgeSwapper = RootedTBRSwap,
                           maxIter = 50, maxHits = 5, verbosity = 0)
  expect_lte(result[[3]], startScore)
})

test_that("NativeBootstrap returns valid edge list", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  tree <- NJTree(Lobo.phy)
  nd <- PrepareNativeData(Lobo.phy)
  tree2 <- RenumberTips(tree, names(Lobo.phy))
  edgeList <- RenumberEdges(tree2[["edge"]][, 1], tree2[["edge"]][, 2])
  bootResult <- NativeBootstrap(edgeList, nd, EdgeSwapper = NNISwap,
                                maxIter = 20, maxHits = 5, verbosity = 0)
  expect_length(bootResult, 2)
  expect_type(bootResult[[1]], "integer")
  expect_equal(length(bootResult[[1]]), length(edgeList[[1]]))
  expect_equal(nd[["weight"]], nd[["original_weight"]])
})

test_that("TreeSearch() accepts concavity parameter", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  tree <- NJTree(Lobo.phy)
  result_ew <- TreeSearch(tree, Lobo.phy, maxIter = 5, maxHits = 2, verbosity = 0)
  expect_s3_class(result_ew, "phylo")
  result_iw <- TreeSearch(tree, Lobo.phy, concavity = 10,
                          maxIter = 5, maxHits = 2, verbosity = 0)
  expect_s3_class(result_iw, "phylo")
})

test_that("Ratchet() accepts concavity parameter", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  tree <- NJTree(Lobo.phy)
  result <- Ratchet(tree, Lobo.phy, concavity = 10,
                    ratchIter = 1, ratchHits = 1,
                    searchIter = 5, searchHits = 2, verbosity = 0)
  expect_s3_class(result, "phylo")
})
