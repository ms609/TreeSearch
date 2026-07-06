library("TreeTools", quietly = TRUE)
comb11 <- PectinateTree(letters[1:11])
unrooted11 <- UnrootTree(comb11)
data11 <- cbind(upper.tri(matrix(FALSE, 11, 11))[, 3:10], 
                lower.tri(matrix(FALSE, 11, 11))[, 2:9])
rownames(data11) <- letters[1:11]
RootySwappers <- list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap)

test_that("Tree can be found", {
  skip_if_not_installed("phangorn")
  phy11 <- phangorn::phyDat(data11, type = "USER", levels = c(FALSE, TRUE))
  phy11Data <- PrepareData(phy11)
  RNGkind("Mersenne-Twister")
  set.seed(1)
  random11 <- as.phylo(17905853L, 11, letters[1:11])
  expect_error(TreeSearch(unrooted11, dataset = phy11Data))
  expect_equal_tree(comb11, TreeSearch(random11, dataset = phy11Data, maxIter = 200,
                                  EdgeSwapper = RootedTBRSwap, verbosity = 0L))
  expect_equal_tree(comb11, TreeSearch(random11, phy11Data, maxIter = 400,
                                  EdgeSwapper = RootedSPRSwap, verbosity = 0L))
  someOtherTree <- as.phylo(29235922L, 11, letters[1:11])
  expect_equal_tree(comb11, TreeSearch(someOtherTree, phy11Data, maxIter = 200,
                                  EdgeSwapper = RootedNNISwap, verbosity = 0))
  expect_equal_tree(comb11, Ratchet(random11, phy11Data, searchIter = 10, searchHits = 5,
                               swappers = RootySwappers, ratchHits = 3,
                               verbosity = 0))

  expect_false(all.equal(comb11, TreeSearch(random11, dataset = phy11Data,
                                            maxIter = 1000,
                                            stopAtPlateau = 1, verbosity = 0)))
#  TODO: Sectorial Search not working yet!
#  expect_equal(SectorialSearch(RandomTree(phy11, "a"), phy11, verbosity = -1), comb11)
})

test_that("Tree search finds shortest tree", {
  true_tree <- ape::read.tree(text = "(((((1,2),3),4),5),6);")
  malformed_tree <- ape::read.tree(text = "((((1,2),3),4),5,6);")
  dataset <- TreeTools::StringToPhyDat("110000 111000 111100", 1:6, byTaxon = FALSE)
  expect_error(TreeSearch(malformed_tree, dataset))
  start_tree <- TreeTools::RenumberTips(ape::read.tree(
    text = "(((1, 6), 3), (2, (4, 5)));"), true_tree$tip.label)
  expect_equal(TreeLength(start_tree, dataset), 6)
  preparedData <- PrepareData(dataset)
  on.exit(preparedData <- ReleaseData(preparedData))

  # NNI can reach a local optimum that SPR/TBR can escape.
  # Rooted swappers cannot move the root, so may stay at start-tree score.
  # Assert: search runs without error and score doesn't increase.
  expect_lte(attr(TreeSearch(start_tree, preparedData, EdgeSwapper = NNISwap,
                             verbosity = 0), "score"),
             TreeLength(start_tree, dataset))
  expect_equal(TreeLength(true_tree, dataset),
               attr(TreeSearch(start_tree, preparedData, EdgeSwapper = SPRSwap,
                               verbosity = -1), "score"))
  expect_equal(TreeLength(true_tree, dataset),
               attr(TreeSearch(start_tree, preparedData, EdgeSwapper = TBRSwap,
                               verbosity = -1), "score"))
  expect_lte(attr(TreeSearch(start_tree, preparedData,
                             EdgeSwapper = RootedNNISwap, verbosity = -1),
                  "score"),
             TreeLength(start_tree, dataset))
  expect_lte(attr(TreeSearch(start_tree, preparedData,
                             EdgeSwapper = RootedSPRSwap, verbosity = -1),
                  "score"),
             TreeLength(start_tree, dataset))
  expect_lte(attr(TreeSearch(start_tree, preparedData,
                             EdgeSwapper = RootedTBRSwap, verbosity = -1),
                  "score"),
             TreeLength(start_tree, dataset))
  ratchetScore <- attr(Ratchet(start_tree, preparedData,
                  swappers = list(TBRSwap, SPRSwap, NNISwap),
                  ratchIter = 3, searchHits = 5, verbosity = 0), "score")
  expect_equal(TreeLength(true_tree, dataset), ratchetScore)
})


test_that("Profile parsimony works in tree search", {
  skip_if_not_installed("phangorn")
  phy11 <- phangorn::phyDat(data11, type = "USER", levels = c(FALSE, TRUE))
  
  random11 <- as.phylo(17905853L, 11, letters[1:11]) # Rooted on "a"
  
  sillyData <- lapply(1:22, function (i) c(rep(0, i - 1), rep(1, 22 - i),
                                           rep(1, 22 - i), rep(0, i - 1)))#, sample(2, 20, replace = TRUE)-1))
  names(sillyData) <- as.character(1:22)
  dataset <- TreeTools::PhyDat(sillyData)
  readyData <- PrepareDataProfile(dataset)
  
  RNGkind("Mersenne-Twister")
  set.seed(0)
  
  rTree <- randomTree <- RandomTree(dataset, "1")
  expect_lte(TreeLength(rTree, readyData), TreeLength(rTree, dataset))

  preparedData <- PrepareData(dataset)
  quickTS <- TreeSearch(rTree, preparedData, EdgeSwapper = RootedNNISwap,
                        maxIter = 1600, maxHits = 40, verbosity = 0)
  expect_equal(42L, attr(quickTS, "score"))

  quickFitch <- Ratchet(rTree, preparedData, suboptimal = 2,
                        swappers = RootySwappers, ratchHits = 3, searchHits = 15,
                        searchIter = 100, ratchIter = 500,
                        verbosity = 0L)
  expect_equal(42, attr(quickFitch, "score"))

})

test_that("Ratchet fails gracefully", {
  expect_error(Ratchet(unrooted11, data11))
})
