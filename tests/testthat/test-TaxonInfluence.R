test_that("TaxonInfluence() works", {
  library("TreeTools") # for phyDat manipulation
  data("congreveLamsdellMatrices", package = "TreeSearch")
  set.seed(0)
  dataset <- congreveLamsdellMatrices[[42]][1:6, ]
  expect_error(TaxonInfluence(dataset, list(list(StarTree(dataset)))), 
               " class \"phylo\"")
  
  inf <- TaxonInfluence(dataset, ratchIter = 0, startIter = 0, verb = 0)
  expect_equal(colnames(inf), names(dataset))
  expect_true(all(inf >= 0))
  expect_true(all(inf <= ClusteringEntropy(BalancedTree(dataset)) * 2))
  expect_true(all(inf["min", ] <= inf["dwMean", ]))
  expect_true(all(inf["max", ] >= inf["dwMean", ]))
  
  # Check distance can be specified
  rf <- TaxonInfluence(dataset, tree = StarTree(dataset),
                       Distance = TreeDist::RobinsonFoulds,
                       calcWeighted = FALSE,
                       ratchIter = 0, startIter = 0, verb = 0)[
                         c("min", "max"), ]
  expect_true(all(rf == as.integer(rf)))
})

test_that("TaxonInfluence() saves intermediate trees", {
  library("TreeTools") # for phyDat manipulation
  data("congreveLamsdellMatrices", package = "TreeSearch")
  set.seed(0)
  dataset <- congreveLamsdellMatrices[[42]][1:5, ]
  
  testDir <- tempdir()
  on.exit(unlink(testDir))
  inf <- TaxonInfluence(dataset, ratchIter = 0, startIter = 0, verb = 0,
                        saveTo = paste0(testDir, "/tmp-"))
  expect_true(file.exists(paste0(testDir, "/tmp-5.nex")))
})
