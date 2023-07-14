test_that("TaxonInfluence() works", {
  library("TreeTools") # for phyDat manipulation
  data("congreveLamsdellMatrices", package = "TreeSearch")
  set.seed(0)
  dataset <- congreveLamsdellMatrices[[42]][1:6, ]
  expect_error(TaxonInfluence(dataset, list(list(StarTree(dataset)))), 
               " class \"phylo\"")
  
  inf <- TaxonInfluence(dataset, ratchIter = 0, startIter = 0, verb = 0)
  expect_equal(names(inf), names(dataset))
  expect_true(all(inf >= 0))
  expect_true(all(inf <= ClusteringEntropy(BalancedTree(dataset)) * 2))
  
  # Check distance can be specified
  rf <- TaxonInfluence(dataset, tree = StarTree(dataset),
                       Distance = TreeDist::RobinsonFoulds,
                       ratchIter = 0, startIter = 0, verb = 0)
  expect_true(all(rf == as.integer(rf)))
})
