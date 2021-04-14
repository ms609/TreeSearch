test_that("ConcordantInformation() works", {
  data(congreveLamsdellMatrices)
  dat <- congreveLamsdellMatrices[[10]]
  tree <- TreeTools::NJTree(dat)
  
  ci <- ConcordantInformation(tree, dat)
  expect_equal(expect_warning(Evaluate(tree, dat)),
               ci)
  expect_equal(TreeLength(tree, dat, concavity = 'prof'),
               unname(ci['noise']))
  expect_equal(Log2Unrooted(22), unname(ci['treeInformation']))
  expect_equal(sum(apply(PhyDatToMatrix(dat), 2, CharacterInformation)),
               unname(ci['informationContent']))
})
