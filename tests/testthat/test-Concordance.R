
dataset <- congreveLamsdellMatrices[[10]][, 1]
tree <- TreeTools::NJTree(dataset)

ConcordantInformation(tree, dataset)['noise']
TreeLength(tree, dataset, concavity = 'prof')

test_that("ConcordantInformation() works", {
  data(congreveLamsdellMatrices)
  dat <- congreveLamsdellMatrices[[10]]
  tree <- TreeTools::NJTree(dat)
  
  ci <- ConcordantInformation(tree, dat)
  expect_equal(expect_warning(Evaluate(tree, dat)), ci)
  expect_equal(TreeLength(tree, dat, concavity = 'prof'),
               unname(ci['noise']))
  expect_equal(Log2Unrooted(22), unname(ci['treeInformation']))
  expect_equal(sum(apply(PhyDatToMatrix(dat), 2, CharacterInformation)),
               unname(ci['informationContent']))
  
  dataset <- MatrixToPhyDat(cbind(setNames(c(rep(1, 11), 2:5), paste0('t', 1:15))))
  tree <- TreeTools::PectinateTree(length(dataset))
  ci <- ConcordantInformation(tree, dataset)
  expect_equal(0, unname(ci['signal']))
  expect_equal(0, unname(ci['noise']))
})
