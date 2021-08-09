
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
  expect_error(ConcordantInformation(tree, dataset))
  # expect_equal(0, unname(ci['signal']))
  # expect_equal(0, unname(ci['noise']))
  
  dataset <- MatrixToPhyDat(c(a = 1, b = 2, c = 1, d = 2, e = 3, f = 3))
  tree <- TreeTools::PectinateTree(dataset)
  ci <- expect_warning(ConcordantInformation(tree, dataset))
  expect_equal(c(signal = log2(3)), ci['signal'])
  expect_equal(c(noise = log2(3)), ci['noise'])
  expect_equal(c(ignored = CharacterInformation(c(0,0,1,1,2,2)) - 
                   log2(3) - log2(3)), ci['ignored'])
  
})
