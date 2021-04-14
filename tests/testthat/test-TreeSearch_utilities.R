test_that("EmptyDataset()", {
  tree <- TreeTools::PectinateTree(8)
  ret <- EmptyPhyDat(tree)
  expect_equal(TipLabels(tree), names(ret))
  expect_true(inherits(ret, 'phyDat'))
})
