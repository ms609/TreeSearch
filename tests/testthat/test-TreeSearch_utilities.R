test_that("EmptyDataset()", {
  tree <- TreeTools::PectinateTree(8)
  ret <- EmptyPhyDat(tree)
  expect_equal(TipLabels(tree), names(ret))
  expect_true(inherits(ret, 'phyDat'))
})

test_that("DoNothing()", {
  expect_equal(DoNothing(c(1, 2, 3)), c(1, 2, 3))
  expect_null(DoNothing())
  
  tree <- TreeTools::PectinateTree(8)
  ret <- EmptyPhyDat(tree)
  
  expect_equal(DoNothing(tree), tree)
  expect_equal(DoNothing(ret), ret)
})
