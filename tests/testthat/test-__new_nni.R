library("TreeTools")

test_that("TBR working", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  # Move single tip
  expect_equal(8, length(all_tbr(tr$edge, 12)))
  expect_equal(8, length(all_tbr(tr$edge, 11)))
  expect_equal(8, length(all_tbr(tr$edge, 10)))
  expect_equal(8, length(all_tbr(tr$edge, 7)))
  expect_equal(8, length(all_tbr(tr$edge, 6)))
  expect_equal(8, length(all_tbr(tr$edge, 3)))
  
  # Move cherry
  expect_equal(6, length(all_tbr(tr$edge, 9)))
  expect_equal(6, length(all_tbr(tr$edge, 5)))
})

