context("tree_properties.R")

test_that('Rooting and partition counting',{
  set.seed(0)
  
  expect_false(TreeIsRooted(ape::read.tree(text='(a, b, c);')))
  expect_true(TreeIsRooted(ape::read.tree(text='(a, (b, c));')))
  
  tree8 <- ape::rtree(8L)
  expect_true(TreeIsRooted(tree8))
  expect_false(TreeIsRooted(ape::unroot(tree8)))
  
  expect_equal(5L, NPartitions(tree8))
  expect_equal(5L, NPartitions(ape::unroot(tree8)))
  
  expect_equal(c(5L, 5L), NPartitions(list(tree8, ape::unroot(tree8))))
  expect_error(NPartitions(5L))
})