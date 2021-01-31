test_that("cNNI()", {
  tr <- Preorder(root(TreeTools::BalancedTree(letters[1:7]), 'a', resolve.root = TRUE))
  expect_equal(ape::read.tree(text="(a, (b, (g, ((c, d), (e, f)))));"), # by observation, not calculation
               cNNI(tr, 1, 1))
  expect_equal(cNNI(tr, 1, 1), cNNI(tr, 1, 3))
  expect_equal(ape::read.tree(text="(a, (b, ((e, f), ((c, d), g))));"), # by observation, not calculation
               cNNI(tr, 1, 2))
  expect_equal(cNNI(tr, 1, 2), cNNI(tr, 1, 0))
  expect_equal(ape::read.tree(text="(a, (b, (d, (c, (g, (e, f))))));"), # by observation, not calculation
               cNNI(tr, 2, 1))
  expect_equal(ape::read.tree(text="(a, ((b, (c, d)), ((e, f), g)));"), # by observation, not calculation
               cNNI(tr, 3, 1))
  expect_equal(ape::read.tree(text="((a, b), ((c, d), ((e, f), g)));"), # by observation, not calculation
               cNNI(tr, 4, 1))
  suppressWarnings(RNGversion('3.5.0'))
  set.seed(0) # sample.int gives 4, 1
  expect_equal(cNNI(tr, 4, 1), cNNI(tr))
})