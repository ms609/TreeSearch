test_that("ImposeConstraint() works", {
  library("TreeTools", quietly = TRUE, warn.conflict = FALSE)
  tips <- letters[1:9]
  tree <- as.phylo(1, 9, tips)
  constraint <- StringToPhyDat('0000?1111 000111111 0000??110', tips, FALSE)
  expect_equal(read.tree(text = '((a, (b, c)), (d, (e, (f, (i, (g, h))))));'),
               ImposeConstraint(tree, constraint))
})
