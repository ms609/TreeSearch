test_that("Groups present/contradicted produces sense", {
  library("TreeTools", quietly = TRUE) # for as.phylo
  forest <- as.phylo(1:100, 8)
  tree <- as.phylo(0, 8)
  expect_equal(PresCont(tree, forest, plot = FALSE),
               c("", "", (c(13, 8, 14, 100, 100) - c(15, 15, 15, 0, 0)) / 100))
})
