context("RMorphy.C[++]")

test_that("NULL pointers don't cause crash", {
  ptr <- mpl_new_Morphy()
  expect_equal(0, mpl_delete_Morphy(ptr))
  expect_true(is.na(mpl_delete_Morphy(ptr)))
})

test_that("Pointers survive garbage collection", {
  ptr <- mpl_new_Morphy()
  gc()
  expect_equal(0, mpl_delete_Morphy(ptr))
})

test_that("preorder_morphy()", {
  library('TreeTools', quietly = TRUE)
  tree <- Preorder(RootTree(BalancedTree(6), 1))
  dat <- MatrixToPhyDat(matrix(c(0, 1, 0, 1, 0, 1,
                                 0, 0, 0, 1, 1, 1), byrow = FALSE, 6,
                               dimnames = list(TipLabels(6), NULL)))
  morphyObj <- PhyDat2Morphy(dat)
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  tree$edge - 1
  expect_equal(4L, preorder_morphy(tree$edge, morphyObj))
})