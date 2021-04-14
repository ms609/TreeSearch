test_that("Errors are handled", {
  tokens <- matrix(c(
   0, '-', '-', 1, 1, 2,
   0, '-', '-', 0, 0, 0), byrow = TRUE, nrow = 2L,
   dimnames = list(letters[1:2], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(3, RandomTreeScore(morphyObj))
  morphyObj <- UnloadMorphy(morphyObj)
  
  expect_error(RandomMorphyTree(-1))
  expect_error(RandomMorphyTree(0))
  expect_error(RandomMorphyTree(1))
  
})

test_that("Random tree score on small trees", {
  
  mo <- mpl_new_Morphy()
  expect_equal(0L, RandomTreeScore(mo))
  mpl_delete_Morphy(mo)
  
  
  tokens <- matrix(c(
    0, '-', '-', 1, 1, 2,
    0, '-', '-', 1, 1, 2,
    0, '-', '-', 0, 0, 0), byrow = TRUE, nrow = 3L,
    dimnames = list(letters[1:3], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(3, RandomTreeScore(morphyObj))
  morphyObj <- UnloadMorphy(morphyObj)
  
})

test_that("Two tip 'random' tree", {
  expect_equal(list(c(2, 2, 2), 0, 1), RandomMorphyTree(2))
})
