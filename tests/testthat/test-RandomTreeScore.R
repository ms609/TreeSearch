test_that("RandomMorphyTree() errors are handled", {
  expect_error(RandomMorphyTree(-1))
  expect_error(RandomMorphyTree(0))
  expect_error(RandomMorphyTree(1))
})

test_that("Two tip 'random' tree", {
  expect_equal(RandomMorphyTree(2), list(c(2, 2, 2), 0, 1))
})

test_that("RandomTreeScore() on small trees", {
  
  mo <- mpl_new_Morphy()
  expect_equal(0L, RandomTreeScore(mo))
  mpl_delete_Morphy(mo)
  
  tokens <- matrix(c(
    0, "-", "-", 1, 1, 2,
    0, "-", "-", 1, 1, 2,
    0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
    dimnames = list(letters[1:3], NULL))
  
  # One leaf
  pd <- TreeTools::MatrixToPhyDat(tokens[1, , drop = FALSE])
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(mpl_get_numtaxa(morphyObj), 1L)
  expect_equal(0, RandomTreeScore(morphyObj))
  morphyObj <- UnloadMorphy(morphyObj)
  
  # Two leaves
  pd <- TreeTools::MatrixToPhyDat(tokens[2:3, , drop = FALSE])
    morphyObj <- PhyDat2Morphy(pd)
  expect_equal(mpl_get_numtaxa(morphyObj), 2L)
  expect_equal(RandomTreeScore(morphyObj), 3L)
  morphyObj <- UnloadMorphy(morphyObj)
  
  # Three leaves
  pd <- TreeTools::MatrixToPhyDat(tokens)
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(RandomTreeScore(morphyObj), 3L)
  morphyObj <- UnloadMorphy(morphyObj)
  
})
