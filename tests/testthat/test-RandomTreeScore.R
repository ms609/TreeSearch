test_that("RandomMorphyTree() errors are handled", {
  expect_error(RandomMorphyTree(-1))
  expect_error(RandomMorphyTree(0))
  expect_error(RandomMorphyTree(1))
})

test_that("Two tip 'random' tree", {
  expect_equal(RandomMorphyTree(2), list(c(2, 2, 2), 0, 1))
})

test_that("RandomTreeScore() with phyDat dataset", {
  tokens <- matrix(c(
    0, "-", "-", 1, 1, 2,
    0, 1, 0, 1, 2, 2,
    0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
    dimnames = list(letters[1:3], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  
  # 1-tip dataset
  pd1 <- TreeTools::MatrixToPhyDat(tokens[1, , drop = FALSE])
  expect_equal(RandomTreeScore(pd1), 0)
  
  # 3-tip dataset: score should be a non-negative number
  set.seed(4812)
  score <- RandomTreeScore(pd)
  expect_true(is.numeric(score))
  expect_gte(score, 0)
  
  # Repeated calls give valid scores (not crashing)
  scores <- replicate(10, RandomTreeScore(pd))
  expect_true(all(is.numeric(scores)))
  expect_true(all(scores >= 0))
})

test_that("RandomTreeScore() with larger dataset", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  set.seed(7391)
  score <- RandomTreeScore(dataset)
  expect_true(is.numeric(score))
  expect_gt(score, 0)
  
  # Score should be consistent with TreeLength on the same tree
  set.seed(2048)
  tree <- TreeTools::RandomTree(dataset, root = TRUE)
  expected <- TreeLength(tree, dataset)
  expect_true(is.numeric(expected))
  expect_gt(expected, 0)
})

test_that("RandomTreeScore() backward compat with morphyObj", {
  tokens <- matrix(c(
    0, "-", "-", 1, 1, 2,
    0, "-", "-", 1, 1, 2,
    0, "-", "-", 0, 0, 0), byrow = TRUE, nrow = 3L,
    dimnames = list(letters[1:3], NULL))
  
  # One leaf
  pd <- TreeTools::MatrixToPhyDat(tokens[1, , drop = FALSE])
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(morphyObj[["nTip"]], 1L)
  expect_equal(RandomTreeScore(morphyObj), 0)
  morphyObj <- UnloadMorphy(morphyObj)
  
  # Three leaves
  pd <- TreeTools::MatrixToPhyDat(tokens)
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(RandomTreeScore(morphyObj), 3L)
  morphyObj <- UnloadMorphy(morphyObj)
})
