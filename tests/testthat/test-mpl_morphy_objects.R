test_that("PhyDat2Morphy() errors", {
  expect_error(PhyDat2Morphy(NA))
})

test_that("UnloadMorphy() errors and no-ops", {
  expect_error(UnloadMorphy(NA))
  pd <- TreeTools::MatrixToPhyDat(matrix(
    c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
    dimnames = list(letters[1:4], NULL)))
  morphyObj <- PhyDat2Morphy(pd)
  expect_true(is.morphyPtr(morphyObj))
  expect_equal(UnloadMorphy(morphyObj), 0L)
})

test_that("PhyDat2Morphy() builds native objects", {
  tokens <- matrix(c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
                   dimnames = list(letters[1:4], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  morphyObj <- PhyDat2Morphy(pd)
  expect_true(is.morphyPtr(morphyObj))
  expect_equal(morphyObj[["nTip"]], 4L)
  expect_equal(0, RandomTreeScore(morphyObj))
})

test_that("Non-inapplicable gap treatments are unsupported", {
  tokens <- matrix(c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
                   dimnames = list(letters[1:4], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  expect_error(PhyDat2Morphy(pd, "ambig"))
  expect_error(PhyDat2Morphy(pd, "eXt"))
  expect_error(PhyDat2Morphy(pd, "ERROR")) # unrecognized treatment
})

test_that("SingleCharMorphy() builds native objects", {
  morphyObj <- SingleCharMorphy("-0-0")
  expect_true(is.morphyPtr(morphyObj))
  expect_equal(morphyObj[["nTip"]], 4L)
  expect_error(SingleCharMorphy("-0-0", "ERROR"))
  expect_error(SingleCharMorphy("-0-0", "eXt"))
})

test_that("Resampling weights change the score", {
  # Regression guard (T-2.0 native resampling): the custom-search scoring path
  # must use the *current* character weights, not weights frozen at object
  # construction.  Doubling every weight must double the score; if weights were
  # frozen the two scores would be identical.
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  set.seed(1)
  tree <- TreeTools::RandomTree(dataset, root = TRUE)
  morphyObj <- PhyDat2Morphy(dataset)

  full <- MorphyTreeLength(tree, morphyObj)
  expect_gt(full, 0)

  w <- TreeSearch:::.MorphyWeight(morphyObj)
  doubled <- TreeSearch:::.SetMorphyWeight(morphyObj, w * 2)
  expect_equal(MorphyTreeLength(tree, doubled), full * 2)

  # Zeroing a weight cannot increase the score
  dropOne <- TreeSearch:::.SetMorphyWeight(morphyObj,
                                           c(0, w[-1]))
  expect_lte(MorphyTreeLength(tree, dropOne), full)
})
