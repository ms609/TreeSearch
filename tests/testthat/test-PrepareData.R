test_that("PrepareData() validates input", {
  expect_error(PrepareData(NA))
  expect_error(PrepareData(structure(list(), class = "phyDat"), concavity = -1))
})

test_that("ReleaseData() / is.ParsimonyData()", {
  pd <- TreeTools::MatrixToPhyDat(matrix(
    c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
    dimnames = list(letters[1:4], NULL)))
  obj <- PrepareData(pd)
  expect_true(is.ParsimonyData(obj))
  expect_false(is.ParsimonyData(pd))
  expect_identical(ReleaseData(obj), obj)
})

test_that("PrepareData() builds ParsimonyData objects", {
  tokens <- matrix(c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
                   dimnames = list(letters[1:4], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  obj <- PrepareData(pd)
  expect_true(is.ParsimonyData(obj))
  expect_equal(obj[["nTip"]], 4L)
  expect_equal(0, RandomTreeScore(obj))
})

test_that("SingleCharData() builds ParsimonyData objects", {
  obj <- SingleCharData("-0-0")
  expect_true(is.ParsimonyData(obj))
  expect_equal(obj[["nTip"]], 4L)
  expect_error(SingleCharData(character(0)))
})

test_that("TreeScore() equals TreeLength() (EW, IW, profile)", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  set.seed(1)
  tree <- TreeTools::RandomTree(dataset, root = TRUE)

  # Equal weights
  expect_equal(TreeScore(tree, PrepareData(dataset)),
               TreeLength(tree, dataset))
  # Implied weights (k = 10), simple (non-extended) IW
  expect_equal(TreeScore(tree, PrepareData(dataset, concavity = 10)),
               TreeLength(tree, dataset, concavity = 10, extended_iw = FALSE))
  # Profile parsimony
  expect_message(
    expect_message(
      expect_equal(TreeScore(tree, PrepareData(dataset, concavity = "profile")),
                   TreeLength(tree, dataset, concavity = "profile")),
      "Inapplicable.*treated as ambiguous"),
    "Inapplicable.*treated as ambiguous")
})

test_that("Resampling weights change the score", {
  # Regression guard: scoring uses the *current* weights, not weights frozen at
  # construction.  Doubling every weight doubles the score; zeroing one cannot
  # increase it.
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  set.seed(1)
  tree <- TreeTools::RandomTree(dataset, root = TRUE)
  obj <- PrepareData(dataset)

  full <- TreeScore(tree, obj)
  expect_gt(full, 0)

  doubled <- obj
  doubled[["weight"]] <- obj[["weight"]] * 2L
  expect_equal(TreeScore(tree, doubled), full * 2)

  dropOne <- obj
  dropOne[["weight"]] <- c(0L, obj[["weight"]][-1])
  expect_lte(TreeScore(tree, dropOne), full)
})

test_that("Deprecated Morphy aliases still work", {
  pd <- TreeTools::MatrixToPhyDat(matrix(
    c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
    dimnames = list(letters[1:4], NULL)))
  expect_warning(obj <- PhyDat2Morphy(pd), "PrepareData")
  expect_true(is.ParsimonyData(obj))
  expect_warning(expect_true(is.morphyPtr(obj)), "is.ParsimonyData")
  expect_warning(expect_equal(UnloadMorphy(obj), 0L), "ReleaseData")
  expect_warning(SingleCharMorphy("-0-0"), "SingleCharData")
  expect_error(suppressWarnings(PhyDat2Morphy(pd, "ambiguous")))
})
