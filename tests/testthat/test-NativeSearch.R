# End-to-end checks of the native custom-search scoring path, including
# `concavity` (equal weights / implied weights / profile parsimony), which is
# now baked into `dataset` by `PrepareData()` rather than threaded through
# `TreeSearch()`/`Ratchet()`/`Jackknife()`.

test_that("TreeSearch()/Ratchet() honour concavity baked into `dataset` (EW / IW / profile)", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  for (k in list(Inf, 10, "profile")) {
    prepared <- PrepareData(dataset, concavity = k)
    res <- TreeSearch(start, prepared, maxIter = 25L, maxHits = 4L, verbosity = 0)
    expect_s3_class(res, "phylo")
    # The reported score is internally consistent with re-scoring under the
    # same criterion.
    expect_equal(attr(res, "score"), TreeScore(res, prepared))
  }
})

test_that("Ratchet() honours concavity baked into `dataset`", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  prepared <- PrepareData(dataset, concavity = 10)
  res <- Ratchet(start, prepared,
                 ratchIter = 2, searchIter = 25, searchHits = 4,
                 ratchHits = 2, verbosity = 0)
  expect_s3_class(res, "phylo")
  expect_equal(attr(res, "score"), TreeScore(res, prepared))
})

test_that("The default TreeScorer requires a ParsimonyData `dataset`", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  expect_error(TreeSearch(start, dataset, maxIter = 5L, verbosity = 0),
               "ParsimonyData")
})

test_that("InitializeData/CleanUpData are deprecated but still functional", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  released <- FALSE
  expect_warning(
    res <- TreeSearch(start, dataset,
                      InitializeData = function(d) PrepareData(d, concavity = 10),
                      CleanUpData = function(d) {
                        released <<- TRUE
                        ReleaseData(d)
                      },
                      maxIter = 20L, maxHits = 3L, verbosity = 0),
    "deprecated")
  expect_true(released)
  expect_equal(attr(res, "score"),
               TreeScore(res, PrepareData(dataset, concavity = 10)))
})
