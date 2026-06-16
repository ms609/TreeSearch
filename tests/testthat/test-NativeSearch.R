# End-to-end checks of the native custom-search scoring path, including the
# `concavity` parameter (equal weights / implied weights / profile parsimony).

test_that("TreeSearch() honours concavity (EW / IW / profile)", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  for (k in list(Inf, 10, "profile")) {
    res <- TreeSearch(start, dataset, concavity = k, maxIter = 25L,
                      maxHits = 4L, verbosity = 0)
    expect_s3_class(res, "phylo")
    # The reported score is internally consistent with re-scoring under the
    # same criterion.
    expect_equal(attr(res, "score"),
                 TreeScore(res, PrepareData(dataset, concavity = k)))
  }
})

test_that("Ratchet() honours concavity", {
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  res <- Ratchet(start, dataset, concavity = 10,
                 ratchIter = 2, searchIter = 25, searchHits = 4,
                 ratchHits = 2, verbosity = 0)
  expect_s3_class(res, "phylo")
  expect_equal(attr(res, "score"),
               TreeScore(res, PrepareData(dataset, concavity = 10)))
})

test_that("Custom InitializeData ignores concavity", {
  # A user-supplied InitializeData receives the dataset unchanged; concavity is
  # not threaded into it.
  dataset <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  start <- TreeTools::NJTree(dataset)
  res <- TreeSearch(start, dataset, concavity = 10,
                    InitializeData = function(d) PrepareData(d, concavity = Inf),
                    maxIter = 20L, maxHits = 3L, verbosity = 0)
  # Scored under equal weights despite concavity = 10, because the custom
  # InitializeData built an EW dataset.
  expect_equal(attr(res, "score"),
               TreeScore(res, PrepareData(dataset)))
})
