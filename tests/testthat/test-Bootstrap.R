test_that("BootstrapTree() avoids the sample() length-1 vector trap", {
  # When only one character carries nonzero weight and that character's index
  # is not 1, deindexedChars is a length-1 vector holding that index (say 2).
  # sample(deindexedChars, ...) then samples from 1:2 rather than always
  # returning 2, so a zero-weight character could spuriously be resampled.
  captured <- new.env()
  mockSearch <- function(edgeList, dataset, ...) {
    captured$dataset <- dataset
    list(edgeList[[1]], edgeList[[2]])
  }
  testthat::local_mocked_bindings(EdgeListSearch = mockSearch,
                                  .package = "TreeSearch")

  dataset <- list(original_weight = c(0, 1))
  edgeList <- list(1, 2, 3)
  set.seed(1)
  reps <- replicate(50, {
    BootstrapTree(edgeList, dataset, maxIter = 1, maxHits = 1)
    captured$dataset[["weight"]]
  })

  expect_true(all(reps[1, ] == 0))
  expect_true(all(reps[2, ] == 1))
})
