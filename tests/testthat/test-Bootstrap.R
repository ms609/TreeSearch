test_that("BootstrapTree() resamples fractional weights without degenerating (#139)", {
  # Regression test: fractional character weights used to be truncated to
  # integer in `original_weight`, which for uniform sub-1 weights floors
  # every value to zero. tabulate(sample(integer(0), ...)) then produced an
  # all-zero resampled weight vector with no error, warning or message, so
  # the search accepted every rearrangement as tied (see #139).
  dataset <- TreeTools::StringToPhyDat(
    "1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000",
    1:7,
    byTaxon = FALSE
  )
  names(dataset) <- c(LETTERS[1:6], "out")
  attr(dataset, "weight") <- rep(0.5, length(attr(dataset, "weight")))

  preparedData <- PrepareData(dataset)
  expect_true(sum(preparedData[["original_weight"]]) > 0)

  start_tree <- ape::read.tree(text = "(((((A,D),B),E),(C,F)),out);")
  start_tree <- TreeTools::RenumberTips(start_tree, names(dataset))
  edgeList <- TreeTools::RenumberEdges(start_tree[["edge"]][, 1],
                                       start_tree[["edge"]][, 2])

  withr::local_rng_version("3.5.0")
  set.seed(0)
  res <- BootstrapTree(edgeList[1:2], preparedData,
                       maxIter = 8L, maxHits = 4L, verbosity = 0L)
  expect_type(res, "list")
  expect_length(res, 2L)
})

test_that("BootstrapTree() and JackknifeTree() refuse a zeroed original_weight", {
  dataset <- TreeTools::StringToPhyDat("1100000 1110000", 1:2, byTaxon = FALSE)
  names(dataset) <- paste0("t", seq_along(dataset))
  obj <- PrepareData(dataset)
  obj[["original_weight"]] <- integer(length(obj[["original_weight"]]))

  edgeList <- list(integer(0), integer(0))
  expect_error(BootstrapTree(edgeList, obj), "sums to zero")
  expect_error(JackknifeTree(edgeList, obj), "sums to zero")
})
