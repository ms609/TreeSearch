library("TreeTools", quietly = TRUE)

test_that("PaintCharacters() returns valid hex colours", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]][, 1:12]
  tree <- referenceTree

  cols <- PaintCharacters(dataset, tree)

  expect_type(cols, "character")
  expect_length(cols, 12L)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
})

test_that("PaintCharacters() returns grey when no concordant signal", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]][, 1:8]
  tree <- referenceTree

  # threshold=Inf excludes all edges → all characters grey
  cols <- PaintCharacters(dataset, tree, threshold = Inf)
  expect_true(all(cols == "#888888"))
})

test_that("PaintCharacters() accepts palette variants", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]][, 1:6]
  tree <- referenceTree

  expect_type(PaintCharacters(dataset, tree, palette = "default"),     "character")
  expect_type(PaintCharacters(dataset, tree, palette = "protanopia"),  "character")
  expect_type(PaintCharacters(dataset, tree, palette = "tritanopia"),  "character")
  grey_pal <- function(h, s) grDevices::grey(1 - s * 0.8)
  expect_type(PaintCharacters(dataset, tree, palette = grey_pal),      "character")
})
