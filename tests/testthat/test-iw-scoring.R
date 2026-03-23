test_that("IW Scoring", {
  library("TreeTools", quietly = TRUE)
  data("Lobo", package = "TreeTools")
  dataset <- Lobo.phy
  tree <- NJTree(dataset)

  concavity <- 4.5

  tree <- Preorder(RenumberTips(tree, names(dataset)))
  weight <- attr(dataset, "weight")
  minLength <- MinimumLength(dataset, compress = TRUE)

  # Verify IW score matches manually computed value.
  # Reference: per-character step counts from C++ engine, IW formula applied in R.
  charSteps <- CharacterLength(tree, dataset, compress = TRUE)
  homoplasies <- charSteps - minLength
  fit <- homoplasies / (homoplasies + concavity)
  manualIW <- sum(fit * weight)

  expect_equal(TreeLength(tree, dataset, concavity = concavity,
                         extended_iw = FALSE), manualIW)

  # Sanity check: IW score should be positive for a non-trivial dataset

  expect_gt(manualIW, 0)
})
