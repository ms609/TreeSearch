context("pp_exact")

# TODO this test was recovered from a stash and requires updating -- 
# or may be obselete.
test_that("Profile score correct for small trees", {
  library("TreeTools", quietly = TRUE)
  tree <- as.phylo(200, 9)
  
  mataset <- matrix(c(
    1, 1, 1, 1, 0, 0, 0, 0, 0, # 3 steps
    1, 0, 0, 1, 0, 0, 1, 0, 0, # 2 steps
    1, 0, 0, 1, 0, 0, 1, 0, 0, # 2 steps again [duplicated]
    0, 1, 0, 0, 0, 0, 0, 1, 1, # 1 step
    2, 1, 1, 1, 1, 1, 1, 1, 1),# 1 step; non-informative
    nrow = 9, dimnames = list(paste0("t", 1:9), NULL))
    

  dataset <- MatrixToPhyDat(mataset)

  # EW score = 3 + 2 + 2 + 1 + 1 = 9
  expect_equal(9, TreeLength(tree, dataset))
  
  # With integer-step profile tables, profile scoring should equal EW scoring
  expect_equal(sum(CharacterLength(tree, dataset, compress = TRUE) *
                     attr(dataset, "weight")),
               TreeLength(tree, dataset))
})


test_that("Profile score can be calculated from real data", {
  data(referenceTree)
  data(congreveLamsdellMatrices)
  tree <- referenceTree
  dataset <- PrepareDataProfile(congreveLamsdellMatrices[[1]])
  expect_equal(TreeLength(tree, dataset), 
               sum(CharacterLength(tree, dataset, compress = TRUE) *
                     attr(dataset, "weight")))
  score <- TreeLength(tree, dataset, "profile")

  # Check score hasn't materially changed:
  # 511.732 is "previous value"; not manually checked.
  expect_equal(511.732, score, tolerance = 0.01)
})
