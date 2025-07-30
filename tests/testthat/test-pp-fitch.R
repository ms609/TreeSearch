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
  
  at <- attributes(dataset)
  characters <- PhyToString(dataset, ps = "", useIndex = FALSE,
                            byTaxon = FALSE, concatenate = FALSE)
  weight <- at$weight
  morphyObjects <- lapply(characters, SingleCharMorphy)
  on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)))
  
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- at$contrast
  simpleCont <- ifelse(rowSums(cont) == 1,
                       apply(cont != 0, 1, function (x) at$levels[x][1]),
                       "?")
  inappLevel <- at$levels == "-"
  
  unlisted <- unlist(dataset, use.names = FALSE)
  charSeq <- seq_len(nChar) - 1L
  
  tokenMatrix <- matrix(simpleCont[unlisted], nChar, 9, byrow = FALSE)
  profileTables <- apply(tokenMatrix, 1, table)
  if (inherits(profileTables, "matrix")) {
    profileTables <- lapply(seq_len(ncol(profileTables)), function (i) profileTables[, i])
  }
  data("profiles", package = "TreeSearch")
  profileCost <- lapply(profileTables, function (x) {
    x <- sort(x[x > 1])
    n <- length(x)
    prof <- switch(n,
                   0,
                   profiles[[sum(x)]][[n]][[x[1] - 1L]]
    )
  })
  profileExtra <- lapply(profileCost, function (x)  x - x[1])
  fixedCost <- -sum(vapply(profileCost, `[[`, 1, 1) * weight)
  maxScore <- sum(Log2Unrooted(vapply(profileTables, sum, 1)))
  pad <- function (x, len) {
    ret <- double(len)
    ret[seq_along(x)] <- x
    ret
  }
  profiles <- vapply(profileExtra, pad, double(4), 4)
  
  TreeSearch:::morphy_profile(tree$edge, morphyObjects, weight, 
                              charSeq, profiles, Inf)
  
  PP <- function (costs) {
    TreeSearch:::morphy_profile(tree$edge, morphyObjects, weight, 
                                charSeq, costs, Inf)
  }
  
  
  # Use integer-step profile tables
  extraSteps <- matrix(1:4, 4, 4)
  expect_equal(TreeLength(tree, dataset), PP(costs = extraSteps))
  expect_equal(3 + 2 + 2 + 1 + 1,
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
