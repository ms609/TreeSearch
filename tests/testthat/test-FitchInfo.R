test_that("Downpass handles overlap vs union correctly", {
  dp <- .DownpassOutcome(4)
  
  # {A}, {A,B} -> {A}
  expect_equal(dp[1, 3], as.raw(1)) 
  
  # {A}, {B} -> {A,B}
  expect_equal(dp[1, 2], as.raw(3)) 
  
  # {A,B}, {A,B,C} -> {A,B}
  expect_equal(dp[3, 7], as.raw(3)) 
  
  # symmetry
  expect_equal(dp[2, 1], dp[1, 2])
})

test_that("Uppass covers all cases", {
  up <- .UppassOutcome(4)
  
  # Case I: ancestor subset of node
  expect_equal(up[1, 3, 1], as.raw(1)) # {A}, {A,B}, ancestor {A}
  
  # Case III: disjoint descendants
  expect_equal(up[1, 2, 4], as.raw(7)) # {A}, {B}, ancestor {C}
  
  # Default case
  expect_equal(up[3, 6, 4], as.raw(6)) # {A,B}, {B,C}, ancestor {C} => {B,C}
  
  # Ancestor irrelevant
  expect_equal(up[1, 1, 8], as.raw(1)) # {A}, {A}, ancestor {D}
})

test_that("Multi-state characters calculated correctly", {
  library("TreeTools", quietly = TRUE)
  tree <- BalancedTree(6)
  character <- MatrixToPhyDat(`rownames<-`(rbind(1, 1, 2, 2, 3, 3), TipLabels(6)))
  all6 <- as.phylo(0:104, 6)
  tl1 <- table(vapply(all6, TreeLength, 1, character))
  fi <- table(vapply(all6, FitchInfo, 0, character))
  expect_equal(unname(fi), unname(rev(tl1)))
  expect_equal(as.numeric(names(fi)),
               unname(rev(-log2(cumsum(tl1 / NUnrooted(6))))))
  
})
