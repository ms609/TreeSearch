test_that("All labellings counted", {
  
  labels <- c(2, 3, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  first <- labels[[1]]
  leftAfterFirst <- nTip - first
  second <- labels[[2]]
  nComb <- choose(nTip, first) * choose(leftAfterFirst, second)
  allComb <- apply(combn(nTip, first), 2, function(i) {
    apply(combn(leftAfterFirst, second), 2, function(j) {
      ret <- rep(1, nTip)
      ret2 <- rep(3, leftAfterFirst)
      ret2[j] <- 2
      ret[-i] <- ret2
      ret
    }, simplify = FALSE)
  }) |>
    unlist() |>
    paste0(collapse = "") |>
    StringToPhyDat(tips = nTip, byTaxon = FALSE)
  expect <- table(CharacterLength(tree, allComb, compress = TRUE))
  
  FixedTreeCount(tree, c(2, 3, 2), 5) / 24
  expect_equal(FixedTreeCount(tree, c(7), 2),
               c("0" = 1, "1" = 0, "2" = 0))
  expect_equal(FixedTreeCount(tree, c(1, 6), 2),
               c("0" = 0, "1" = 7, "2" = 0))
               
  expect_equal(FixedTreeCount(tree, c(2, 3, 2), 5),
               c("0" = 0, "1" = 0, expect, "5" = 0))
  
})
