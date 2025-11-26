test_that("All labellings counted", {
  
  labels <- c(2, 3, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  allComb <- combinat::permn(rep(seq_along(labels), labels)) |>
    unique() |>
    unlist() |>
    paste0(collapse = "") |>
    StringToPhyDat(tips = nTip, byTaxon = FALSE)
  expect <- table(CharacterLength(tree, allComb, compress = TRUE)) / 
    factorial(nTip) * prod(factorial(labels))
  FixedTreeCount(tree, c(2, 3, 2), 5) / 24
  expect_equal(FixedTreeCount(tree, c(7), 2),
               c("0" = 1, "1" = 0, "2" = 0))
  expect_equal(FixedTreeCount(tree, c(1, 6), 2),
               c("0" = 0, "1" = 1, "2" = 0))
               
  expect_equal(FixedTreeCount(tree, c(2, 3, 2), 5),
               c("0" = 0, "1" = 0, expect, "5" = 0))
  
})
