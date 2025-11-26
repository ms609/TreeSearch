test_that("All labellings counted", {
  
  # https://stackoverflow.com/questions/22569176/how-to-generate-permutations-or-combinations-of-object-in-r
  labels <- c(2, 3, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  allComb <- combinat::permn(rep(seq_along(labels), labels)) |>
    unlist() |>
    paste0(collapse = "") |>
    StringToPhyDat(tips = nTip, byTaxon = FALSE)
  expect <- table(CharacterLength(tree, allComb, compress = TRUE))
  
  expect_equal(as.numeric(FixedTreeSteps(tr, 4, c(2, 3, 2))[3:5]),
               as.numeric(expect))
  
  
})