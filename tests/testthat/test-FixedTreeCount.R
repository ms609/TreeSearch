test_that("All labellings counted", {
  
  # https://stackoverflow.com/questions/22569176/how-to-generate-permutations-or-combinations-of-object-in-r
  labels <- c(2, 3, 2)
  nTip <- sum(labels)
  tree <- TreeTools::BalancedTree(nTip)
  allComb <- combinat::permn(rep(seq_along(labels), labels)) |>
    unlist() |>
    paste0(collapse = "") |>
    StringToPhyDat(tips = nTip, byTaxon = FALSE)
  
  
})