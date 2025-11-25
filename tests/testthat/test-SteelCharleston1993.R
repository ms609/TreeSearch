test_that("SteelCharleston results stack up", {
  expect_warning(expect_equal(SteCha1("10", 4), 0), "numeric")
  raw <- SteCha1(10, 4)
  expect_equal(exp(SteCha1(10, 4, TRUE)), raw)
  expect_equal(2 ^ (SteCha1(TreeTools::StarTree(10), 4, 2)), raw)
  
  nTip <- 8
  char <- R.utils::intToBin(0:(2 ^ nTip - 1)) |>
    as.numeric() |>
    formatC(width = nTip, flag = "0", format = "d") |>
    paste0(collapse = "")
  expected <- tabulate(
    1 + CharacterLength(TreeTools::BalancedTree(nTip),
                        StringToPhyDat(char, nTip, byTaxon = FALSE))
  )
  
  for (i in seq_along(expected)) {
    expect_equal(SteCha1(nTip, i - 1), expected[[i]])
  }
})
