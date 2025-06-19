library("TreeTools")

test_that("SPR handles root rearrangement", {
  # Soft tests.
  expect_true(inherits(SPR(BalancedTree(4), edgeToBreak = 1), "phylo"))
  expect_true(inherits(SPR(BalancedTree(5), edgeToBreak = 1), "phylo"))
  
  expect_equal(sum(sapply(c(
    # hard-coded 1 %in% 2 due to 
    # https://github.com/r-lib/testthat/issues/1661
    SPR(PectinateTree(4), edgeToBreak = 1, mergeEdge = 3),
    SPR(PectinateTree(4), edgeToBreak = 1, mergeEdge = 4),
    SPR(PectinateTree(4), edgeToBreak = 1, mergeEdge = 5),
    SPR(PectinateTree(4), edgeToBreak = 1, mergeEdge = 6)),
    all.equal, SPR(PectinateTree(4), edgeToBreak = 1))), 1)
  expect_warning(SPR(PectinateTree(4), edgeToBreak = 2),
                 "No rearrangement possible with this root position")
})
