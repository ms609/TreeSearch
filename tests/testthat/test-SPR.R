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


test_that("SPR samples uniformly", {
  set.seed(5)
  t5 <- BalancedTree(5)
  # ub(SPR(t5), times = 1000) = ~160 ms
  expect_gt(chisq.test(
    table(replicate(1000, as.integer(as.TreeNumber(SPR(t5))))),
    p = rep(1, 13) / 13)$p.value, 0.001)
})
