test_that("CI & RI calculated correctly", {
  owch3 <- "0102220333"
  tr3 <- ape::read.tree(text=("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
  if (interactive()) {
    PlotCharacter(tr3, StringToPhyDat(owch3, TipLabels(tr3)))
  }
  .AsPhyDat <- function(string, tree) {
    StringToPhyDat(string, TipLabels(tree))
  }
  expect_equal(
    Consistency(.AsPhyDat(owch3, tr3), tr3),
    c(ci = 0.6, ri = NA, rc = NA)
  )
  expect_equal(
    Consistency(.AsPhyDat(rep(owch3, 2), tr3), tr3),
    rbind(c(ci = 0.6, ri = NA, rc = NA), c(ci = 0.6, ri = NA, rc = NA))
  )
})
