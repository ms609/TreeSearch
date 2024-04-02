test_that("CI & RI calculated correctly", {
  tree <- ape::read.tree(text=("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
  char <- "0102220333"
  charDat <- StringToPhyDat(char, TipLabels(tree))
  if (interactive()) {
    PlotCharacter(tree, charDat)
  }
  m <- 3
  expect_equal(MinimumLength(char, tree), m)
  s <- 5
  expect_equal(TreeLength(tree, charDat), s)
  h <- s - m
  g <- 7
  expect_equal(MaximumLength(char, tree), g)
  r <- (g - s) / (g - m)
  expect_equal(
    Consistency(StringToPhyDat(char, TipLabels(tree)), tree),
    c(ci = m / s, ri = r, rc = r * m / s)
  )
})
