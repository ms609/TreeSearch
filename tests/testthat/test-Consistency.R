test_that("Consistency() fails gracefully with unrooted trees", {
  tree <- TreeTools::RandomTree(8, root = FALSE)
  char <- "00112222"
  expect_error(Consistency(StringToPhyDat(char, TipLabels(tree)), tree),
               "tree. must be rooted")
})

test_that("Consistency() notes tree-leaf mismatch", {
  tree <- TreeTools::BalancedTree(10)
  char <- "00112222"
  expect_error(Consistency(StringToPhyDat(char, TipLabels(tree)[-c(1:2)]), tree),
               "Tip label mismatch")
})

test_that("CI & RI calculated correctly", {
  tree <- ape::read.tree(
    text = ("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
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

test_that("Consistency() handles `-`", {
  tree <- ape::read.tree(
    text = ("((a1, a2), (((b1, b2), (i1, (i2, (i3, (c, d))))), ((e1, e2), (f, g))));"))
  char <- "0102---220333"
  charDat <- StringToPhyDat(char, TipLabels(tree))
  if (interactive()) {
    PlotCharacter(tree, charDat)
  }
  m <- 3
  expect_equal(MinimumLength(char, tree), m)
  s <- 5
  expect_equal(TreeLength(tree, charDat), s)
  h <- s - m
  g <- 7 + 1
  expect_equal(MaximumLength(char, tree), g)
  r <- (g - s) / (g - m)
  expect_equal(
    Consistency(StringToPhyDat(c(char, char), TipLabels(tree)), tree),
    rbind(c(ci = m / s, ri = r, rc = r * m / s),
          c(ci = m / s, ri = r, rc = r * m / s))
  )
})
