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

test_that(".SortChar() works", {
  contrast <- structure(c(0, 0, 1, 1, 0, 0, 0,
                          1, 0, 1, 0, 0, 0, 0,
                          0, 1, 1, 0, 0, 0, 0,
                          0, 0, 1, 0, 1, 0, 0,
                          0, 0, 1, 0, 0, 1, 0,
                          0, 0, 1, 0, 0, 0, 1), dim = 7:6, 
                        dimnames = list(NULL, c("-", "0", "1", "2", "3", "4")))
  cont <- apply(contrast, 1, .Bin)
  # Simplest
  expect_equal(.SortChar(rep(1:2, 5:4), 1:2, NA), rep(1:2, 5:4))
  expect_equal(.SortChar(rep(1:2, 4:5), 1:2, NA), rep(2:1, 4:5))
  expect_equal(.SortChar(rep(1:3, 4:6), 1:2, NA), rep(c(2, 1, 3), 4:6))
  
  # Straightforward, no inapp
  expect_equal(.SortChar(rep(2 ^ c(1, 2, 4), c(4, 2, 3)), cont, inapp = NA),
               rep(2 ^ c(0, 2, 1), c(4, 2, 3)))
  
  # Straightforward
  expect_equal(.SortChar(rep(2^c(1, 2, 4), c(4, 2, 3)), cont, inapp = 1),
               rep(2 ^ c(1, 3, 2), c(4, 2, 3)))
  
  # Inapplicable
  expect_equal(.SortChar(rep(2^c(1, 2, 0, 4), c(4, 2, 3, 3)), cont, inapp = 1),
               rep(2^c(1, 3, 0, 2), c(4, 2, 3, 3)))
  
  # Ambiguity: 8->2; 2->4; 1->1
  expect_equal(.SortChar(c(8, 8, 8 + 2, 8 + 1, 8, 2 + 1, 2, 2, 1, 1, 1), cont, inapp = 1),
               c(2, 2, 2 + 4, 2 + 1, 2, 4 + 1, 4, 4, 1, 1, 1))
})
