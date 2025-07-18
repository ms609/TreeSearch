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

test_that(".SortTokens() handles edge cases", {
  expect_equal(.SortTokens(
    c(1, 1, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1),
    c(7, 1, 2, 4, 3)),
    c(14, 14, 2, 2, 2, 2, 2, 2, 2, 4, 4, 14)
  )
  
  expect_equal(.SortTokens(
    c(1, 2, 3, 4, 5, 5, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 8, 8, 8),
    c(1, 2, 192, # 192 gives ambiguity between two unobserved tokens
      16, 8, 32, 255, 4, 12, 5)),
    c(8, 16, 192 * 2, 32, 4, 4, 64, rep(510, 7), rep(2, 5))
  )
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
    Consistency(StringToPhyDat(char, TipLabels(tree)), tree, nRelabel = 0),
    c(ci = m / s, ri = r, rc = r * m / s, rhi = NA)
  )
  set.seed(1)
  byChar <- Consistency(StringToPhyDat(char, TipLabels(tree)), tree,
                        nRelabel = 10, byChar = TRUE)
  set.seed(1)
  byTree <- Consistency(StringToPhyDat(char, TipLabels(tree)), tree,
                        nRelabel = 10, byChar = FALSE)
  expect_equal(byChar, byTree)
})

test_that("RHI calculated okay", {
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
  
  null <- 6
  # calculated slightly cheekily using
  # median(replicate(10000,
  #                  TreeLength(RandomTree(tree, root = TRUE), charDat)))
  # RHI uses leaf rearrangement, not randomization
  expect_equal(
    Consistency(StringToPhyDat(char, TipLabels(tree)), tree, nRelabel = 100),
    c(ci = m / s, ri = r, rc = r * m / s, rhi = h / (null - m))
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
  null <- 6
  # calculated slightly cheekily using
  # median(replicate(10000,
  #                  TreeLength(RandomTree(tree, root = TRUE), charDat)))
  # RHI uses leaf rearrangement, not randomization
  
  exp <- c(ci = m / s, ri = r, rc = r * m / s, rhi = h / (null - m))
  expect_equal(
    Consistency(StringToPhyDat(c(char, char), TipLabels(tree)), tree,
                nRelabel = 42),
    rbind(exp, exp, deparse.level = 0)
  )
})

test_that(".SortTokens() works", {
  contrast <- structure(c(0, 0, 1, 1, 0, 0, 0, 1, 0,
                          1, 0, 1, 0, 0, 0, 0, 0, 1, 
                          0, 1, 1, 0, 0, 0, 0, 1, 1, 
                          0, 0, 1, 0, 1, 0, 0, 0, 0,
                          0, 0, 1, 0, 0, 1, 0, 0, 0,
                          0, 0, 1, 0, 0, 0, 1, 0, 0), dim = c(9, 6), 
                        dimnames = list(NULL, c("-", "0", "1", "2", "3", "4")))
  cont <- apply(contrast, 1, .Bin)
  # Simplest
  expect_equal(.SortTokens(rep(1:2, 5:4), 1:2, NA), rep(c(2, 4), 5:4))
  expect_equal(.SortTokens(rep(1:2, 4:5), 1:2, NA), rep(c(4, 2), 4:5))
  expect_equal(.SortTokens(rep(1:3, 4:6), 1:3, NA), rep(c(4, 2, 6), 4:6))
  
  # Straightforward, no inapp
  expect_equal(.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = NA),
               rep(c(2, 8, 4), c(4, 2, 3)))
  
  # Straightforward
  expect_equal(.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = 2),
               rep(c(1, 4, 2), c(4, 2, 3)))
  expect_equal(.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = 4),
               rep(c(2, 1, 4), c(4, 2, 3)))
  
  # Inapplicables with ambiguity
  # TODO it would be nice to return 7 in place of 63, but
  # unnecessarily complex to implement at the moment
  expect_equal(.SortTokens(rep(c(1, 2, 3, 4, 8, 9), 
                               c(2, 3, 4, 5, 1, 1)), cont, inapp = 1),
               rep(c(4, 2, 63, 1, 3, 6), c(2, 3, 4, 5, 1, 1)))
  
})
