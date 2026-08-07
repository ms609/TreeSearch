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
    Consistency(StringToPhyDat(char, TipLabels(tree)), tree, nRelabel = 0),
    c(ci = m / s, ri = r, rc = r * m / s, rhi = NA)
  )
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

test_that("ExpectedLength() handles a state only seen within a polymorphism", {
  # A state that never appears on its own -- only ever inside an ambiguous
  # (polymorphic) token -- must not crash .SortTokens()'s wholes/ambiguity
  # remapping.  Regression test for a crash reported downstream of the
  # #88/#87/#94/#112 fix:
  # Error in names(object) <- nm :
  #   'names' attribute [4] must be the same length as the vector [2]
  tree <- TreeTools::BalancedTree(paste0("t", 1:4))
  dat <- StringToPhyDat("00(12)(12)", TipLabels(tree))
  expect_silent(el <- ExpectedLength(dat, tree, nRelabel = 20))
  expect_type(el, "double")
  expect_length(el, 1)
})

test_that(".SortTokens() works", {
  contrast <- structure(c(0, 0, 1, 1, 0, 0, 0, 1, 0,
                          1, 0, 1, 0, 0, 0, 0, 0, 1, 
                          0, 1, 1, 0, 0, 0, 0, 1, 1, 
                          0, 0, 1, 0, 1, 0, 0, 0, 0,
                          0, 0, 1, 0, 0, 1, 0, 0, 0,
                          0, 0, 1, 0, 0, 0, 1, 0, 0), dim = c(9, 6), 
                        dimnames = list(NULL, c("-", "0", "1", "2", "3", "4")))
  cont <- apply(contrast, 1, TreeSearch:::.Bin)
  # Simplest
  expect_equal(TreeSearch:::.SortTokens(rep(1:2, 5:4), 1:2, NA), rep(c(2, 4), 5:4))
  expect_equal(TreeSearch:::.SortTokens(rep(1:2, 4:5), 1:2, NA), rep(c(4, 2), 4:5))
  expect_equal(TreeSearch:::.SortTokens(rep(1:3, 4:6), 1:3, NA), rep(c(4, 2, 6), 4:6))
  
  # Straightforward, no inapp
  expect_equal(TreeSearch:::.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = NA),
               rep(c(2, 8, 4), c(4, 2, 3)))
  
  # Straightforward
  expect_equal(TreeSearch:::.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = 2),
               rep(c(1, 4, 2), c(4, 2, 3)))
  expect_equal(TreeSearch:::.SortTokens(rep(c(1, 2, 4), c(4, 2, 3)), cont, inapp = 4),
               rep(c(2, 1, 4), c(4, 2, 3)))
  
  # Inapplicables with ambiguity
  # TODO it would be nice to return 7 in place of 63, but
  # unnecessarily complex to implement at the moment
  expect_equal(TreeSearch:::.SortTokens(rep(c(1, 2, 3, 4, 8, 9),
                               c(2, 3, 4, 5, 1, 1)), cont, inapp = 1),
               rep(c(4, 2, 63, 1, 3, 6), c(2, 3, 4, 5, 1, 1)))

  # States that only ever occur within a polymorphism (never on their own)
  # must not be silently mapped to zero -- regression test for a crash in
  # ExpectedLength() when a character like "00(12)(12)" is scored.
  expect_equal(TreeSearch:::.SortTokens(rep(1:2, 2:2), c(1, 6), NA),
               rep(c(2, 12), 2:2))
})
