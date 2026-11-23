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
    rbind(c(ci = m / s, ri = r, rc = r * m / s, rhi = NA), deparse.level = 0)
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
    rbind(c(ci = m / s, ri = r, rc = r * m / s, rhi = h / (null - m)),
          deparse.level = 0)
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
})

test_that(".SortTokens() keeps a present-only partial-ambiguity token", {
  # "-" = 1, "0" = 2, "1" = 4, "2" = 8, "?" = 15 (fully ambiguous),
  # "(01)" = 6 (ambiguous over states 0 and 1 only)
  contr <- c(1, 2, 4, 8, 15, 6)
  # Character uses only "-", "0", "1" and "(01)" -- token 5 ("?") never
  # appears, so the dataset-wide ambiguous set {15, 6} is broader than the
  # ambiguous tokens actually present in this character ({6})
  char <- rep(c(2, 3, 6), c(3, 3, 3))

  # "(01)" must be rewritten to the union of its own two present states'
  # new codes (0 -> 2, 1 -> 4; union = 6), not corrupted by "?"
  expect_equal(TreeSearch:::.SortTokens(char, contr, inapp = 1),
               rep(c(2, 4, 6), c(3, 3, 3)))

  # A second contrast in which the only ambiguous token present ("?", fully
  # ambiguous) sits at a different position from the unused ambiguous token
  # ("(01)", contr[5])
  expect_equal(
    TreeSearch:::.SortTokens(c(1, 1, 3, 3, 3, 3, 3, 3, 3, 2, 2, 1),
                              c(7, 1, 2, 4, 3)),
    c(14, 14, 2, 2, 2, 2, 2, 2, 2, 4, 4, 14)
  )
})

test_that("ExpectedLength() cache does not collide across tree topologies", {
  tips <- paste0("t", 1:16)
  bal <- TreeTools::BalancedTree(tips)
  pec <- TreeTools::PectinateTree(tips)
  charDat <- StringToPhyDat("0000000000011111", tips)

  set.seed(999)
  balLength <- ExpectedLength(charDat, bal, 500)
  # Scoring `bal` first populates .CharLengthCache; the pectinate query below
  # must not silently reuse `bal`'s cache entry
  set.seed(999)
  pecLength <- ExpectedLength(charDat, pec, 500)

  expect_equal(balLength, 4)
  expect_equal(pecLength, 5)

  # Re-querying `bal` (without resetting the seed) must still return its own
  # cached value, confirming the cache is actually being hit and not merely
  # avoiding collisions by chance
  expect_equal(ExpectedLength(charDat, bal, 500), balLength)
})

test_that("ExpectedLength() cache key is invariant to edge order and labels", {
  tips <- paste0("t", 1:10)
  tree <- TreeTools::BalancedTree(tips)
  charDat <- StringToPhyDat("0000011111", tips)

  set.seed(101)
  postLength <- ExpectedLength(charDat, TreeTools::Postorder(tree), 200)
  nKeys <- length(ls(TreeSearch:::.CharLengthCache))

  # Each of these presents the same rooted shape, so each must reuse the
  # existing entry rather than add one.  Asserting the key count, rather than
  # the returned value, is what makes this a regression test: the median is
  # stable enough that a recomputation would return the same number.
  expect_equal(ExpectedLength(charDat, TreeTools::Preorder(tree), 200),
               postLength)
  expect_equal(ExpectedLength(charDat, ape::rotate(tree, length(tips) + 2L),
                              200), postLength)
  # A different labelling of the same shape samples the same distribution, so
  # it shares the entry too
  relabelled <- TreeTools::RenumberTips(
    TreeTools::BalancedTree(sample(tips)), tips)
  expect_equal(ExpectedLength(charDat, relabelled, 200), postLength)

  expect_equal(length(ls(TreeSearch:::.CharLengthCache)), nKeys)
})

test_that(".ShapeKey() identifies rooted shapes", {
  tips <- paste0("t", 1:12)
  bal <- TreeTools::BalancedTree(tips)
  pec <- TreeTools::PectinateTree(tips)

  # Invariant to edge order, node rotation and labelling; distinguishes shape
  expect_equal(TreeSearch:::.ShapeKey(TreeTools::Preorder(bal)),
               TreeSearch:::.ShapeKey(TreeTools::Postorder(bal)))
  expect_equal(TreeSearch:::.ShapeKey(ape::rotate(bal, 15L)),
               TreeSearch:::.ShapeKey(bal))
  expect_equal(TreeSearch:::.ShapeKey(TreeTools::RenumberTips(
                 TreeTools::BalancedTree(rev(tips)), tips)),
               TreeSearch:::.ShapeKey(bal))
  expect_false(TreeSearch:::.ShapeKey(pec) == TreeSearch:::.ShapeKey(bal))

  # Agrees with TreeTools' independent enumeration as an equivalence relation
  set.seed(2)
  trees <- lapply(1:60, function(i) TreeTools::RandomTree(9, root = TRUE))
  mine <- vapply(trees, TreeSearch:::.ShapeKey, character(1))
  theirs <- vapply(trees, function(tr) {
    as.character(TreeTools::RootedTreeShape(tr))
  }, character(1))
  expect_equal(as.integer(factor(mine, levels = unique(mine))),
               as.integer(factor(theirs, levels = unique(theirs))))

  # Unlike RootedTreeShape(), no leaf-count ceiling
  expect_error(TreeTools::RootedTreeShape(TreeTools::BalancedTree(56)))
  expect_type(TreeSearch:::.ShapeKey(TreeTools::BalancedTree(56)), "character")
  # Leaf counts whose codes pad to the same length stay distinct
  expect_false(TreeSearch:::.ShapeKey(TreeTools::BalancedTree(55)) ==
                 TreeSearch:::.ShapeKey(TreeTools::BalancedTree(56)))
})

test_that("Consistency() returns a matrix, not a vector, for one character", {
  tree <- ape::read.tree(
    text = ("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
  charDat <- StringToPhyDat("0102220333", TipLabels(tree))

  res <- Consistency(charDat, tree, nRelabel = 0)
  expect_true(is.matrix(res))
  expect_equal(dim(res), c(1, 4))
  expect_equal(colnames(res), c("ci", "ri", "rc", "rhi"))
})

test_that("Consistency() returns the documented NaN for degenerate chars", {
  tree <- TreeTools::BalancedTree(6)
  tips <- TipLabels(tree)

  # Constant character: no informative variation, so observed and minimum
  # length are both zero -> ci is 0/0
  constDat <- StringToPhyDat("000000", tips)
  constRes <- Consistency(constDat, tree, nRelabel = 0)
  expect_true(is.nan(constRes[, "ci"]))
  expect_true(is.nan(constRes[, "ri"]))
  expect_true(is.nan(constRes[, "rc"]))

  # Autapomorphy: a single tip differs, so maximum and minimum length
  # coincide -> ri and rc are 0/0
  autDat <- StringToPhyDat("000001", tips)
  autRes <- Consistency(autDat, tree, nRelabel = 0)
  expect_true(is.nan(autRes[, "ri"]))
  expect_true(is.nan(autRes[, "rc"]))

  # Median null length equals the minimum length -> rhi is 0/0
  set.seed(1)
  rhiRes <- Consistency(autDat, tree, nRelabel = 50)
  expect_true(is.nan(rhiRes[, "rhi"]))

  # States that only ever occur within a polymorphism (never on their own)
  # must not be silently mapped to zero -- regression test for a crash in
  # ExpectedLength() when a character like "00(12)(12)" is scored.
  expect_equal(TreeSearch:::.SortTokens(rep(1:2, 2:2), c(1, 6), NA),
               rep(c(2, 12), 2:2))
})
