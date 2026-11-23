library("TreeTools")

test_that("TBR errors", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  expect_warning(r1 <- TreeSearch:::all_tbr(tr$edge, -1));   expect_equal(0, length(r1))
  expect_warning(r2 <- TreeSearch:::all_tbr(tr$edge, 1));    expect_equal(0, length(r2))
  expect_warning(r3 <- TreeSearch:::all_tbr(tr$edge, 111));  expect_equal(0, length(r3))
})

test_that("SPR errors", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  expect_warning(r1 <- TreeSearch:::all_spr(tr$edge, -1));   expect_equal(0, length(r1))
  expect_warning(r2 <- TreeSearch:::all_spr(tr$edge, 1));    expect_equal(0, length(r2))
  expect_warning(r3 <- TreeSearch:::all_spr(tr$edge, 111));  expect_equal(0, length(r3))
})

test_that("TBR working", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))

  # Move single tip
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 12)))
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 11)))
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 10)))
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 7)))
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 6)))
  expect_equal(8, length(x <- TreeSearch:::all_tbr(tr$edge, 3)))

  # Move tip 1, by bisecting the root edge.  Edges 1 and 2 are the two halves
  # of tip 1's pendant edge, so this is a bisection like any other: the
  # fragment has n - 1 = 6 leaves and hence 2 * 6 - 3 = 9 places to re-root,
  # one of which recreates `tr`.
  expect_equal(2 * (7 - 1) - 3 - 1, length(TreeSearch:::all_tbr(tr$edge, 2)))
  # The severed tip 1 admits no re-rooting, so TBR and SPR coincide here.
  expect_identical(TreeSearch:::all_tbr(tr$edge, 2),
                   TreeSearch:::.all_spr(tr$edge, 2))

  # Move cherry
  expect_equal(6, length(x <- TreeSearch:::all_tbr(tr$edge, 9)))
  expect_equal(6, length(x <- TreeSearch:::all_tbr(tr$edge, 5)))
  expect_equal(6, length(TBRMoves(tr, 5)))

  # Move more
  expect_equal(6, length(unique(x <- TreeSearch:::all_tbr(tr$edge, 4))))
  expect_equal(3 * 4 + 2, length(unique(x <- TreeSearch:::all_tbr(tr$edge, 8))))

  # All moves: seven single-leaf bisections (six pendant edges plus the root
  # edge) at eight apiece, then the two cherries, the third-from-root edge,
  # and the deepest internal edge.
  expect_equal(7*8 + 12+ 6 + 14, length(x <- TreeSearch:::all_tbr(tr$edge, integer(0))))
  # 64 = size of the complete TBR neighbourhood of this tree, from the
  # independent oracle in dev/tbr-root-edge/, which agrees with an exhaustive
  # TBRSwap() sweep and with every unrooted seven-leaf topology.
  expect_equal(64, length(unique(x <- TreeSearch:::all_tbr(tr$edge, integer(0)))))
  expect_equal(64, length(TBRMoves(tr)))

  # TBR contains SPR by definition; omitting a break edge from either
  # enumerator breaks this (agent-issues/TreeSearch#147).  Key on split
  # membership by tip label, which is injective over unrooted tree space and
  # blind to where a tree happens to be rooted.
  Key <- function (trees) {
    unique(vapply(trees, function (tree) {
      splits <- TreeTools::as.Splits(tree)
      members <- as.logical(splits)
      if (is.null(dim(members))) members <- matrix(members, nrow = 1)
      colnames(members) <- attr(splits, "tip.label")
      members <- members[, order(colnames(members)), drop = FALSE]
      tips <- colnames(members)
      paste(sort(apply(members, 1, function (inSplit) {
        # Complement so that the first tip is always outside the split
        if (inSplit[[1]]) inSplit <- !inSplit
        paste0(tips[inSplit], collapse = ",")
      })), collapse = "|")
    }, character(1)))
  }
  expect_true(all(Key(SPRMoves(tr)) %in% Key(TBRMoves(tr))))

  tr <- Preorder(root(TreeTools::BalancedTree(14), 't1', resolve.root = TRUE))
  desc <- TreeTools::CladeSizes(tr)

  external <- c(3, 6, 7, 11, 12, 13, 17, 18, 20, 21, 24:26)
  # Move single
  for (leaf in external) {
    expect_equal(22, length(x <- TreeSearch:::all_tbr(tr$edge, leaf)))
  }
  # Moving tip 1 by bisecting the root edge costs the same, for the same
  # reason: 2 * 13 - 3 re-rootings of the 13-leaf fragment, less the identity.
  expect_equal(2 * (14 - 1) - 3 - 1, length(TreeSearch:::all_tbr(tr$edge, 2)))
  expect_identical(TreeSearch:::all_tbr(tr$edge, 2),
                   TreeSearch:::.all_spr(tr$edge, 2))

  Test <- function (edge) {
    nDesc <- desc[tr$edge[edge, 2]]
    expected <- (2 * nDesc - 3) * (22 - (2 * nDesc - 3)) - 1
    expect_equal(expected, length(TreeSearch:::all_tbr(tr$edge, edge)))
  }
  for (internal in which(!1:26 %in% external)[-(1:2)]) {
    Test(internal)
  }
})

test_that("SPR fails gracefully", {
  # `.all_spr` guards the conditions that ASAN dislikes seeing Rcpp::stop() on.
  expect_error(TreeSearch:::.all_spr(as.phylo(1, 3)$edge, integer(0)),
               "< 5 edges")
  expect_error(TreeSearch:::.all_spr(Postorder(as.phylo(1, 6))$edge, integer(0)),
               "must connect root to leaf")
  expect_error(TreeSearch:::.all_spr(SortTree(as.phylo(1, 6))$edge, integer(0)),
               "must connect root to leaf")
})

test_that("SPR works", {
  t2 <- as.phylo(518, 7) # (t1, ((t2, t3), ((t4, t5), (t6, t7))))
  expect_equal(8, length(TreeSearch:::all_spr(t2$edge, 2)))
  
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))

  # Move single tip
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 12)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 11)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 10)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 7)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 6)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 3)))
  expect_equal(8, length(TreeSearch:::all_spr(tr$edge, 2)))
  
  # Move cherry
  expect_equal(6, length(TreeSearch:::all_spr(tr$edge, 9)))
  expect_equal(6, length(TreeSearch:::all_spr(tr$edge, 5)))
  expect_equal(12, length(TreeSearch:::all_spr(tr$edge, c(9, 5))))
  
  # Move more
  expect_equal(0, length(unique(TreeSearch:::all_spr(tr$edge, 4))))
  expect_equal(4, length(unique(TreeSearch:::all_spr(tr$edge, 8))))
  
  # All moves
  expect_equal(7*8 + 2*6 + 4, length(TreeSearch:::all_spr(tr$edge, integer(0))))
  uniqueMoves <- length(unique(TreeSearch:::all_spr(tr$edge, integer(0))))
  expect_equal(54, # Not formally calculated
               uniqueMoves)
  expect_equal(uniqueMoves, length(SPRMoves(tr)))
  
  tr <- Preorder(root(TreeTools::BalancedTree(14), 't1', resolve.root = TRUE))
  tr$edge
  desc <- TreeTools::CladeSizes(tr)
  
  external <- c(3, 6, 7, 11, 12, 13, 17, 18, 20, 21, 24:26)
  # Move single
  for (leaf in external) {
    expect_equal(22, length(x <- TreeSearch:::all_spr(tr$edge, leaf)))
  }
  
  Test <- function (edge) {
    nDesc <- desc[tr$edge[edge, 2]]
    expected <- (22 - (2 * nDesc - 3)) - 1
    expect_equal(expected, length(TreeSearch:::all_spr(tr$edge, edge)))
  }
  for (internal in which(!1:26 %in% external)[-(1:2)]) {
    Test(internal)
  }
  
  expect_equal(SPRMoves(tr)[[428]]$edge, SPRMoves(tr$edge)[[428]])
  
  tr <- BalancedTree(7)
  expect_equal(SPRMoves(tr)[[54]]$edge, SPRMoves(tr$edge)[[54]])
})

if (FALSE) test_that("SPR works", {
  testTree <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  plot(testTree); nodelabels(); edgelabels()
  edge <- testTree$edge
  
  t2 <- testTree
  #t2$edge = root_on_node(edge, 11)
  plot(t2)
  
  1L + tbr_moves(edge)
  
  Test <- function (m, p1, r1) {
    test.tr <- testTree
    test.tr$edge <- spr(edge, m)
    plot(test.tr)
    
    oldWay <- SortTree(root(SPR(testTree, p1, r1), 't1', resolve.root = TRUE))
    expect_equal(oldWay, SortTree(test.tr))
  }
  Test(0, 1, 5)
  Test(1, 1, 6)
  Test(2, 1, 7)
  Test(3, 1, 8)
  Test(4, 1, 9)
  Test(5, 1, 10)
  Test(6, 1, 11)
  Test(7, 1, 12)
  
  Test(8 , 3, 5)
  Test(9 , 3, 6)
  Test(10, 3, 7)
  Test(11, 3, 8)
  Test(12, 3, 9)
  Test(13, 3, 10)
  Test(14, 3, 11)
  Test(15, 3, 12)
  

})