library("TreeTools")

test_that("TBR errors", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, -1))))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, 1))))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, 111))))
})

test_that("SPR errors", {
  skip_if(TRUE)
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  expect_equal(0, length(expect_warning(all_spr(tr$edge, -1))))
  expect_equal(0, length(expect_warning(all_spr(tr$edge, 1))))
  expect_equal(0, length(expect_warning(all_spr(tr$edge, 111))))
})

test_that("TBR working", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))

  # Move single tip
  expect_equal(8, length(x <- all_tbr(tr$edge, 12)))
  expect_equal(8, length(x <- all_tbr(tr$edge, 11)))
  expect_equal(8, length(x <- all_tbr(tr$edge, 10)))
  expect_equal(8, length(x <- all_tbr(tr$edge, 7)))
  expect_equal(8, length(x <- all_tbr(tr$edge, 6)))
  expect_equal(8, length(x <- all_tbr(tr$edge, 3)))
  
  # Move cherry
  expect_equal(6, length(x <- all_tbr(tr$edge, 9)))
  expect_equal(6, length(x <- all_tbr(tr$edge, 5)))
  expect_equal(6, length(TBRMoves(tr, 5)))
  
  # Move more
  expect_equal(6, length(unique(x <- all_tbr(tr$edge, 4))))
  expect_equal(3 * 4 + 2, length(unique(x <- all_tbr(tr$edge, 8))))
  
  # All moves
  expect_equal(6*8 + 12+ 6 + 14, length(x <- all_tbr(tr$edge, integer(0))))
  expect_equal(58, length(unique(x <- all_tbr(tr$edge, integer(0))))) # 58 not formally calculated
  expect_equal(58, length(TBRMoves(tr)))
  
  tr <- Preorder(root(TreeTools::BalancedTree(14), 't1', resolve.root = TRUE))
  desc <- TreeTools::CladeSizes(tr)
  
  external <- c(3, 6, 7, 11, 12, 13, 17, 18, 20, 21, 24:26)
  # Move single
  for (leaf in external) {
    expect_equal(22, length(x <- all_tbr(tr$edge, leaf)))
  }
  
  Test <- function (edge) {
    nDesc <- desc[tr$edge[edge, 2]]
    expected <- (2 * nDesc - 3) * (22 - (2 * nDesc - 3)) - 1
    expect_equal(expected, length(all_tbr(tr$edge, edge)))
  }
  for (internal in which(!1:26 %in% external)[-(1:2)]) {
    Test(internal)
  }
})

test_that("SPR fails gracefully", {
  expect_error(all_spr(as.phylo(1, 3)$edge, integer(0)))
  expect_error(all_spr(Postorder(as.phylo(1, 6))$edge, integer(0)))
  expect_error(all_spr(SortTree(as.phylo(1, 6))$edge, integer(0)))
})

test_that("SPR works", {
  dput(" - SPR1 ")
  t2 <- as.phylo(518, 7) # (t1, ((t2, t3), ((t4, t5), (t6, t7))))
  expect_equal(8, length(all_spr(t2$edge, 2)))
  
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))

  dput(" - SPR Single tip ")
  # Move single tip
  expect_equal(8, length(all_spr(tr$edge, 12)))
  expect_equal(8, length(all_spr(tr$edge, 11)))
  expect_equal(8, length(all_spr(tr$edge, 10)))
  expect_equal(8, length(all_spr(tr$edge, 7)))
  expect_equal(8, length(all_spr(tr$edge, 6)))
  expect_equal(8, length(all_spr(tr$edge, 3)))
  expect_equal(8, length(all_spr(tr$edge, 2)))
  
  dput(" - SPR Cherry ")
  # Move cherry
  expect_equal(6, length(all_spr(tr$edge, 9)))
  expect_equal(6, length(all_spr(tr$edge, 5)))
  expect_equal(12, length(all_spr(tr$edge, c(9, 5))))
  
  dput(" - SPR Bush ")
  # Move more
  expect_equal(0, length(unique(all_spr(tr$edge, 4))))
  expect_equal(4, length(unique(all_spr(tr$edge, 8))))
  
  dput(" - SPR All ")
  # All moves
  expect_equal(7*8 + 2*6 + 4, length(all_spr(tr$edge, integer(0))))
  uniqueMoves <- length(unique(all_spr(tr$edge, integer(0))))
  expect_equal(54, # Not formally calculated
               uniqueMoves)
  expect_equal(uniqueMoves, length(SPRMoves(tr)))
  
  dput(" - SPR Clear ")
  tr <- Preorder(root(TreeTools::BalancedTree(14), 't1', resolve.root = TRUE))
  tr$edge
  desc <- TreeTools::CladeSizes(tr)
  
  external <- c(3, 6, 7, 11, 12, 13, 17, 18, 20, 21, 24:26)
  # Move single
  for (leaf in external) {
    expect_equal(22, length(x <- all_spr(tr$edge, leaf)))
  }
  
  Test <- function (edge) {
    nDesc <- desc[tr$edge[edge, 2]]
    expected <- (22 - (2 * nDesc - 3)) - 1
    expect_equal(expected, length(all_spr(tr$edge, edge)))
  }
  for (internal in which(!1:26 %in% external)[-(1:2)]) {
    Test(internal)
  }
  
  expect_equal(SPRMoves(tr)[[428]]$edge, SPRMoves(tr$edge)[[428]])
  
  tr <- BalancedTree(7)
  expect_equal(SPRMoves(tr)[[54]]$edge, SPRMoves(tr$edge)[[54]])
})

# TODO Restore or delete
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