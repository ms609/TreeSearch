library("TreeTools")

test_that("TBR errors", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, -1))))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, 1))))
  expect_equal(0, length(expect_warning(all_tbr(tr$edge, 111))))
})

test_that("SPR errors", {
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

test_that("SPR working", {
  tr <- Preorder(root(TreeTools::BalancedTree(7), 't1', resolve.root = TRUE))

  # Move single tip
  expect_equal(8, length(x <- all_spr(tr$edge, 12)))
  expect_equal(8, length(x <- all_spr(tr$edge, 11)))
  expect_equal(8, length(x <- all_spr(tr$edge, 10)))
  expect_equal(8, length(x <- all_spr(tr$edge, 7)))
  expect_equal(8, length(x <- all_spr(tr$edge, 6)))
  expect_equal(8, length(x <- all_spr(tr$edge, 3)))
  
  # Move cherry
  expect_equal(6, length(x <- all_spr(tr$edge, 9)))
  expect_equal(6, length(x <- all_spr(tr$edge, 5)))
  
  # Move more
  expect_equal(0, length(unique(x <- all_spr(tr$edge, 4))))
  expect_equal(4, length(unique(x <- all_spr(tr$edge, 8))))
  
  # All moves
  expect_equal(6*8 + 2*6 + 4, length(x <- all_spr(tr$edge, integer(0))))
  expect_equal(48, length(unique(x <- all_spr(tr$edge, integer(0))))) # 48 not formally calculated
  expect_equal(48, length(SPRMoves(tr)))
  
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