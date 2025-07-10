library("TreeTools")

context("Tree rearrangements")
tree5a <- read.tree(text = '(a, (b, (c, (d, e))));')
tree5b <- read.tree(text = '((a, b), (c, (d, e)));')
tree6  <- Preorder(read.tree(text = "((a, (b, (c, d))), (e, f));"))
tree6b <- Preorder(read.tree(text = "((a, (b, c)), (d, (e, f)));"))
tree8  <- read.tree(text = "(((a, (b, (c, d))), (e, f)), (g, h));")
tree11 <- read.tree(text = "((((a, b), (c, d)), e), ((f, (g, (h, i))), (j, k)));")
attr(tree5a, 'order') <- attr(tree5b, 'order') <- attr(tree8, 'order') <- attr(tree11, 'order') <- 'preorder'

test_that("Malformed trees don't crash anything", {
  treeDoubleNode <- read.tree(text = "((((((1,2)),3),4),5),6);")
  treePolytomy   <- read.tree(text = "((((1,2,3),4),5),6);")
  treeDoublyPoly <- read.tree(text = "(((((1,2,3)),4),5),6);")

  expect_error(NNI(treeDoubleNode))
  expect_error(NNI(treePolytomy))
  expect_error(NNI(treeDoublyPoly))
  
  expect_error(SPR(treeDoubleNode))
  expect_error(SPR(treePolytomy))
  expect_error(SPR(treeDoublyPoly))
  
  expect_error(TBR(treeDoubleNode))
  expect_error(TBR(treePolytomy))
  expect_error(TBR(treeDoublyPoly))
  
})

test_that("NNI works", {
  trComb <- read.tree(text = "(((((1,2),3),4),5),6);")
  edge <- trComb$edge
  Test <- function (e, r, e1, e2) {
    edge1 <- edge
    edge1[c(e1, e2), 2] <- edge1[c(e2, e1), 2]
    edge1 <- do.call(cbind, RenumberEdges(edge1[, 1], edge1[, 2]))
    expect_equal(edge1, nni(trComb$edge, e, r))
  }
  Test(0, 0, 5, 7)
  Test(0, 2, 5, 7)
  Test(3, 0, 5, 7) # Option 0 == option 3.
  Test(0, 1, 6, 7)
  Test(1, 0, 4, 8)
  Test(1, 1, 7, 8)
  Test(2, 0, 3, 9)
  Test(2, 1, 8, 9)
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(0)
  nniComb <- NNI(trComb)
  expect_equal(nniComb$tip.label, trComb$tip.label)
  expect_equal(nniComb$Nnode, trComb$Nnode)
  expect_equal(nniComb, read.tree(text = "(((((3,2),1),4),5),6);"))  
})


test_that("SPR works", {
  testTree <- Preorder(root(BalancedTree(7), 1, resolve.root = TRUE))
  edge <- testTree[["edge"]]
  expect_equal(spr(edge, 66), cSPR(testTree, 66)$edge)
  
  Test <- function (m, p1, r1) {
    test.tr <- testTree
    test.tr$edge <- spr(edge, m)
    
    oldWay <- SortTree(root(SPR(testTree, p1, r1), "t1", resolve.root = TRUE))
    expect_equal(oldWay, SortTree(test.tr))
  }
  Test(0, 1, 5)
  Test(64, 1, 5) # Modulo 64!
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
  
  Test(16, 5, 3)
  Test(17, 5, 9)
  Test(18, 5, 10)
  Test(19, 5, 11)
  Test(20, 5, 12)
  
  Test(28, 7, 3)
  Test(29, 7, 4)
  Test(30, 7, 8)
  Test(31, 7, 9)
  Test(32, 7, 10)
  Test(33, 7, 11)
  Test(34, 7, 12)
  
  Test(35,  8,  3)
  Test(36,  8,  6)
  Test(37,  8,  7)
  Test(38,  9,  3)
  Test(39,  9,  4)
  Test(40,  9,  5)
  Test(41,  9,  6)
  Test(42,  9,  7)
  Test(43, 10,  3)
  Test(44, 10,  4)
  Test(45, 10,  5)
  Test(46, 10,  6)
  Test(47, 10,  7)
  Test(48, 10,  8)
  Test(49, 10, 12)
  Test(50, 11,  3)
  Test(51, 11,  4)
  Test(52, 11,  5)
  Test(53, 11,  6)
  Test(54, 11,  7)
  Test(55, 11,  8)
  Test(56, 11, 12)
  Test(57, 12,  3)
  Test(58, 12,  4)
  Test(59, 12,  5)
  Test(60, 12,  6)
  Test(61, 12,  7)
  Test(62, 12, 10)
  Test(63, 12, 11)
})

test_that("TBR can swap over root", {
  expect_equal(TBR(tree5a, 1, c(7, 1)), read.tree(text = '(a, (d, (e, (c, b))));'))
  expect_equal(TBR(tree5a, 2, c(5, 1)), read.tree(text = '(a, (c, (b, (d, e))));'))
  expect_equal(TBR(tree5b, 1, c(7, 1)), read.tree(text = '((a, b), (d, (c, e)));'))
  expect_equal(TBR(tree5b, 4, c(7, 1)), read.tree(text = '((a, b), (d, (c, e)));'))
})

test_that("TBR works", {
  tree <- tree8
  ### expect_equal(TBR(tree, 3, 1 ), read.tree(text = "((a, ((b, (c, d)), (e, f))), (g, h));"))
  ### expect_warning(expect_identical(TBR(tree, 3, 2), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 3), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 4), tree))
  ### expect_warning(expect_identical(TBR(tree, 3, 44), tree))
  ### expect_equal(TBR(tree, 3, 5 ), read.tree(text = "((((a, b), (c, d)), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 6 ), read.tree(text = "(((b, (a, (c, d))), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 7 ), read.tree(text = "(((b, ((a, c), d)), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 8 ), read.tree(text = "(((b, (c, (a, d))), (e, f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 9 ), read.tree(text = "(((b, (c, d)), (a, (e, f))), (g, h));"))
  ### expect_equal(TBR(tree, 3, 10), read.tree(text = "(((b, (c, d)), ((a, e), f)), (g, h));"))
  ### expect_equal(TBR(tree, 3, 11), read.tree(text = "(((b, (c, d)), (e, (a, f))), (g, h));"))
  ### expect_equal(TBR(tree, 3, 12), read.tree(text = "(((b, (c, d)), (e, f)), (a, (g, h)));"))
  ### expect_equal(TBR(tree, 3, 13), read.tree(text = "(((b, (c, d)), (e, f)), ((g, a), h));"))
  ### expect_equal(TBR(tree, 3, 14), read.tree(text = "(((b, (c, d)), (e, f)), (g, (a, h)));"))
  
  tree <- tree8
  expect_equal(TBR(tree, 6, c(1 , 6)), read.tree(text = "((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 7)), read.tree(text = "((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(1 , 8)), read.tree(text = "((((a, b), (e, f)), (c, d)), (g, h));"))
  expect_equal(TBR(tree, 6, c(2 , 6)), TBR(tree, 6, c(2 , 7)))
  expect_equal(TBR(tree, 6, c(2 , 6)), TBR(tree, 6, c(2 , 8)))
  expect_equal(TBR(tree, 6, c(2 , 6)), read.tree(text = "((((a, b), (c, d)), (e, f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(3 , 6)), read.tree(text = "(((((c, d), a), b), (e, f)), (g, h));"))
  expect_warning(expect_identical(TBR(tree, 6, c(4 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 8, c(6 , 8)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(5 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 6)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 7)), tree))
  expect_warning(expect_identical(TBR(tree, 6, c(6 , 8)), tree))
  expect_equal(TBR(tree, 6, c(9 , 6)), read.tree(text = "(((a, b), ((c, d), (e, f))), (g, h));"))
  expect_equal(TBR(tree, 6, c(10, 6)), read.tree(text = "(((a, b), (((c, d), e), f)), (g, h));"))
  expect_equal(TBR(tree, 6, c(11, 6)), read.tree(text = "(((a, b), (((c, d), f), e)), (g, h));"))
  expect_equal(TBR(tree, 6, c(12, 6)), read.tree(text = "(((a, b), (e, f)), ((c, d), (g, h)));"))
  expect_equal(TBR(tree, 6, c(13, 6)), read.tree(text = "(((a, b), (e, f)), (((c, d), g), h));"))
  expect_equal(TBR(tree, 6, c(14, 6)), read.tree(text = "(((a, b), (e, f)), (((c, d), h), g));"))
  expect_warning(expect_identical(TBR(tree, 6, c(6, 15)), tree))
  
  expect_equal(TBR(tree, 4, c(1, 5)),  read.tree(text = "(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 6)),  read.tree(text = "(((a, (e, f)), (b, (c, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 7)),  read.tree(text = "(((a, (e, f)), (c, (b, d))), (g, h));"))
  expect_equal(TBR(tree, 4, c(1, 8)),  read.tree(text = "(((a, (e, f)), (d, (b, c))), (g, h));"))
  
  tree <- tree11 
  tree[["edge.length"]] = rep.int(1, 20) 
  expect_equal(TBR(tree11, 11, c(8, 17)), read.tree(text = '((j, k), (e, ((a, b), (c, (d, (i, (h, (g, f))))))));'))
  expect_equal(TBR(tree11, 11, c(2, 11)), read.tree(text = '((j, k), (e, (((a, b), (c, d)), (f, (g, (i, h))))));'))
  expect_warning(TBR(tree11, 10, c(2, 11)))
  expect_equal(TBR(tree11, 10, c(3, 11)), read.tree(text = '(e, ((c, d), ((a, b), ((j, k), (f, (g, (h, i)))))));'))
    
})

test_that("RootedTBR fails", {
  #  tree8 <- read.tree(text = "(((a, (b, (c, d))), (e, f)), (g, h));")
  #  tree11 <- read.tree(text = "((((a, b), (c, d)), e), ((f, (g, (h, i))), (j, k)));")

  expect_equal(TBR(tree8, 4, c(3, 7)), RootedTBR(tree8, 4, c(3, 7)))
  expect_equal(TBR(tree8, 4, c(1, 5)), RootedTBR(tree8, 4, c(1, 5)))
  expect_warning(RootedTBR(tree5a, edgeToBreak = 1))
  expect_warning(RootedTBR(tree5a, edgeToBreak = 2))
  expect_equal(RootedTBR(tree5a, edgeToBreak = 3, mergeEdges=6), read.tree(text = '(a, (c, (b, (d, e))));'))
  expect_silent(replicate(100, RootedTBR(tree5a)))
  expect_warning(RootedTBR(tree8, 4, c(13, 6)))
  expect_warning(RootedTBR(read.tree(text = '((a, b), (c, d));')))
})

test_that("RootedSPR fails", {
  expect_warning(RootedSPR(read.tree(text = '((a, b), (c, d));')))
  expect_warning(RootedSPR(tree8, edgeToBreak=1))
  expect_warning(RootedSPR(tree8, edgeToBreak=13))
  expect_warning(RootedSPR(tree8, edgeToBreak=14))
  warnTree1 <- read.tree(text = '((a, (b, (c, d))), (e, (f, (g, h))));')
  warnTree2 <- read.tree(text = '((a, (b, (c, d))), (((e, f), g), h));')
  attr(warnTree1, 'order') <- attr(warnTree2, 'order') <- 'preorder'
  expect_warning(RootedSPR(warnTree1, 3))
  expect_warning(RootedSPR(warnTree1, 10))
  expect_warning(RootedSPR(warnTree2, 9))
  expect_warning(RootedSPR(warnTree2, 8))
})

test_that("SPR is special case of TBR", {
  expect_equal(SPR(tree11, 3, 9), TBR(tree11, 3, c(3, 9)))
  expect_equal(SPR(tree11, 12, 9), TBR(tree11, 12, c(12, 9)))
  expect_equal(root(SPR(tree11, 1, 14), letters[1:5], resolve.root=TRUE), TBR(tree11, 1, c(1, 14)))
  expect_error(SPR(tree11, 1, 6))
})

#' @template MRS
CheckTreeSanity <- function (tree) {
  nTip <- length(tree[["tip.label"]])
  nNode <- tree[["Nnode"]]
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  aok <- TRUE
  expect_true(all(parent > nTip),
              info=paste0("Parent nodes on edge(s) ", paste(which(parent <= nTip), collapse=', '), 
                          " are tips (nTip = ", nTip, ')')
  )
  expect_equal(min(parent), nTip + 1,
               info=paste0("Root is numbered ", min(parent), "; expecting ", nTip + 1)
  )
  expect_false(min(parent) %in% child, 
               info=paste0("Root node (", min(parent), ") is child of edge ", paste0(which(min(parent) == child), collapse = ", "))
  )
  expect_true(all(seq_len(nTip) %in% child)) # No missing tips
  expect_equal(max(parent), nTip + nNode)
  tips <- child <= nTip
  expect_equal(sum(tips), nTip)
  expect_true(all(child[!tips] > parent[!tips]), info="Parent nodes must be > child nodes")
}

suppressWarnings(RNGversion("3.5.0"))
set.seed(0)
small_tree <- rtree(8)
large_tree <- rtree(80)  
test_that("NNI trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- NNI(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- NNI(large_tree))
})
test_that("SPR trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- SPR(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- SPR(large_tree))
})
test_that("TBR trees conform to phylo expectations", {
  for (i in 1:60)  CheckTreeSanity(small_tree <- TBR(small_tree))
  for (i in 1:250) CheckTreeSanity(large_tree <- TBR(large_tree))
})
