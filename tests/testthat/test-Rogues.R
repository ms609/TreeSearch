test_that("Rogues found", {
  
  library("TreeTools", quietly = TRUE)
  trees <- AddTipEverywhere(BalancedTree(8), 'Rogue')
  instab <- TipInstability(trees)
  expect_equal('Rogue', names(which.max(instab)))
  
  ci <- TipVolatility(trees)
  expect_equal('Rogue', names(which.max(ci)))
  
  dists <- TreeDist::PhylogeneticInfoDistance(trees, normalize = TRUE)
  expect_equal(mean(dists) - 0, max(ci))
               
  expect_equal(8L, NTip(BestConsensus(trees)))
  expect_equal(8L, NTip(Roguehalla(trees, 1)))
  
  
  trees[] <- lapply(trees, AddTip, 'Rogue', 'Rogue2')
  ci <- TipVolatility(trees)
  expect_equal(c('Rogue', 'Rogue2'), names(ci[ci == max(ci)]))
  
  # Interesting aside: Majority rule consensus favours balanced splits!
  bc <- BestConsensus(trees)
  expect_equal(10L, NTip(bc))
  expect_equal(9L, bc$Nnode)
  
  bc <- BestConsensus(trees[-11])
  expect_equal(8L, NTip(bc))
  expect_equal(7L, bc$Nnode)
  expect_equal(10L, NTip(Roguehalla(trees[-11], 1)))
  expect_equal(8L, NTip(Roguehalla(trees[-11], 2)))
})

test_that("Wilkinson & Crotti's examples are satisfied", {
  scaffold <- BalancedTree(c(6:4, 1:3))
  fig2 <- list(AddTip(scaffold, '3', 'X'),
               AddTip(scaffold, '4', 'X'))
  trees <- fig2
  expect_equal(match('X', TipLabels(fig2)), 
               unname(which.max(TipVolatility(fig2))))
  
  fig2b <- fig2[rep(1:2, c(67, 33))]
  expect_equal('X', names(which.max(TipVolatility(fig2b))))
  
  fig3 <- lapply(list(AddTip(scaffold, '1', 'X'),
                      AddTip(scaffold, '6', 'X')), AddTip, 'X', 'Y')
  
  trees <- fig3
  tr3 <- TipVolatility(fig3)
  expect_equal(c('X', 'Y'), names(tr3[tr3 == max(tr3)]))
  
  fig3b <- fig3[rep(1:2, c(60, 40))]
  tr3b <- TipVolatility(fig3b)
  expect_equal(c('X', 'Y'), names(tr3b[tr3b == max(tr3b)]))
  
  fig3c <- lapply(fig3b, drop.tip, names(tr3b[tr3b == max(tr3b)])) 
  expect_true(all(TipVolatility(fig3c) == 0))
  
  Tree <- function (txt) ape::read.tree(text = txt)
  fig4 <- list(Tree('((1, 2)68, 3, W, X, Y, Z, 4, (5, 6)70);'),
               Tree('((1, 2)74, 3, X, Y, Z, 4, (5, 6)74);'),
               Tree('(((1, 2)80, 3)55, Y, Z, (4, (5, 6)74)54);'),
               Tree('(((1, 2)100, 3)62, Z, (4, (5, 6)82)61);'),
               Tree('(((1, 2)100, 3)100, (4, (5, 6)100)100);'))
  siScores <- vapply(fig4, SplitwiseInfo, 0)
  expect_equal(3, which.max(siScores))
})
