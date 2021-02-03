context('Jackknife.R')

test_that("Jackknife supports are correct", {
  true_tree <-  ape::read.tree(text = "((((((A,B),C),D),E),F),out);")
  start_tree <- ape::read.tree(text = "(((((A,D),B),E),(C,F)),out);")
  dataset <- TreeTools::StringToPhyDat('1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000',
                            1:7, byTaxon = FALSE)
  names(dataset) <- c(LETTERS[1:6], 'out')
  
  suppressWarnings(RNGversion("3.5.0")) # Until we can require R3.6.0
  set.seed(0)
  
  strict <- TreeSearch(start_tree, dataset, verbosity = 0)
  expect_equal(1, length(unique(list(true_tree), list(start_tree)))) # Right tree found
  jackTrees <- Jackknife(strict, dataset, resampleFreq = 4/7, searchIter = 24L,
                         searchHits = 7L, EdgeSwapper=RootedTBRSwap, jackIter = 20L,
                         verbosity = 0L)
  # Note: one cause of failure could be a change in characters sampled, due to randomness
  expect_true(length(unique(jackTrees)) > 2L)
})

test_that("Jackknife ouputs good for node.labels", {
  library('TreeTools') # for as.phylo
  
  # jackTrees will usually be generated with Jackknife(), but for simplicity:
  jackTrees <- as.phylo(1:100, 8)
  
  tree <- as.phylo(0, 8)
  expect_equal(c('', '', '0.13', '0.08', '0.14', '1', '1'),
               JackLabels(tree, jackTrees, plot = FALSE))
  
  tree <- RootTree(as.phylo(0, 8), c('t1', 't4'))
  expect_equal(c('', '0.08', '0.13', '', '0.14', '1', '1'),
               JackLabels(tree, jackTrees, plot = FALSE))
  })
