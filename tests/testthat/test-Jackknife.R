context("Jackknife.R")

test_that("Jackknife supports are correct", {
  true_tree <-  ape::read.tree(text = "((((((A,B),C),D),E),F),out);")
  start_tree <- ape::read.tree(text = "(((((A,D),B),E),(C,F)),out);")
  dataset <- TreeTools::StringToPhyDat("1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000",
                            1:7, byTaxon = FALSE)
  names(dataset) <- c(LETTERS[1:6], "out")
  
  expect_error(Jackknife(unroot(true_tree), dataset))
  expect_error(Jackknife(start_tree, dataset, resampleFreq = 0))
  expect_error(Jackknife(start_tree, dataset, resampleFreq = 9/10))
  
  suppressWarnings(RNGversion("3.5.0"))
  set.seed(0)
  
  strict <- TreeSearch(start_tree, dataset, verbosity = 0)
  expect_equal(1, length(unique(list(true_tree), list(start_tree)))) # Right tree found
  jackTrees <- Jackknife(strict, dataset, resampleFreq = 4/7, searchIter = 24L,
                         searchHits = 7L, EdgeSwapper=RootedTBRSwap, 
                         jackIter = 20L, verbosity = 0L)
  
  # Note: one cause of failure could be a change in characters sampled, due to randomness
  expect_true(length(unique(jackTrees)) > 2L)
})

test_that("Jackknife ouputs good for node.labels", {
  library("TreeTools", quietly = TRUE) # for as.phylo
  
  # jackTrees will usually be generated with Jackknife(), but for simplicity:
  jackTrees <- as.phylo(1:100, 8)
  
  tree <- as.phylo(0, 8)
  expect_equal(JackLabels(tree, jackTrees, plot = FALSE),
               c(NA_real_, NA_real_, 0.13, 0.08, 0.14, 1, 1))
  
  tree <- RootTree(as.phylo(0, 8), c("t1", "t4"))
  expect_equal(JackLabels(tree, jackTrees, plot = FALSE),
               c(NA_real_, 0.08, 0.13, NA_real_, 0.14, 1, 1))
  
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("plot-jackknife", function() {
    expect_equal(as.double(JackLabels(tree, jackTrees, plot = FALSE)[-c(1, 4)]),
                 unname(JackLabels(tree, jackTrees)))
  })
})

test_that("JackLabels() handles multiple trees per iteration", {
  tree <- BalancedTree(5)
  plot(tree)
  nodelabels()
  dispute8 <- ape::read.tree(text = "(((t1, t3), t2), (t4, t5));")
  disagree <- ape::read.tree(text = "(((t5, t2), t3), (t4, t1));")
  jackTrees <- list(
    c(dispute8, dispute8),
    c(tree, tree),
    c(dispute8, tree),
    c(disagree, disagree, disagree),
    BalancedTree(5)
  )
  expect_equal(
    JackLabels(tree, jackTrees),
    structure(c("7" = 4 / 5, "8" = 2 / 4), decisive = c("7" = 5, "8" = 4))
  )
})
