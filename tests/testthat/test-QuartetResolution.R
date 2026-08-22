test_that("QuartetResolution()", {
  expect_equal(
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")),
    c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
      2, 2, 2, 1, 3, 2, 3, 3, 2, 1, 3, 3, 3, 2, 3, 2, 2, 2, 3, 3, 2, 
      2, 1, 2, 3, 2, 1, 1, 2, 1, 3, 2, 3, 3, 2, 1, 1, 1, 3, 1, 2, 1, 
      2, 1, 3, 3, 2, 1, 2, 1, 2, 2, 3))
  expect_equal(
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")),
    QuartetResolution(inapplicable.trees[["Vinther2008"]],
                      c("Nemertean", "Halkieria", "Wiwaxia", "Acaenoplax"))
  )
})

test_that("QuartetResolution() handles an unresolved (star) quartet", {
  library("TreeTools", quietly = TRUE)
  # A collapsed polytomy across all four focal tips, the shape
  # MaximizeParsimony(collapse = TRUE) produces when no character resolves
  # their relationship.
  tips <- c("Lingula", "Halkieria", "Wiwaxia", "Acaenoplax")
  tree <- as.phylo(0, 4)
  tree$tip.label <- tips
  internalNode <- tree$edge[tree$edge[, 2] > 4, 2]
  starTree <- CollapseNode(tree, internalNode)
  trees <- structure(list(starTree), class = "multiPhylo")

  expect_equal(QuartetResolution(trees, tips), NA_real_)
})
