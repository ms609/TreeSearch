test_that("Failures are graceful", {
  library("TreeTools", quietly = TRUE)
  data("inapplicable.datasets")
  dat <- inapplicable.phyData[[1]]
  unrooted <- RandomTree(dat, root = FALSE)
  expect_error(TreeLength(unrooted, dat), "must be rooted")
  expect_error(TreeLength(CollapseNode(RandomTree(dat, root = TRUE), 77), mo),
               "must be binary")
  
  # Same number of leaves, but different labels
  bal <- BalancedTree(dat)
  bal$tip.label[1:2] <- c("no1", "no2")
  expect_error(TreeLength(bal, dat), "Missing in `dataset`: no1, no2")
  
  mo <- PhyDat2Morphy(dat)
  on.exit(mo <- UnloadMorphy(mo))
  
  sparse <- DropTip(RandomTree(dat, root = FALSE), 10)
  expect_error(MorphyTreeLength(sparse, mo),
               "Number of taxa .* not equal to number of tips")
  expect_error(MorphyTreeLength(sparse, NA),
               "a valid Morphy pointer")
  
  expect_error(MorphyLength(sparse$edge[, 1], sparse$edge[, 2], mo, nTaxa = 0),
               mpl_translate_error(0))
  expect_error(MorphyLength(sparse$edge[, 1], sparse$edge[, 2], dat),
               "must be a Morphy pointer")
  
  expect_null(TreeLength(NULL))
  
})

test_that("Deprecations throw warning", {
  data("inapplicable.datasets")
  dat <- inapplicable.phyData[[1]]
  tree <- TreeTools::RandomTree(dat, root = TRUE)
  expect_equal(TreeLength(tree, dat),
               expect_warning(Fitch(tree, dat)))
  expect_equal(CharacterLength(tree, dat, compress = TRUE),
               expect_warning(FitchSteps(tree, dat)))
  
})

test_that("Morphy generates correct lengths", {
  ## Tree
  tree <- ape::read.tree(text = "((((((1,2),3),4),5),6),(7,(8,(9,(10,(11,12))))));")
  relabel <- ape::read.tree(text = "((6,(5,(4,(3,(2,1))))),(7,(8,(9,(10,(11,12))))));")
  trees <- list(tree, relabel)
  ## Test cases designed by Thomas Guillerme
  characters <- c("23--1??--032", # 0,  expect score = 5 
                  "1---1111---1", # 1,  expect score = 2
                  "1100----1100", # 2,  expect score = 3
                  "11-------100", # 3,  expect score = 2
                  "----1111---1", # 4,  expect score = 1
                  "01----010101", # 5,  expect score = 5
                  "01---1010101", # 6,  expect score = 5
                  "1??--??--100", # 7,  expect score = 2
                  "21--3??--032", # 8,  expect score = 5
                  "11--1??--111", # 9,  expect score = 2
                  "11--1000001-", # 10, expect score = 2
                  "01------0101", # 11, expect score = 4
                  "110--?---100", # 12, expect score = 3
                  "11--1??--111", # 13, expect score = 2
                  "210--100--21", # 14, expect score = 5
                  "????----1???", # 15, expect score = 0
                  "23--1----032", # 16, expect score = 5
                  "1----1----1-", # 17, expect score = 2
                  "-1-1-1--1-1-", # 18, expect score = 4
                  "23--1??--032", # 19, expect score = 5
                  "--------0101", # 20, expect score = 2
                  "10101-----01", # 21, expect score = 4
                  "011--?--0011", # 22, expect score = 3
                  "110--??--100", # 23, expect score = 3
                  "11--1000001-", # 24, expect score = 2
                  "21--1----012", # 25, expect score = 5
                  "11----111111", # 26, expect score = 1
                  "10101-----01", # 27, expect score = 4
                  "210210------", # 28, expect score = 4
                  "----1111----", # 29, expect score = 0
                  "230--??1--32", # 30, expect score = 5
                  "023--??1--32", # 31, expect score = 5
                  "023-???1--32", # 32, expect score = 4
                  "23--1?1--023", # 33, expect score = 5
                  "----1010----", # 34, expect score = 2
                  "------11---1", # 35, expect score = 1
                  "10----11---1", # 36, expect score = 3
                  "320--??3--21", # 37, expect score = 5
                  "000011110000"  # 38, expect score = 2
                  ) 
  ## Results
  expected_results <- c(5, 2, 3, 2, 1, 5, 5, 2, 5, 2, 2, 4, 3, 2, 5, 0, 5, 2,
                        4, 5, 2, 4, 3, 3, 2, 5, 1, 4, 4, 0, 5, 5, 4, 5, 2, 1, 
                        3, 5, 2)
  expected_minLength <- c(3, 0, 1, 1, 0, 1, 1, 1, 3, 0, 1, 1, 1, 0, 2, 0, 3, 0,
                          0, 3, 1, 1, 1, 1, 1, 2, 0, 1, 2, 0, 3, 3, 3, 3, 1, 0, 
                          1, 3, 1)
  expected_homoplasies <- expected_results - expected_minLength

  ##plot(tree); nodelabels(12:22); tiplabels(0:11)
  ## Run the tests
  for(test in seq_along(characters)) {
    morphyObj <- SingleCharMorphy(characters[test])
    tree_length <- MorphyTreeLength(tree, morphyObj)
    morphyObj <- UnloadMorphy(morphyObj)
    #if (tree_length != expected_results[test]) message("Test case", test - 1, characters[test], "unequal: Morphy calcluates",
    #  tree_length, "instead of", expected_results[test],"\n")
    expect_equal(tree_length, expected_results[test])
  }
  
  ## Test combined matrix
  bigPhy <- TreeTools::StringToPhyDat(paste0(characters, collapse = "\n"),
                                      tree$tip.label, 
                                      byTaxon = FALSE)
  profPhy <- TreeTools::StringToPhyDat(paste0(characters[-c(15, 29, 34)],
                                              collapse = "\n"),
                                       tree$tip.label, 
                                       byTaxon = FALSE)
  expect_identical(characters,
                   TreeTools::PhyToString(bigPhy, byTaxon = FALSE,
                                          concatenate = FALSE))
  expect_identical(paste0(collapse = '', 
                          vapply(characters, substr, start = 0, stop = 1,
                                 character(1))),
                   substr(TreeTools::PhyToString(bigPhy, ps = ";",
                                                 useIndex = TRUE,
                                                 byTaxon = TRUE,
                                                 concatenate = TRUE),
                    start = 0, stop = length(characters)))
  
  morphyObj <- PhyDat2Morphy(bigPhy)
  moSummary <- summary(morphyObj)
  expect_equal(c(length(bigPhy), attr(bigPhy, "nr"), length(bigPhy) - 1),
               c(moSummary$nTax, moSummary$nChar, moSummary$nInternal))
  tree_length <- MorphyTreeLength(tree, morphyObj)
  morphyObj <- UnloadMorphy(morphyObj)
  
  expect_equal("0123", moSummary$allStates)
  expect_equal(tree_length, sum(expected_results))
  expect_equal(tree_length, TreeLength(tree, bigPhy))
  expect_equal(tree_length, TreeLength(relabel, bigPhy))
  expect_equal(rep(tree_length, 2), TreeLength(trees, bigPhy))
  
  expected_fit <- expected_homoplasies / (expected_homoplasies + 6)
  tree_score_iw <- TreeLength(tree, bigPhy, concavity = 6)
  expect_equal(sum(expected_fit), tree_score_iw)
  expect_equal(tree_score_iw, TreeLength(relabel, bigPhy, concavity = 6))
  expect_equal(vapply(trees, TreeLength, double(1), bigPhy, concavity = 6),
               TreeLength(trees, bigPhy, concavity = 6))
  
  expect_equal(vapply(trees, TreeLength, double(1), profPhy, concavity = "p"),
               TreeLength(trees, profPhy, concavity = "profile"))
  

  ## Run the bigger tree tests
  bigTree <- ape::read.tree(
    text = "((1,2),((3,(4,5)),(6,(7,(8,(9,(10,((11,(12,(13,(14,15)))),(16,(17,(18,(19,20))))))))))));")
  bigChars <- c("11111---111---11---1")
  ## Results
  expected_results <- c(3)

  ## Run the tests
  for(test in 1:length(bigChars)) {
    phy <- TreeTools::StringToPhyDat(bigChars[test], bigTree$tip.label)
    # Presently a good test to confirm that PhyDat2Morphy works with single-character phys
    morphyObj <- PhyDat2Morphy(phy)
    on.exit(morphyObj <- UnloadMorphy(morphyObj))
    tree_length <- MorphyTreeLength(bigTree, morphyObj)
    
    expect_equal(tree_length, expected_results[test])
  }
})

test_that("(random) lists of trees are scored", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  mat <- congreveLamsdellMatrices[[42]]
  
  # Expected values calculated from 100k samples
  expect_gt(t.test(TreeLength(100, mat), mu = 318.5877)$p.val, 0.001)
  expect_gt(t.test(TreeLength(100, mat, 10L), mu = 17.16911)$p.val, 0.001)
  expect_gt(t.test(TreeLength(100, mat, "profile"), mu = 830.0585)$p.val, 0.001)
})

test_that("TreeLength() handles unrooted / non-preorder trees", {
  library("TreeTools", quietly = TRUE)
  data("congreveLamsdellMatrices", package = "TreeSearch")
  mat <- congreveLamsdellMatrices[[42]]
  
  set.seed(0)
  rand <- RandomTree(mat, root = TRUE)
  scores <- TreeLength(c(rand, Postorder(rand)), mat)
  expect_equal(scores[[1]], scores[[2]])
  expect_equal(TreeLength(RootTree(Postorder(rand), 1), mat), scores[[1]])
  expect_equal(TreeLength(RootTree(Preorder(rand), 1), mat), scores[[1]])
  
  set.seed(0)
  unrooted <- RandomTree(mat, root = FALSE)
  
  expect_equal(expect_warning(TreeLength(c(unrooted), mat),
                              "rooted on tip 1"),
               TreeLength(c(RootTree(unrooted, 1)), mat))
  
  expect_equal(TreeLength(RootTree(Postorder(unrooted), 1), mat),
               TreeLength(RootTree(unrooted, 1), mat))
  scores <- expect_warning(
    TreeLength(c(unrooted, Postorder(unrooted)), mat),
    "rooted on tip 1")
  expect_equal(scores[[1]], scores[[2]])
  expect_equal(TreeLength(RootTree(unrooted, 1), mat), scores[[1]])
})

test_that("TreeLength() handles subsetted trees", {
  data("inapplicable.datasets")
  dat <- inapplicable.phyData[[1]] 
  t8 <- as.phylo(1:4, 8, tipLabels = names(dat)[1:8])
  expect_equal(4, length(TreeLength(t8, dat)))
})

test_that("Profile scoring is reported correctly", {
  data("congreveLamsdellMatrices")
  dataset <- congreveLamsdellMatrices[[42]]
  prepDataset <- PrepareDataProfile(dataset)
  tree <- NJTree(prepDataset)
  edge <- Preorder(tree)$edge
  at <- attributes(prepDataset)
  profiles <- attr(prepDataset, "info.amounts")
  charSeq <- seq_along(prepDataset[[1]]) - 1L
  
  characters <- PhyToString(prepDataset, ps = '', useIndex = FALSE,
                            byTaxon = FALSE, concatenate = FALSE)
  startWeights <- at$weight
  morphyObjects <- lapply(characters, SingleCharMorphy)
  on.exit(morphyObjects <- vapply(morphyObjects, UnloadMorphy, integer(1)),
          add = TRUE)
  
  expect_equal(TreeLength(tree, dataset, "profile"),
               TreeLength(tree, prepDataset, "profile"))
  expect_equal(TreeLength(tree, dataset, "profile"),
               morphy_profile(edge, morphyObjects, startWeights, charSeq, 
                              profiles, Inf))
})

test_that("CharacterLength() fails gracefully", {
  expect_error(CharacterLength(as.phylo(1, 8), 1))
  
  data("inapplicable.datasets")
  dataset <- inapplicable.phyData[[12]]
  # Unlabelled leaves
  expect_error(CharacterLength(structure(list(), class = "phylo"), dataset))
  
  # Missing leaves
  expect_error(CharacterLength(as.phylo(1, 4), dataset))
  tMinus1 <- as.phylo(1, 42, tipLabels = names(dataset)[-1])
  expect_equal(CharacterLength(tMinus1, dataset[-1]),
               CharacterLength(tMinus1, dataset))
  expect_error(CharacterLength(as.phylo(1, 43), dataset))
  tPlus1 <- as.phylo(1, 44, tipLabels = c("extra", names(dataset)))
  expect_equal(CharacterLength(DropTip(tPlus1, "extra"), dataset),
               CharacterLength(tPlus1, dataset))
  expect_error(CharacterLength(as.phylo(1:2, 43, tipLabels = names(dataset)),
                               dataset))
  # no error:
  CharacterLength(as.phylo(1, 43, tipLabels = names(dataset)), dataset)
  
  skip_if_not_installed("phangorn")
  library("phangorn") # for phyDat subsetting
  expect_equal(c(53, 59, 6),
               as.numeric(table(CharacterLength(NJTree(dataset[1:4, ]),
                                                dataset[1:4], compress = TRUE))))
  
})

test_that("Character compression works", {
  data("inapplicable.datasets")
  dataset <- inapplicable.phyData[[12]]
  tree <- TreeTools::NJTree(dataset)
  expect_equal(length(CharacterLength(tree, dataset)), 137)
  expect_equal(length(MinimumLength(dataset)), 137)
  expect_equal(dim(Consistency(dataset, tree, nRelabel = 0)), c(137, 4))
  expect_equal(length(CharacterLength(tree, dataset, compress = TRUE)), 118)
  expect_equal(length(MinimumLength(dataset, compress = TRUE)), 118)
  expect_equal(dim(Consistency(dataset, tree, nRelabel = 0, compress = TRUE)),
               c(118, 4))
})

test_that("X_MorphyLength", {
  dataset <- congreveLamsdellMatrices[[42]]
  morphyObj <- PhyDat2Morphy(dataset)
  on.exit(UnloadMorphy(morphyObj))
  nTaxa <- mpl_get_numtaxa(morphyObj)
  
  tree <- NJTree(dataset)
  edgeList <- Preorder(tree$edge)
  edgeList <- edgeList[PostorderOrder(edgeList), ]
  parent <- edgeList[, 1]
  child <- edgeList[, 2]

  maxNode <- nTaxa + mpl_get_num_internal_nodes(morphyObj)
  rootNode <- nTaxa + 1L
  allNodes <- rootNode:maxNode
  
  parentOf <- parent[match(seq_len(maxNode), child)]
  parentOf[rootNode] <- rootNode # Root node's parent is a dummy node
  leftChild <- child[length(parent) + 1L - match(allNodes, rev(parent))]
  rightChild <- child[match(allNodes, parent)]

  expected <- MorphyLength(parent, child, morphyObj)
  
  expect_equal(expected,
               C_MorphyLength(parentOf, leftChild, rightChild, morphyObj))
  expect_equal(expected,
               GetMorphyLength(parentOf - 1, leftChild - 1, rightChild - 1,
                              morphyObj))
})
