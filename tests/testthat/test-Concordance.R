library("TreeTools", quietly = TRUE)

test_that("_Concordance() handles null input", {
  expect_warning(expect_null(QuartetConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(PhylogeneticConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(ClusteringConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(MutualClusteringConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(SharedPhylogeneticConcordance(BalancedTree(8), NULL)))
})

test_that("_Concordance() handles tip mismatch", {
  char <- MatrixToPhyDat(cbind(c(a = 0, b = 0, c = 0, d = 1, e = 1)))
  tree <- BalancedTree(5)
  expect_warning(expect_null(QuartetConcordance(tree, char)),
                 "No overlap between tree labels and dataset.")
})

test_that(".TreeQuarters() gives complete list", {
  tree <- ape::read.tree(text = "(a, (b, (c, (d, ((e, f), (g, h))))));")
  coll <- CollapseNode(Preorder(tree), 12:13)
  expect_error(.TreeQuarters(tree), "must be in preorder")
  tree <- Preorder(tree)
  
  expect_equal(
    .TreeQuarters(tree),
    cbind("11" = c(a = 0, b = 1, c = 2, d = 3, e = 3, f = 3, g = 3, h = 3),
          "12" = c(0, 0, 1, 2, 3, 3, 3, 3),
          "13" = c(0, 0, 0, 1, 2, 2, 3, 3),
          "14" = c(rep(0, 4), 2, 3, 1, 1),
          "15" = c(rep(0, 4), 1, 1, 2, 3))
  )
  expect_equal(unname(.TreeQuarters(coll)),
               unname(.TreeQuarters(Preorder(tree))[, c(1, 4, 5)]))
})

test_that("QuartetConcordance(method = minhq)", {
  tree <- ape::read.tree(text = "(a, (b, (c, (d, ((e, f), (g, h))))));")
  tree <- Preorder(tree)
  mataset <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,  0,
                      0, 0, 1, 1, 1, 1, 1, 1,  0,
                      0, 0, 1, 1, 1, 1, 2, 2,  0,
                      1, 0, 0, 0, 1, 1, 1, 1,  0,
                      1, 0, 0, 0, 0, 1, 1, 1,  0,
                      0, 0, 0, 0, 1, 1, 2, 2,  0,
                      0, 0, 1, 1, 2, 2, 3, 3,  0,
                      0, 1, 2, 3, 0, 1, 2, 3,  0,
                      0, 1, 2, 0, 0, 2, 2, 3,  0), 9,
                    dimnames = list(letters[1:9], NULL))
  dat <- MatrixToPhyDat(mataset)
  
  # plot(tree); nodelabels();
  expect_concordance <- function(i, expectation, iq = "minhq") {
    expect_equal(
      QuartetConcordance(tree, dat[, i],
                         method = ifelse(iq == "iq", "iqtree", "minhq")),
      setNames(expectation, 11:15))
  }
  # Expectations computed by working through tables manually
  expect_concordance(1, c(NaN, NaN, 1, NaN, NaN))
  expect_concordance(2, c(1, NaN, NaN, NaN, NaN))
  expect_concordance(3, c(1, NaN, NaN, NaN, 1))
  expect_concordance(4, c(0, 0, 1, NaN, NaN))
  expect_concordance(5, c(0, 0, 4 / 6, 0, 1))
  expect_concordance(6, rep(NaN, 5))
  expect_concordance(7, c(1, NaN, NaN, NaN, NaN))
  expect_concordance(8, c(NaN, NaN, 0, NaN, NaN))
  expect_concordance(9, c(NaN, 0, 1 / 2, 0, NaN))
  
  # Values calculated from summing results above
  expect_equal(unname(QuartetConcordance(tree, dat, method = "minhq")), 
               c(5 + 3 + 1,
                 0,
                 12 + 8 + 4 + 1,
                 0,
                 4 + 3) /
                c(5 + 3 + 4 + 3 + 1,
                  4 + 3 + 2,
                  12 + 8 + 6 + 2 + 2,
                  6 + 2,
                  4 + 3))
  
  # Expectations computed by iq-tree
  expect_concordance(iq = "iq", 1, c(0, 0, 1, 0, 0))
  expect_concordance(iq = "iq", 2, c(1, 0, 0, 0, 0))
  expect_concordance(iq = "iq", 3, c(0.6, 0, 0, 0, 0.5))
  expect_concordance(iq = "iq", 4, c(0, 0, 2 / 3, 0, 0))
  expect_concordance(iq = "iq", 5, c(0, 0, 1 / 3, 0, 3 / 8))
  expect_concordance(iq = "iq", 6, c(0, 0, 0, 0, 0))
  expect_concordance(iq = "iq", 7, c(1 / 5, 0, 0, 0, 0))
  expect_concordance(iq = "iq", 8, rep(0, 5))
  expect_concordance(iq = "iq", 9, c(0, 0, 1 / 12, 0, 0))
  expect_equal(unname(QuartetConcordance(tree, dat, method = "iqtree")), 
               c(56.7, 0, 85.4, 0, 62.5) / 100, tolerance = 0.01)
  
  collapsed <- CollapseNode(tree, c(12, 13))
  expect_equal(
    unname(QuartetConcordance(collapsed, dat, method = "iqtree")),
    unname(QuartetConcordance(tree, dat, method = "iqtree"))[-c(2, 3)]
  )
})

test_that("QuartetConcordance(method = minh)", {
  tree <- Preorder(
    ape::read.tree(text = "(a, (b, (c, (d, ((e, f), (g, h))))));"))
  mataset <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,  0,
                      0, 0, 1, 1, 1, 1, 1, 1,  0,
                      0, 0, 1, 1, 1, 1, 2, 2,  0,
                      1, 0, 0, 0, 1, 1, 1, 1,  0,
                      1, 0, 0, 0, 0, 1, 1, 1,  0,
                      0, 0, 0, 0, 1, 1, 2, 2,  0,
                      0, 0, 1, 1, 2, 2, 3, 3,  0,
                      0, 1, 2, 3, 0, 1, 2, 3,  0,
                      0, 1, 2, 0, 0, 2, 2, 3,  0), 9,
                    dimnames = list(letters[1:9], NULL))
  dat <- MatrixToPhyDat(mataset)
  
  # plot(tree); nodelabels();
  expect_concordance <- function(i, expectation) {
    expect_equal(
      QuartetConcordance(tree, dat[, i], method = "minh", n = 1250),
      setNames(expectation, 11:15), tolerance = 0.05)
  }

  # Expectations computed by iq-tree
  expect_concordance(1, c(0, 0, 1, 0, 0))
  expect_concordance(2, c(1, 0, 0, 0, 0))
  expect_concordance(3, c(0.6, 0, 0, 0, 0.5))
  expect_concordance(4, c(0, 0, 2 / 3, 0, 0))
  expect_concordance(5, c(0, 0, 1 / 3, 0, 3 / 8))
  expect_concordance(6, c(0, 0, 0, 0, 0))
  expect_concordance(7, c(1 / 5, 0, 0, 0, 0))
  expect_concordance(8, rep(0, 5))
  expect_concordance(9, c(0, 0, 1 / 12, 0, 0))
  expect_equal(tolerance = 0.05,
    unname(QuartetConcordance(tree, dat, method = "minh", n = 1234)),
    c(56.7, 0, 85.4, 0, 62.5) / 100)
  
  collapsed <- CollapseNode(tree, c(12, 13))
  expect_equal(tolerance = 0.05,
    QuartetConcordance(collapsed, dat, method = "minh", n = 1234),
    QuartetConcordance(tree, dat, method = "minh", n = 1234)[-(2:3)]
  )
})

test_that("QuartetConcordance() calculates correct values - weighting", {
  labels <- letters[1:5]
  tr <- PectinateTree(labels)
  # plot(tr); nodelabels(adj = c(4, 0.5))
  
  char <- MatrixToPhyDat(cbind(
    c(a = 0, b = 0, c = 0, d = 1, e = 1),
    c(0, 0, 1, 0, 1),
    c(0, 1, 0, 0, 1)
  )[, c(1:3, 1)])
  
  # Is quartet concordant with character?
  # Quartets labelled by omitted leaf
  C <- 1    # Concordant
  D <- 0    # Discordant
  I <- NaN  # Irrelevant
  expected_concordance <- cbind(
    ch1 = c(a = C, b = C, c = C, d = I, e = I),
    ch2 = c(D, D, I, C, I),
    ch3 = c(D, I, D, D, I),
    ch4 = c(a = C, b = C, c = C, d = I, e = I))
  
  conc <- vapply(labels, function(lab) {
    quartet <- DropTip(tr, lab)
    vapply(1:4, function(i) QuartetConcordance(quartet, char[, i]), double(1))
  }, c(ch1 = NaN, ch2 = NaN, ch3 = NaN, ch4 = NaN))
  expect_equal(conc, t(expected_concordance))
  
  split_relevant <- list(
    "8" = letters[3:5],
    "9" = letters[1:3])
  
  expect_equal(
    QuartetConcordance(tr, char, method = "split"), # 8 = 0.75; 9 = 0.5
    vapply(split_relevant, function(splitQ) {
      sum(expected_concordance[splitQ, ], na.rm = TRUE) /
        sum(!is.nan(expected_concordance[splitQ, ]))
    }, double(1))
  )
  
  # Same proportions, different number of chars -> weights
  char10 <- char[, rep(1:4, each = 10)]
  expect_equal(
    QuartetConcordance(tr, char, method = "split"),
    QuartetConcordance(tr, char10, method = "split"))
})

test_that("QuartetConcordance() calculates correct values", {
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,  0,
                      0, 1, 0, 1, 0, 1, 0, 1,  0,
                      0, 0, 0, 1, 0, 1, 1, 1,  0,
                      0, 0, 0, 0, 1, 1, 2, 2,  0,
                      0, 0, 1, 1, 2, 2, 3, 3,  0,
                      0, 1, 2, 3, 0, 1, 2, 3,  0), 9,
                    dimnames = list(paste0("t", 1:9), NULL))
  expect_error(QuartetConcordance(tree, mataset),
               "`dataset` must be a phyDat object")
  dat <- MatrixToPhyDat(mataset)
  expect_equal(unname(QuartetConcordance(tree, dat[, 1])), rep(1, 5))

  # plot(tree); nodelabels();
  expect_equal(QuartetConcordance(tree, dat[, 2]),
               c("10" = 1/9, "11" = 0, "12" = 0,
                 "13" = 1/9, "14" = 0, "15" = 0)[names(as.Splits(tree))])
  
  allQuartets <- combn(8, 4)
  for (charI in seq_len(ncol(mataset))) {
    qc <- QuartetConcordance(tree, dat[, charI])
    for (splitI in seq_along(splits)) {
      split <- splits[[splitI]]
      logiSplit <- as.logical(split)
      case <- apply(allQuartets, 2, function (q) {
        qSplit <- logiSplit[q]
        qChar <- mataset[q, charI]
        if (identical(unique(table(qSplit)), 2L) &&
            identical(unique(table(qChar)), 2L)) {
          tbl <- table(qSplit, qChar)
          tab <- paste0(sort(tbl[tbl > 0]), collapse = "")
          switch(tab,
                 "1111" = FALSE,
                 "112" = NA,
                 "13" = NA,
                 "22" = TRUE,
                 "4" = NA,
                 stop(q, ": ", tab)
          )
        } else {
          NA
        }
      })
      expect_equal(sum(case, na.rm = TRUE) / sum(!is.na(case)),
                   unname(qc[as.character(names(split))]))
    }
  }
  
  expect_equal(QuartetConcordance(tree, dat[, c(1:4, 6)]),
               c("10" = (36 + 2 + 9 + 12) / (36 + 18 + 18 + 12 + 6),
                 "11" = ( 6 + 0 + 6 +  2) / ( 6 +  9 +  6 +  2 + 1),
                 "12" = ( 6 + 0 + 0 +  2) / ( 6 +  9 +  9 +  2 + 1),
                 "13" = (36 + 2 + 9 + 12) / (36 + 18 + 18 + 12 + 6),
                 "14" = ( 6 + 0 + 0 +  7) / ( 6 +  9 +  9 +  7 + 1),
                 "15" = ( 6 + 0 + 6 +  7) / ( 6 +  9 +  6 +  7 + 1))[
                   names(as.Splits(tree))]
  )
})

test_that("QuartetConcordance() handles ambiguity", {
  tree <- BalancedTree(12)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, "{01}", 0, 0, "{01}", 1, 1, "-", 1, 1, "-",
                      0, 1, "?", 0, 1, "?", 0, 1, "(01)", 0, 1, "(01)",
                      0, 0, "?", 0, 1, "(12)", 0, 1, "(12)", 1, 1, "(12)",
                      0, 0, "?", 0, 0, "?", 1, 1, "?", 2, 2, "?",
                      0, 0, "?", 0, 0, "?", 0, 0, "-", 0, 0, "-",
                      rep("?", 12),
                      0, 1, "?", 2, 3, "?", 0, 1, "-", 2, 3, "-"), 12,
                    dimnames = list(paste0("t", 1:12), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  expectation <- unname(QuartetConcordance(tree, dat)[
    c("14", "16", "18", "19", "21", "23")])
  expect_equal(
    unname(QuartetConcordance(DropTip(tree, paste0("t", 3 * 1:4)), dat)),
    expectation[!is.na(expectation)]
  )
  
  expectation <- unname(QuartetConcordance(tree, dat)[
    c("14", "15", "17", "19", "20", "22")])
  expect_equal(
    unname(QuartetConcordance(DropTip(tree, paste0("t", 3 * 1:4)), dat)),
    expectation[!is.na(expectation)]
  )
})

test_that("QuartetConcordance() handles incomplete data", {
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      rep("?", 8)), 8,
                    dimnames = list(paste0("t", 1:8), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  expect_equal(unname(QuartetConcordance(tree, dat)), rep(NA_real_, 5))
})

dataset <- congreveLamsdellMatrices[[10]][, 1]
tree <- TreeTools::NJTree(dataset)

ConcordantInformation(tree, dataset)["noise"]
TreeLength(tree, dataset, concavity = "prof")

test_that("ConcordantInformation() works", {
  data(congreveLamsdellMatrices)
  dat <- congreveLamsdellMatrices[[10]]
  tree <- TreeTools::NJTree(dat)
  
  ci <- ConcordantInformation(tree, dat)
  expect_equal(expect_warning(Evaluate(tree, dat)), ci)
  expect_equal(TreeLength(tree, dat, concavity = "prof"),
               unname(ci["noise"]))
  expect_equal(Log2Unrooted(22), unname(ci["treeInformation"]))
  expect_equal(sum(apply(PhyDatToMatrix(dat), 2, CharacterInformation)),
               unname(ci["informationContent"]))
  
  dataset <- MatrixToPhyDat(cbind(setNames(c(rep(1, 11), 2:5), paste0("t", 1:15))))
  tree <- TreeTools::PectinateTree(length(dataset))
  expect_error(ConcordantInformation(tree, dataset))
  # expect_equal(0, unname(ci["signal"]))
  # expect_equal(0, unname(ci["noise"]))
  
  dataset <- MatrixToPhyDat(c(a = 1, b = 2, c = 1, d = 2, e = 3, f = 3))
  tree <- TreeTools::PectinateTree(dataset)
  ci <- expect_warning(ConcordantInformation(tree, dataset))
  expect_equal(c(signal = log2(3)), ci["signal"])
  expect_equal(c(noise = log2(3)), ci["noise"])
  expect_equal(c(ignored = CharacterInformation(c(0,0,1,1,2,2)) - 
                   log2(3) - log2(3)), ci["ignored"])
  
})
