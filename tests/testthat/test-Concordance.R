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

test_that("QuartetConcordance() works", {
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
