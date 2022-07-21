library("TreeTools", quietly = TRUE)

test_that("Profile fails gracefully", {
  dataset <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 3, f = 3))
  expect_warning(PrepareDataProfile(dataset))
  expect_warning(MaximizeParsimony(dataset, concavity = "pr"))
})

test_that("Constraints work", {
  constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
  characters <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  set.seed(0)
  ewResults <- MaximizeParsimony(characters,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)
  expect_equal(PectinateTree(letters[1:6]), ewResults[[1]])
  expect_equal(c(seed = 0, start = 1, final = 0),
               attr(ewResults, "firstHit"))
  expect_equal(names(ewResults), "start_1")
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, concavity = "p",
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)[[1]])
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, concavity = 10,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint)[[1]])
  # Start tree not consistent with constraint
  dataset <- characters
  tree <- PectinateTree(c("a", "c", "f", "d", "e", "b"))
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, 
                                 PectinateTree(c("a", "c", "f", "d", "e", "b")),
                                 ratchIter = 0, constraint = constraint)[[1]])
  
  
  dataset <- MatrixToPhyDat(matrix(c(0, 0, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 0, 0, 0), ncol = 2,
                                   dimnames = list(letters[1:7], NULL)))
  constraint <- MatrixToPhyDat(matrix(c(0, 0, 1, "?", 1, 1,
                                        1, 1, 1,   1, 0, 0), ncol = 2,
                                      dimnames = list(letters[1:6], NULL)))
  cons <- consensus(MaximizeParsimony(dataset, constraint = constraint))
  expect_true(as.Splits(as.logical(c(0, 0, 1, 1, 1)), letters[c(1:3, 5:6)]) %in% 
                as.Splits(DropTip(cons, c("d", "g"))))
  
  expect_true(as.Splits(as.logical(c(0, 0, 0, 0, 1, 1)), letters[1:6]) %in% 
                as.Splits(DropTip(cons, "g")))
  
})

test_that("Inconsistent constraints fail", {
  constraint <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_error(MaximizeParsimony(constraint,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint))
})

test_that("MaximizeParsimony() times out", {
  # Do not run on CRAN: Writing R Extensions discourages testing timings
  skip_if(Sys.getenv("GITHUB_PAT") == "") # Run only on GH Actions
  
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[42]]
  startTime <- Sys.time()
  MaximizeParsimony(dataset, ratchIter = 10000, tbrIter = 1, maxHits = 1,
                    maxTime = 0)
  expect_gt(as.difftime(5, units = "secs"), Sys.time() - startTime)
})

test_that("Seed trees retained", {
  tree1 <- read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  tree2 <- read.tree(text = "(a, (b, (c, (f, (e, d)))));")
  badTree <- read.tree(text = "(f, (b, (c, (a, (e, d)))));")
  dat <- StringToPhyDat("110000 110000 111000 111000 111100 111001",
                        letters[1:6], byTaxon = FALSE)
  results <- MaximizeParsimony(dataset = dat, 
                               tree = c(tree1, tree2, badTree),
                               ratchIter = 0, verbosity = 4)
  expect_equal(attr(results, "firstHit"),
               c(seed = 2, start = 0, final = 0))
})

test_that("Mismatched tree/dataset handled with warnings", {
  treeAf <- read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  treeBg <- read.tree(text = "(g, (b, (c, (d, (e, f)))));")
  datAf <- StringToPhyDat("110000 110000 111100 111000",
                              letters[1:6], byTaxon = FALSE)
  datAe <- StringToPhyDat("11000 11000 11110 11100",
                              letters[1:5], byTaxon = FALSE)
  datAg <- StringToPhyDat("1100000 1100000 1111000 1110000",
                              letters[1:7], byTaxon = FALSE)
  
  QP <- function (...) MaximizeParsimony(..., ratchIter = 0, maxHits = 1,
                                         verbosity = 0)
  
  expect_equal(5, unname(NTip(expect_warning(QP(datAf, treeBg)))))
  expect_equal(5, unname(NTip(expect_warning(QP(datAe, treeAf)))))
  expect_equal(6, unname(NTip(expect_warning(QP(datAg, treeAf)))))
  expect_equal(5, unname(NTip(expect_warning(QP(datAf, treeBg, constraint = datAe)))))
  expect_equal(6, unname(NTip(QP(datAf, treeAf, constraint = datAe))))
  expect_equal(6, unname(NTip(expect_warning(QP(datAf, treeAf, constraint = datAg)))))
})

test_that("Root retained if not 1", {
  tr <- RootTree(BalancedTree(8), "t5")
  dataset <- StringToPhyDat("11000000 11100000 11110000 11111000",
                            paste0("t", 1:8), byTaxon = FALSE)
  
  mpt <- MaximizeParsimony(dataset, tr)
  expect_equal(5, mpt[[1]]$edge[14, 2])
})

test_that("Resample() fails and works", {
  expect_error(Resample(0))
  dataset <- MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0, 0, 0),
    b = c(0, 0, 0, 0, 0, 0),
    c = c(1, 1, 0, 0, 0, 1),
    d = c(1, 1, 0, 0, 1, 0),
    e = c(1, 1, 1, 1, 1, 1),
    f = c(1, 1, 1, 1, 1, 1)))
  
  expect_error(Resample(dataset, method = "ERROR"))
  expect_error(Resample(dataset, proportion = 0))
  expect_error(Resample(dataset, proportion = 6 / 7))

  nRep <- 42L # Arbitrary number to balance runtime vs false +ves & -ves
  bal <- as.Splits(BalancedTree(dataset))
  
  skip_if_not_installed("TreeTools", "1.4.5.9003") # postorder / as.Splits order
  jackTrees <- replicate(nRep, Resample(dataset, NJTree(dataset), verbosity = 0L))
  jackSplits <- as.Splits(unlist(jackTrees, recursive = FALSE))
  jackSupport <- rowSums(vapply(jackSplits, function(sp) bal %in% sp,
                                logical(3)))
  
  skip_if_not_installed("TreeTools", "1.6.0.9002") # names
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(jackSupport, tolerance = 0.2,
               c("8" = 1/2, "9" = 1, "10" = 1/2, "11" = 0)[names(bal)] *
                 sum(vapply(jackTrees, length, 1L)))
  
  bootTrees <- replicate(nRep, Resample(dataset, method = "bootstrap",
                                        verbosity = 0))
  #bootSupport <- rowSums(vapply(lapply(bootTrees, `[[`, 1),
  bootSupport <- rowSums(vapply(unlist(bootTrees, recursive = FALSE),
                                function(tr) bal %in% as.Splits(tr),
                                logical(3)))
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(bootSupport, tolerance = 0.2,
               c("8" = 1/2, "9" = 1, "10" = 1/2, "11" = 0)[names(bal)] * 
                 sum(vapply(bootTrees, length, 1L)))
    
})

test_that(".CombineResults() handles duplicates", {
  x <- structure(
    array(c(
      rep(1L, 8),
      rep(2L, 8),
      rep(3L, 8),
      rep(2L, 8),
      rep(1L, 8)
      ),
      dim = c(4, 2, 5)),
    firstHit = c(start = 5, test = 0, end = 0)
  )
  y <- array(c(rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8)),
          dim = c(4, 2, 5)
          )
  expect_warning(.CombineResults(x, y, stage = "test"))
  uX <- structure(unique(x, MARGIN = 3L),
                  firstHit = c(start = 3, test = 0, end = 0))
  expect_equal(attr(.CombineResults(uX, y, stage = "test"), "firstHit"),
               c(start = 3, test = 1, end = 0))
               
})
