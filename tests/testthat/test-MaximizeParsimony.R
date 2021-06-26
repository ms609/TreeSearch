library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

test_that("Profile fails gracefully", {
  dataset <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 3, f = 3))
  expect_warning(PrepareDataProfile(dataset))
  expect_warning(MaximizeParsimony(dataset, concavity = 'pr'))
})

test_that("Constraints work", {
  constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
  characters <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  set.seed(0)
  ewResults <- MaximizeParsimony(characters,
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint)
  expect_equal(PectinateTree(letters[1:6]), ewResults[[1]])
  expect_equal(c(seed = 0, start = 1, final = 0),
               attr(ewResults, 'firstHit'))
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, concavity = 'p',
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint)[[1]])
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, concavity = 10,
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint)[[1]])
  # Start tree not consistent with constraint
  dataset <- characters
  tree <- PectinateTree(c('a', 'c', 'f', 'd', 'e', 'b'))
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, 
                                 PectinateTree(c('a', 'c', 'f', 'd', 'e', 'b')),
                                 ratchIter = 0, constraint = constraint)[[1]])
  
})

test_that("Constrained NJ trees work", {
  dataset <- MatrixToPhyDat(matrix(
   c(0, 1, 1, 1, 0, 1,
     0, 1, 1, 0, 0, 1), ncol = 2,
   dimnames = list(letters[1:6], NULL)))
  constraint <- MatrixToPhyDat(c(a = 0, b = 0, c = 0, d = 0, e = 1, f = 1))
  expect_equal(ape::read.tree(text = "(a, (d, ((c, b), (e, f))));"),
               ConstrainedNJ(dataset, constraint))
  expect_equal(NJTree(dataset), ConstrainedNJ(dataset, dataset))
})

test_that("Inconsistent constraints fail", {
  constraint <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  #skip_if(T)
  expect_error(MaximizeParsimony(constraint,
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint))
})

test_that("MaximizeParsimony() times out", {
  data('congreveLamsdellMatrices', package = 'TreeSearch')
  dataset <- congreveLamsdellMatrices[[42]]
  startTime <- Sys.time()
  MaximizeParsimony(dataset, ratchIter = 10000, tbrIter = 1, maxHits = 1,
                    maxTime = 0)
  expect_gt(as.difftime(5, units = 'secs'), Sys.time() - startTime)
})

test_that("Mismatched tree/dataset handled with warnings", {
  treeAf <- read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  treeBg <- read.tree(text = "(g, (b, (c, (d, (e, f)))));")
  datAf <- StringToPhyDat('110000 110000 111100 111000',
                              letters[1:6], byTaxon = FALSE)
  datAe <- StringToPhyDat('11000 11000 11110 11100',
                              letters[1:5], byTaxon = FALSE)
  datAg <- StringToPhyDat('1100000 1100000 1111000 1110000',
                              letters[1:7], byTaxon = FALSE)
  
  QP <- function (...) MaximizeParsimony(..., ratchIter = 0, maxHits = 1,
                                         verbosity = 0)
  
  expect_equal(5, NTip(expect_warning(QP(datAf, treeBg))))
  expect_equal(5, NTip(expect_warning(QP(datAe, treeAf))))
  expect_equal(6, NTip(expect_warning(QP(datAg, treeAf))))
  expect_equal(5, NTip(expect_warning(QP(datAf, treeBg, constraint = datAe))))
  expect_equal(6, NTip(QP(datAf, treeAf, constraint = datAe)))
  expect_equal(6, NTip(expect_warning(QP(datAf, treeAf, constraint = datAg))))
})

test_that("Root retained if not 1", {
  tr <- RootTree(BalancedTree(8), 't5')
  dataset <- StringToPhyDat('11000000 11100000 11110000 11111000',
                            paste0('t', 1:8), byTaxon = FALSE)
  
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
  
  expect_error(Resample(dataset, method = 'ERROR'))
  expect_error(Resample(dataset, proportion = 0))
  expect_error(Resample(dataset, proportion = 6 / 7))

  nRep <- 42L # Arbitrary number to balance runtime vs false +ves & -ves
  bal <- as.Splits(BalancedTree(dataset))
  
  jackTrees <- replicate(nRep, Resample(dataset, NJTree(dataset), verbosity = 0L))
  jackSplits <- as.Splits(unlist(jackTrees, recursive = FALSE))
  jackSupport <- rowSums(vapply(jackSplits, function (sp) in.Splits(bal, sp),
                                logical(3)))
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(c(1/2, 1, 0) * sum(vapply(jackTrees, length, 1L)), jackSupport,
               tolerance = 0.2)
  
  bootTrees <- replicate(nRep, Resample(dataset, method = 'bootstrap',
                                        verbosity = 0))
  #bootSupport <- rowSums(vapply(lapply(bootTrees, `[[`, 1),
  bootSupport <- rowSums(vapply(unlist(bootTrees, recursive = FALSE),
                                function (tr) in.Splits(bal, as.Splits(tr)),
                                logical(3)))
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(c(1/2, 1, 0) * sum(vapply(bootTrees, length, 1L)), bootSupport,
               tolerance = 0.2)
    
})
