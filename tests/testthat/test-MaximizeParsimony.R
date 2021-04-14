library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)

test_that("constraints work", {
  constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
  characters <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters,
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint)[[1]])
  # Start tree not consistent with constraint
  expect_equal(PectinateTree(letters[1:6]),
               MaximizeParsimony(characters, 
                                 PectinateTree(c('a', 'c', 'f', 'd', 'e', 'b')),
                                 ratchIter = 0, constraint = constraint)[[1]])
  
})

test_that("inconsistent constraints fail", {
  constraint <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  expect_error(MaximizeParsimony(constraint,
                                 PectinateTree(c('a', 'b', 'f', 'd', 'e', 'c')),
                                 ratchIter = 0, constraint = constraint))
})

test_that("Root retained if not 1", {
  tr <- RootTree(BalancedTree(8), 't5')
  dataset <- StringToPhyDat('11000000 11100000 11110000 11111000',
                            paste0('t', 1:8), byTaxon = FALSE)
  
  mpt <- MaximizeParsimony(dataset, tr)
  expect_equal(5, mpt[[1]]$edge[14, 2])
})

test_that("Resample() fails and works", {
  # Not sure why this is necessary:
  library("TreeTools", quietly = TRUE, warn.conflicts = FALSE)
  
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
  expect_error(Resample(dataset, proportion = 6/7))

  nRep <- 42L # Arbitrary number to balance runtime vs false +ves & -ves
  bal <- as.Splits(BalancedTree(dataset))
  jackTrees <- replicate(nRep, Resample(dataset, NJTree(dataset), verbosity = 0L))
  jackSplits <- as.Splits(unlist(jackTrees, recursive = FALSE))
  jackSupport <- rowSums(vapply(jackSplits, function (sp) bal %in% sp,
                                logical(3)))
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(c(1/2, 1, 0) * sum(vapply(jackTrees, length, 1L)), jackSupport,
               tolerance = 0.2)
  
  bootTrees <- replicate(nRep, Resample(dataset, method = 'bootstrap',
                                        verbosity = 0))
  #bootSupport <- rowSums(vapply(lapply(bootTrees, `[[`, 1),
  bootSupport <- rowSums(vapply(unlist(bootTrees, recursive = FALSE),
                                function (tr) {bal %in% as.Splits(tr)},
                                logical(3)))
  # This test could be replaced with a more statistically robust alternative!
  expect_equal(c(1/2, 1, 0) * sum(vapply(bootTrees, length, 1L)), bootSupport,
               tolerance = 0.2)
    
})
  