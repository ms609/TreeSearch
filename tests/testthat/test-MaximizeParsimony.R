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
