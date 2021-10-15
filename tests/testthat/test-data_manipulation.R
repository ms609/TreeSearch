test_that("PrepareDataProfile() handles empty matrices", {
  dat <- TreeTools::MatrixToPhyDat(matrix(c(0, 1, rep('?', 5)),
                                          dimnames = list(letters[1:7], NULL)))
  expectation <- dat[0]
  attr(expectation, 'info.amounts') <- numeric(0)
  expect_equal(expectation, PrepareDataProfile(dat))
})

test_that("PrepareDataProfile()", {
  
  # Easy one
  mtx <- cbind(c('0', '0', 1,1,1,1),
               c(0,0,1,1,1,1),# again
               c(0,0,0,1,1,'?'))
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  phy1 <- TreeTools::MatrixToPhyDat(mtx)
  expect_equivalent(phy1, PrepareDataProfile(phy1))
  expect_equal(attributes(phy1), attributes(PrepareDataProfile(phy1))[1:10])
  
  # Easy one
  mtx <- cbind(c('0', '0', 1,1,1,1),
               c(1,1,0,0,0,0),# flipped
               c(0,0,0,1,1,'{012}'))
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  phy2 <- TreeTools::MatrixToPhyDat(mtx)
  expect_equivalent(phy1, PrepareDataProfile(phy2))
  expect_equal(attributes(PrepareDataProfile(phy1)),
               attributes(PrepareDataProfile(phy2)))
  
  
  mtx <- cbind(c('0', '0', 1,1,1, '2', '2', 3,3,3,3),
               c('?', '?', 1,1,1, '?', '?', 0,0,0,0),
               c(0,0,1,1,1,2,2,3,3,3,3),# again
               c(rep('?', 5), '2', '2', 0,0,0,0),
               c('?', '?', 1,1,1, 1,1, 0,0,0,0),
               c('0', '1', rep('?', 9))
               )
  rownames(mtx) <- letters[seq_len(nrow(mtx))]
  dataset <- TreeTools::MatrixToPhyDat(mtx)
  
  q <- '?'
  decomposed <- matrix(c(0,0,q,q,q,q,q,1,1,1,1,
                         q,q,0,0,0,q,q,1,1,1,1,
                         q,q,q,q,q,0,0,1,1,1,1,
                         
                         q,q,0,0,0,q,q,1,1,1,1,
                         
                         0,0,q,q,q,q,q,1,1,1,1,
                         q,q,0,0,0,q,q,1,1,1,1,
                         q,q,q,q,q,0,0,1,1,1,1,
                         
                         q,q,q,q,q,0,0,1,1,1,1,
                         q,q,0,0,0,0,0,1,1,1,1),
                       ncol = 9, dimnames = list(letters[1:11], NULL))
                         
                         
  expect_warning(pd <- PrepareDataProfile(dataset))
  expect_equal(decomposed, PhyDatToMatrix(pd))
  expect_equal(c(1, 2, 3, 2, 1, 2, 3, 3, 4), attr(pd, 'index'))
  expect_equal(c(2, 3, 3, 1), attr(pd, 'weight'))
  
  dataset2 <- TreeTools::MatrixToPhyDat(mtx[!mtx[, 1] %in% c(0, 2), ])
  expect_equal(attr(PrepareDataProfile(dataset2), 'info.amounts'),
               attr(pd, 'info.amounts')[1:3, 2, drop = FALSE])
  
  
  data('Lobo', package = "TreeTools")
  expect_warning(prep <- PrepareDataProfile(Lobo.phy))
  expect_equal(c(17, attr(prep, 'nr')),
               dim(attr(prep, 'info.amounts')))
  
})
