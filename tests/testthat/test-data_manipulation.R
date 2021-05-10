context("data_manipulation.R")

test_that("Deprecation", {
  expect_equal(MinimumLength(1:3), expect_warning(MinimumSteps(1:3)))
})

test_that("Minimum step counts are correctly calculated", {
  expect_equal(1, MinimumLength(1:3))
  expect_equal(1, MinimumLength(c(1:3, 5)))
  expect_equal(0, MinimumLength(c(6, 7, 14)))
  expect_equal(1, MinimumLength(0:3)) # 0 representing the inapplicable token
  
  # ++++, .++., ..++
  expect_equal(0, MinimumLength(c(2046, 384, 1152)))
  
  # ++++, +..., .++., ..++
  expect_equal(1, MinimumLength(c(15, 8, 6, 3)))
  
  # ++++++, +....., .++..., .+.+.., ...++.
  expect_equal(2, MinimumLength(c(63, 32, 24, 20, 6)))
  
  dudDat <- TreeTools::StringToPhyDat('----{-,1}22', letters[1:7])
  expect_equal('----<-,1>22', TreeTools::PhyDatToString(dudDat, '>', ','))
  expect_equal(0, attr(PrepareDataIW(dudDat), 'min.length'))
  
  dudTwo <- TreeTools::StringToPhyDat('{-1}{-2}{-3}2233', letters[1:7])
  expect_equal('{-1}{-2}{-3}2233', TreeTools::PhyDatToString(PrepareDataIW(dudTwo)))
  
  tr <- ape::read.tree(text='(((a, b), c), (d, (e, ((f, g), (h, (i, (j, k)))))));')
  expect_equal(CharacterLength(tr,
                               TreeTools::StringToPhyDat('11---22--33', letters[1:11])),
               MinimumLength(c(0, 0, 0, 0, 0, 0, 2, 2, 4, 4, 8, 8)))

  # 04, 14, 24, 34, 05, 16, 27, 38, 9A
  # In this case, chosing the most common state (4) means that we have to choose 567&8 too
  # 012&3 is a better solution
  # We also have to choose one of 9 or A, but it doesn't matter which.
  expect_equal(4, MinimumLength(c(
    2^0 + 2^4,
    2^1 + 2^4,
    2^2 + 2^4,
    2^3 + 2^4,
    2^0 + 2^5,
    2^1 + 2^6,
    2^2 + 2^7,
    2^3 + 2^8,
    2^9 + 2^10
  )))
  
  data('inapplicable.datasets')
  expect_equal(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 
                 1, 2, 1, 1, 4, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 
                 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
               MinimumLength(inapplicable.phyData[[4]]))
  
})

test_that("PrepareDataProfile() handles empty matrices", {
  dat <- TreeTools::MatrixToPhyDat(matrix(c(0, 1, rep('?', 5)),
                                          dimnames = list(letters[1:7], NULL)))
  expect_warning(expect_error(PrepareDataProfile(dat)))
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
