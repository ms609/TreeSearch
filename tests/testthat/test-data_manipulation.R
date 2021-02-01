context("data_manipulation.R")
test_that("Minimum step counts are correctly calculated", {
  expect_equal(1, MinimumLength(1:3))
  expect_equal(1, MinimumLength(c(1:3, 5)))
  expect_equal(0, MinimumLength(c(6, 7, 14)))
  expect_equal(1, MinimumLength(0:3)) # 0 representing the inapplicable token
  
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
  
})
