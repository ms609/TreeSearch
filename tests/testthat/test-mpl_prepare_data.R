context('mpl_prepare_data.R')
test_that('PhyToString works', {
  longLevels <- phyDat(rbind(x = c('-', '?', 0:12), y = c(12:0, '-', '?')), 
                       type='USER', levels=c(0:12, '-'))
  expect_equal("-?0123456789ABCCBA9876543210-?", PhyToString(longLevels))
})