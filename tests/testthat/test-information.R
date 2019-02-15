context("Information.R")

test_that("Information functins work", {
  expect_equal(NUnrooted(9), TreesMatchingSplit(0, 9))
  expect_equal(LnUnrooted(9), LogTreesMatchingSplit(0, 9))
  expect_equal(LnUnrooted.int(9), LogTreesMatchingSplit(9, 0))
  expect_equal(log(315/10395)/-log(2), SplitInformation(3, 5))
  expect_equal(JointInformation(8, 3, 3), SplitInformation(3, 5))
  expect_equal(-log2(315/10395), JointInformation(8, 3, 0))
  expect_equal(-log2(135/135135), JointInformation(9, 5, 3))
})

test_that("UnrootedTreesMatchingSplit works", {
  expect_equal(NRooted(3) * NRooted(5), UnrootedTreesMatchingSplit(c(3, 5)))
  expect_equal(NRooted(30) * NRooted(50), UnrootedTreesMatchingSplit(c(30, 50)))
})

test_that("TreesConsistentWithTwoSplits works", {
  
  Test <- function (n, a, b, score) {
    logScore <- log(score)
    
    expect_equal(score, TreesConsistentWithTwoSplits(n, a, b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, b, a))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(score, TreesConsistentWithTwoSplits(n, n - b, n - a))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, a, b))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, b, a))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, n - a, n - b))
    expect_equal(logScore, LogTreesConsistentWithTwoSplits(n, n - b, n - a))
  }
  
  Test(8, 3, 0, 315)
  Test(8, 8, 0, NUnrooted(8))
  Test(8, 3, 3, TreesMatchingSplit(3, 5))
  Test(10, 5, 2, 1575)
  Test(9, 5, 3, 135)
  Test(8, 7, 3, 315)
})
