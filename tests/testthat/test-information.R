context("Information.R")

test_that("Entropy is calculated correctly", {
  expect_equal(1, Entropy(rep(0.5, 2)))
  expect_equal(2, Entropy(c(1/4, 1/4, 0, 1/4, 0, 1/4)))
})

test_that("Trees matching splits calculated correctly", {
  expect_equal(NUnrooted(9), TreesMatchingSplit(0, 9))
  expect_equal(LnUnrooted(9), LogTreesMatchingSplit(0, 9))
  expect_equal(LnUnrooted.int(9), LogTreesMatchingSplit(9, 0))
  expect_equal(log(315/10395)/-log(2), SplitInformation(3, 5))
})

test_that("UnrootedTreesMatchingSplit works", {
  expect_equal(NRooted(3) * NRooted(5), UnrootedTreesMatchingSplit(c(3, 5)))
  expect_equal(NRooted(30) * NRooted(50), UnrootedTreesMatchingSplit(c(30, 50)))
})

test_that("Joint information calculated correctly", {
  # Identical splits: ABCDE:FGH, ABCDE:FGH
  expect_equal(-log2(315/10395), JointInformation(5, 0, 0, 3))
  expect_equal(-log2(315/10395), JointInformation(3, 0, 0, 5))
  expect_equal(-log2(315/10395), JointInformation(0, 5, 3, 0))
  expect_equal(-log2(315/10395), JointInformation(0, 3, 5, 0))
  
  # Agreeable splits: ABCDE:FGHI, ABC:DEFGHI
  expected_3204 <- SplitInformation(5, 4) + SplitInformation(3, 6) --log2(135/135135)
  expect_equal(expected_3204, JointInformation(3, 2, 0, 4))
  expect_equal(expected_3204, JointInformation(2, 3, 4, 0))
  expect_equal(expected_3204, JointInformation(0, 4, 3, 2))
  expect_equal(expected_3204, JointInformation(4, 0, 2, 3))
  
  # Perfect contradiction: AB:CDEFG, AC:BDEFG
  expect_equal(SplitInformation(2, 5) * 2, JointInformation(1, 1, 1, 4))
  
  # Compatible splits: AB:CDEFGH, CD:ABEFGH
  expected_0224 <- SplitInformation(2, 6) - SplitInformation(2, 5)
  expect_equal(expected_0224, JointInformation(0, 2, 2, 4))
  expect_equal(expected_0224, JointInformation(2, 0, 4, 2))
  expect_equal(expected_0224, JointInformation(2, 4, 0, 2))
  expect_equal(expected_0224, JointInformation(4, 2, 2, 0))
  
})

test_that("AllSplitPairings counted correctly", {
  expect_error(AllSplitPairings(3))
  for (n in 4:10) {
    totalSplits <- sum(choose(n, 2:(n-2)))
    expect_equal(totalSplits * totalSplits, sum(AllSplitPairings(n)))
  }
})

