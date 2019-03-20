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
  expect_equal(SplitMutualInformation(8, 5, 5), FullMutualInformation(0, 5, 3, 0))
  # Agreeable splits: ABCDE:FGHI, ABC:DEFGHI
  expect_equal(SplitInformation(5, 4) + SplitInformation(3, 6) --log2(135/135135),
               JointInformation(3, 2, 0, 4))
  expect_equal(SplitMutualInformation(9, 5, 3),
               FullMutualInformation(3, 2, 0, 4))
  # Perfect contradiction: AB:CDEFG, AC:BDEFG
  expect_equal(SplitInformation(2, 5) * 2, JointInformation(1, 1, 1, 4))
  expect_equal(0, FullMutualInformation(1, 1, 1, 4))
  # Contradiction with common information: AB CD : EF GHI, AB EF : CD GHI
  expect_equal(SplitInformation(2, 2) + SplitInformation(2, 3),
               FullMutualInformation(2, 2, 2, 3))
})

test_that("SplitMatchProbability returns expected probabilities", {

  splitAB   <- c(rep(TRUE, 2), rep(FALSE, 7))
  splitABC  <- c(rep(TRUE, 3), rep(FALSE, 6))
  splitABI  <- c(rep(TRUE, 2), rep(FALSE, 6), TRUE)
  splitBCD  <- c(FALSE, rep(TRUE, 3), rep(FALSE, 5))
  splitAEF  <- c(TRUE, rep(FALSE, 3), rep(TRUE, 2), rep(FALSE, 3))
  splitABCD <- c(rep(TRUE, 4), rep(FALSE, 5))
  splitABCE <- c(rep(TRUE, 3), FALSE, TRUE, rep(FALSE, 4))
  splitCDEF <- c(rep(FALSE, 2), rep(TRUE, 4), rep(FALSE, 3))
  splitABEF <- c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 2), rep(FALSE, 3))
  
  splitAI <- c(TRUE, rep(FALSE, 7), TRUE)
  splitBC <- c(FALSE, TRUE, TRUE, rep(FALSE, 6))
  splitCD <- c(FALSE, FALSE, TRUE, TRUE, rep(FALSE, 5))
  
  # Possible matches to ABCD:....
  expect_true(SplitMatchProbability(splitABC, splitABCD) < 
                SplitMatchProbability(splitAB, splitABCD))
  
  expect_true(SplitMatchProbability(splitABCD, splitABCD) < 
                SplitMatchProbability(splitABC, splitABCD))
  
  expect_true(SplitMatchProbability(splitBCD, splitABCD) ==
                SplitMatchProbability(splitABC, splitABCD))
  
  expect_true(SplitMatchProbability(splitABC, splitABCD) < 
                SplitMatchProbability(splitAEF, splitABCD))
  
  expect_true(SplitMatchProbability(splitABC, splitABCD) < 
                SplitMatchProbability(splitABI, splitABCD))
  
  expect_true(SplitMatchProbability(splitABI, splitABCD) < 
                SplitMatchProbability(splitAEF, splitABCD))
  
  expect_true(SplitMatchProbability(splitABC, splitABCD) < 
                SplitMatchProbability(splitABCE, splitABCD))
  
  expect_true(SplitMatchProbability(splitCDEF, splitABCD) ==
                SplitMatchProbability(splitABEF, splitABCD))
  
  expect_true(SplitMatchProbability(splitABCE, splitABCD) < 
                SplitMatchProbability(splitABEF, splitABCD))
  
  # Two splits of AB:...
  expect_true(SplitMatchProbability(splitAB, splitAB) < 
              SplitMatchProbability(splitCD, splitAB))
  expect_true(SplitMatchProbability(splitCD, splitAB) < 
              SplitMatchProbability(splitAI, splitAB))
  
  
  Test <- function (score, split1, split2) {
    expect_equivalent(score, SplitMatchProbability(split1, split2))
    expect_equivalent(score, SplitMatchProbability(split2, split1))

    expect_equivalent(score, SplitMatchProbability(split1, !split2))
    expect_equivalent(score, SplitMatchProbability(split2, !split1))
  
    expect_equivalent(score, SplitMatchProbability(!split1, !split2))
    expect_equivalent(score, SplitMatchProbability(!split2, !split1))
    
    expect_equivalent(score, SplitMatchProbability(!split1, split2))
    expect_equivalent(score, SplitMatchProbability(!split2, split1))
    
    score
  }
  
  Test(1, splitAB, splitAI)
  Test(1, splitAB, splitBC)
  Test(1, splitBC, splitCD)
  Test(1, rev(splitAB), splitAI)
  
  Test(1/36, splitAB, splitAB)
  Test(1/36, splitBC, splitBC)
  
  Test(1/126, splitABCD, splitABCD)
  Test(1/126, splitABEF, splitABEF)
  Test(1/126, splitCDEF, splitCDEF)
  Test(1, splitABCD, splitABEF)
  
  Test(1/12, splitAB, splitABC)
  Test(1/12, splitBC, splitABC)
  Test(1/12, splitBC, splitBCD)
  
  Test(1, splitAEF, splitABCD)
  Test(66 / 126, splitABC, splitABEF)
  Test(66 / 126, splitBCD, splitCDEF)
  
  Test(4/84, splitABC, splitABCD)
  Test(4/84, splitBCD, splitABCD)
  
})

test_that("AllSplitPairings counted correctly", {
  expect_error(AllSplitPairings(3))
  for (n in 4:10) {
    totalSplits <- sum(choose(n, 2:(n-2)))
    expect_equal(totalSplits * totalSplits, sum(AllSplitPairings(n)))
  }
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
