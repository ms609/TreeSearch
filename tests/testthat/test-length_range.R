test_that("Deprecation", {
  expect_equal(MinimumLength(1:3), expect_warning(MinimumSteps(1:3)))
})

test_that("Step counts are correctly calculated", {
  library("TreeTools")
  expect_equal(MinimumLength(rep(1, 3)), 0)
  expect_equal(MaximumLength(rep(1, 3)), 0)
  expect_equal(MinimumLength(1:3), 1)
  expect_equal(MaximumLength(1:3), 1)
  
  # ..+, .+., .++, +.+ 
  expect_equal(MinimumLength(c(1:3, 5)), 1)
  expect_equal(MinimumLength(c(6, 7, 14)), 0)
  expect_equal(MinimumLength(0:3), 1) # 0 representing the inapplicable token
  
  # ++++, .++., ..++
  expect_equal(0, MinimumLength(c(2046, 384, 1152)))
  
  # ++++, +..., .++., ..++
  expect_equal(1, MinimumLength(c(15, 8, 6, 3)))
  
  # ++++++, +....., .++..., .+.+.., ...++.
  expect_equal(2, MinimumLength(c(63, 32, 24, 20, 6)))
  
  dudDat <- TreeTools::StringToPhyDat("----{-,1}22", letters[1:7])
  expect_equal("----<-,1>22", TreeTools::PhyDatToString(dudDat, ">", ","))
  expect_equal(0, attr(PrepareDataIW(dudDat), "min.length"))
  
  dudTwo <- TreeTools::StringToPhyDat("{-1}{-2}{-3}2233", letters[1:7])
  expect_equal("{-1}{-2}{-3}2233", TreeTools::PhyDatToString(PrepareDataIW(dudTwo)))
  
  morphyObj <- SingleCharMorphy("{-1}{-2}{-3}2233")
  expect_equal(MorphyTreeLength(TreeTools::PectinateTree(7), morphyObj), 1)
  morphyObj <- UnloadMorphy(morphyObj)
  
  owch2 <- "{-1}{-2}22{-3}33"
  tr2 <- ape::read.tree(text=("(a, ((b, (c, d)), (e, (f, g))));"))
  # PlotCharacter(tr2, StringToPhyDat(owch2, letters[1:7]))
  
  
  morphyObj <- SingleCharMorphy(owch2)
  expect_equal(MorphyTreeLength(TreeTools::PectinateTree(7), morphyObj), 1)
  morphyObj <- UnloadMorphy(morphyObj)
  
  owch3 <- "-1-222-333"
  tr3 <- ape::read.tree(text=("((a1, a2), (((b1, b2), (c, d)), ((e1, e2), (f, g))));"))
  # PlotCharacter(tr3, StringToPhyDat(owch3, TipLabels(tr3)))
  
  morphyObj <- SingleCharMorphy(owch3)
  expect_equal(MorphyTreeLength(TreeTools::PectinateTree(10), morphyObj), 2)
  expect_equal(MorphyTreeLength(tr3, morphyObj), 2)
  morphyObj <- UnloadMorphy(morphyObj)
  
  expect_equal(2, MinimumLength("-{-1}{-2}{-3}2233"))
  expect_equal(1, MinimumLength("--{-1}{-2}{-3}2233"))
  
  expect_equal(0, attr(PrepareDataIW(dudDat), "min.length"))
  
  tr <- ape::read.tree(text="(((a, b), c), (d, (e, ((f, g), (h, (i, (j, k)))))));")
  expect_equal(CharacterLength(tr, compress = TRUE,
                               TreeTools::StringToPhyDat("11---22--33", letters[1:11])),
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
  
  data("inapplicable.datasets")
  expect_equal(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
                 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 
                 1, 2, 1, 1, 4, 3, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 4, 1, 
                 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 
                 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
               MinimumLength(inapplicable.phyData[[4]], compress = TRUE))
  
})

test_that("MaximumLength() edge cases are handled correctly", {
  expect_equal(MaximumLength("00000000000111"), 3)
  expect_equal(MaximumLength("001122{12}"), 4)
  expect_equal(MaximumLength("0123 0123 0123 ????"), 3 * 3)
  expect_equal(MaximumLength( "00112233{01}{23}{012}?"), 3 + 3 + 1)
})

test_that("MaximumLength() handles inapplicable tokens", {
  # Number of regions = number of inapplicable tokens - 1
  # One extra step allowable for each extra region
  expect_equal(MaximumLength("11111--"), 0)
  
  if (interactive()) {
    tree <- ape::read.tree(
      text = ("(a, (b, (c, (i1, (i2, (i3, e))))));"))
    char <- "1110001"
    charDat <- StringToPhyDat(char, TipLabels(tree))
    PlotCharacter(tree, charDat)
  }
  expect_equal(MaximumLength("111---1"), 1)
  
  if (interactive()) {
    tree <- ape::read.tree(
      text = ("(a, (b, (c, ((i1, (i2, d)), (i3, (i4, e))))));"))
    char <- "111001001"
    charDat <- StringToPhyDat(char, TipLabels(tree))
    PlotCharacter(tree, charDat)
  }
  expect_equal(MaximumLength("111--1--1"), 2)
  
  if (interactive()) {
    tree <- ape::read.tree(
      text = ("(a, (b, (i0, (((i1, d), (i2, e)), ((i3, c), (i4, f))))));"))
    char <- "11001010101"
    charDat <- StringToPhyDat(char, TipLabels(tree))
    PlotCharacter(tree, charDat) # Treated as independent losses
    tree <- ape::read.tree(
      text = ("(a, (b, (c, (i0, (((i1, d), (i2, e)), (i3, (i4, f)))))));"))
    char <- "11100101001"
    charDat <- StringToPhyDat(char, TipLabels(tree))
    PlotCharacter(tree, charDat) # Can support four distinct regions
  }
  expect_equal(MaximumLength("111--1---1"), 3)
  expect_equal(MaximumLength("111-----"), 2)
  expect_equal(MinimumLength("011-----"), 1)
  expect_equal(MaximumLength("011-----"), 2)
  expect_equal(MaximumLength("0111-----"), 3)
  twoPair <- StringToPhyDat("0011-----", tips = letters[1:9])
  colnames(attr(twoPair, "contrast")) <- NULL
  expect_equal(MaximumLength(twoPair), 3)
  expect_equal(MaximumLength("1--1---0"), 2)
  expect_equal(MaximumLength("--1---1"), 1)
  expect_equal(MinimumLength("--1---0"), 1)
  expect_equal(MinimumLength("--1---0"), 1)
  expect_equal(MaximumLength("--1---"), 0)
  expect_equal(MaximumLength("-----"), 0)
  
  expect_equal(MaximumLength("001122{12}---"), 4 + 1)
  expect_equal(MaximumLength("0123 0123 0123 ----"), 3 * 3 + 2)
  expect_equal(MaximumLength("00112233{01}{23}{012}----"), 3 + 3 + 1 + 2)
})

test_that("MaximumLength() handles many states (>8) without overflow", {
  # Regression test for crash with datasets having more than 8 character states
  # When nToken > 255, as.raw() would overflow; now uses bitwAnd/bitwOr instead
  
  # Test with 10 states (nToken = 1023)
  # This should complete without error (the specific result is secondary)
  # 10 distinct character states (powers of 2 for bit representation)
  manyStates <- c(1, 2, 4, 8, 16, 32, 64, 128, 256, 512)
  expect_silent(result <- MaximumLength.numeric(manyStates))
  expect_equal(result, 9)  # max steps for 10 distinct tokens
  
  # Test with combination of states including inapplicable (0)
  manyStatesWithInapp <- c(0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 0)
  expect_silent(result <- MaximumLength.numeric(manyStatesWithInapp))
  expect_equal(result, 9)  # 9 steps + 0 extra for inapplicable regions
  
  # Test with large combined token value (all 10 states present in one tip)
  # This is the scenario that triggers the crash with Vinther2008 dataset
  # 1 = only state 0, 1023 = all states 0-9 combined (2^10 - 1)
  largeToken <- c(1, 1023)
  expect_silent(result <- MaximumLength.numeric(largeToken))
  expect_equal(result, 0)  # state 0 is subset of 1023, so no extra steps needed
})

test_that("MaximumLength() works with Vinther2008 dataset", {
  # Regression test for the original crash report
  # This dataset has 10 states (0-9) which caused overflow
  data("inapplicable.datasets")
  expect_silent(result <- MaximumLength(inapplicable.phyData[["Vinther2008"]]))
  expect_true(is.numeric(result))
  expect_true(all(result >= 0))
})



