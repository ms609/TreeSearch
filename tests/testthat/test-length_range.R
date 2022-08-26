test_that("Deprecation", {
  expect_equal(MinimumLength(1:3), expect_warning(MinimumSteps(1:3)))
})

test_that("Minimum step counts are correctly calculated", {
  library("TreeTools")
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