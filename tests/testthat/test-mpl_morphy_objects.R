test_that("PhyDat2Morphy() errors", {
  suppressWarnings(expect_error(PhyDat2Morphy(NA)))
})

test_that("UnloadMorphy() errors", {
  suppressWarnings(expect_error(UnloadMorphy(NA)))
})

test_that("GapHandler()", {
  expect_error(GapHandler(0))
  tokens <- matrix(c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
                   dimnames = list(letters[1:4], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  
  morphyObj <- suppressWarnings(PhyDat2Morphy(pd))
  expect_equal(0, RandomTreeScore(morphyObj))
  expect_equal("Inapplicable", GapHandler(morphyObj))
  suppressWarnings(UnloadMorphy(morphyObj))
  
  morphyObj <- suppressWarnings(PhyDat2Morphy(pd, "ambigu"))
  expect_equal(0, RandomTreeScore(morphyObj))
  expect_equal("Missing data", GapHandler(morphyObj))
  suppressWarnings(UnloadMorphy(morphyObj))
  
  morphyObj <- suppressWarnings(PhyDat2Morphy(pd, "eXt"))
  expect_lt(0, RandomTreeScore(morphyObj))
  expect_equal("Extra state", GapHandler(morphyObj))
  suppressWarnings(UnloadMorphy(morphyObj))
  
  morphyObj <- SingleCharMorphy("-0-0", "eXt")
  expect_lt(0, RandomTreeScore(morphyObj))
  expect_equal("Extra state", GapHandler(morphyObj))
  UnloadMorphy(morphyObj)
  
  expect_error(SingleCharMorphy("-0-0", "ERROR"))
  suppressWarnings(expect_error(GapHandler(morphyObj)))
})

test_that("morphy_profile fails nicely", {
  morphyObj <- SingleCharMorphy("1")
  on.exit(UnloadMorphy(morphyObj))
  expect_error(TreeSearch:::morphy_profile(matrix(NA, 10, 2), list(morphyObj),
                              1, 1L, matrix(1), 1),
               "Number of edges does not match Morphy object dimensions")
})
