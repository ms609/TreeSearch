test_that("PhyDat2Morphy() errors", {
  expect_error(PhyDat2Morphy(NA))
})

test_that("UnloadMorphy() errors", {
  expect_error(UnloadMorphy(NA))
})

test_that("GapHandler()", {
  expect_error(GapHandler(0))
  tokens <- matrix(c("-", "-", 0, 0), byrow = TRUE, nrow = 4L,
                   dimnames = list(letters[1:4], NULL))
  pd <- TreeTools::MatrixToPhyDat(tokens)
  
  morphyObj <- PhyDat2Morphy(pd)
  expect_equal(0, RandomTreeScore(morphyObj))
  expect_equal("Inapplicable", GapHandler(morphyObj))
  UnloadMorphy(morphyObj)
  
  morphyObj <- PhyDat2Morphy(pd, "ambigu")
  expect_equal(0, RandomTreeScore(morphyObj))
  expect_equal("Missing data", GapHandler(morphyObj))
  UnloadMorphy(morphyObj)
  
  morphyObj <- PhyDat2Morphy(pd, "eXt")
  expect_lt(0, RandomTreeScore(morphyObj))
  expect_equal("Extra state", GapHandler(morphyObj))
  UnloadMorphy(morphyObj)
  
  morphyObj <- SingleCharMorphy("-0-0", "eXt")
  expect_lt(0, RandomTreeScore(morphyObj))
  expect_equal("Extra state", GapHandler(morphyObj))
  UnloadMorphy(morphyObj)
  
  expect_error(SingleCharMorphy("-0-0", "ERROR"))
  expect_error(GapHandler(morphyObj))
})
