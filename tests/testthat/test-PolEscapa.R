test_that("LengthAdded() errors", {
  trees <- inapplicable.trees[["Vinther2008"]]
  dataset <- inapplicable.phyData[["Vinther2008"]]
  
  expect_error(
    LengthAdded(trees, "dataset"),
    "`char` must be a character of class `phyDat`"
  )
  
  expect_error(
    LengthAdded(trees, dataset),
    "`char` must comprise a single character"
  )
  
  # Error when a used token has a zero-sum contrast row
  char51 <- dataset[, 51]
  usedTokens <- unique(unlist(char51, use.names = FALSE))
  attr(char51, "contrast")[usedTokens[[1L]], ] <- 0
  expect_error(
    LengthAdded(trees, char51),
    "`char` contrast matrix lacks levels for token"
  )

  # No error when only unused tokens have zero-sum contrast rows; also
  # verifies that downstream scoring does not choke on the stale row
  char51b <- dataset[, 51]
  cont51b <- attr(char51b, "contrast")
  usedTokens2 <- unique(unlist(char51b, use.names = FALSE))
  unusedRows <- setdiff(seq_len(nrow(cont51b)), usedTokens2)
  if (length(unusedRows) > 0L) {
    attr(char51b, "contrast")[unusedRows[[1L]], ] <- 0
    expect_no_error(LengthAdded(trees, char51b))
  }
})

test_that("LengthAdded()", {
  trees <- inapplicable.trees[["Vinther2008"]]
  dataset <- inapplicable.phyData[["Vinther2008"]]
  
  pe10 <- LengthAdded(trees, dataset[, 10])
  expect_equal(pe10["Neopilina"], c(Neopilina = 1))
  expect_equal(sum(pe10), 1)
  
  # Single tree
  expect_equal(LengthAdded(trees[[1]], dataset[, 10]), pe10)
  
  # No inapplicables
  appData <- dataset
  colnames(attr(appData, "contrast"))[1] <- "x"
  attr(appData, "levels")[1] <- "x"
  attr(appData, "allLevels")[4] <- "x"
  pe10 <- LengthAdded(trees, appData[, 10])
  expect_equal(pe10["Neopilina"], c(Neopilina = 1))
  expect_equal(sum(pe10), 1)
  
  # Implied weighting
  expect_equal(
    unname(PolEscapa(trees, dataset[, 11], concavity = 5)["Neopilina"]),
    as.numeric(TreeLength(trees[[1]], dataset[, 11], concavity = 5))
  )
  
  # minLength changes when only occurrence of 1 -> ?
  wiwaxia <- LengthAdded(trees, dataset[, 39], concavity = 10)
  expect_true(all(wiwaxia >= 0))
  
})
