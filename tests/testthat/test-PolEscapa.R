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

  # Error message names the specific token index.
  # Asher2005 char 67 has 6 used tokens (1..6); zeroing row 6 must produce
  # a message that explicitly says "6".  The contrast check fires before any
  # tree/data compatibility check, so mismatched trees are harmless here.
  char6tok <- inapplicable.phyData[["Asher2005"]][, 67]
  attr(char6tok, "contrast")[6L, ] <- 0
  expect_error(
    LengthAdded(trees, char6tok),
    "`char` contrast matrix lacks levels for token.s. 6"
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

test_that("LengthAdded() qmApp scalar-unwrap: ≥2 fully-applicable-ambiguous rows", {
  # Regression for T-302/T-305: when ≥2 contrast rows satisfy
  # apply(contApp == 1, 1, all) & !inapp, the pre-fix code left qmApp as a
  # vector; assigning it to charQm[[leaf]] silently corrupted the phyDat
  # passed to TreeLength(), producing large negative deltas.
  skip_if_not_installed("phangorn")
  nTips     <- 8L
  tipLabels <- paste0("t", seq_len(nTips))
  tipCodes  <- c("-", "0", "1", "0", "-", "1", "0", "1")
  names(tipCodes) <- tipLabels

  # Levels "{01}" and "{01}dup" both satisfy apply(contApp == 1, 1, all) &
  # !inapp: two qmApp candidates where the pre-fix code produced a vector.
  levs <- c("-", "0", "1", "?", "{01}", "{01}dup")
  cont <- matrix(
    c(1, 0, 0,
      0, 1, 0,
      0, 0, 1,
      1, 1, 1,   # "?" — qm; excluded from qmApp because inapp = TRUE
      0, 1, 1,   # "{01}" — first qmApp candidate
      0, 1, 1),  # "{01}dup" — second candidate; pre-fix qmApp was c(5, 6)
    nrow = 6L, ncol = 3L, byrow = TRUE,
    dimnames = list(levs, c("-", "0", "1"))
  )
  char <- phangorn::phyDat(
    setNames(as.list(tipCodes), tipLabels),
    type = "USER", levels = levs, contrast = cont
  )

  set.seed(42L)
  trees <- c(TreeTools::RandomTree(char, root = TRUE))

  result <- LengthAdded(trees, char)

  # All deltas non-negative (violated pre-fix when multiple qmApp rows existed).
  expect_true(all(result >= 0))

  # Independent check: tip t3 (coded "1") — manually set to the first qmApp
  # token (row 5 = "{01}") and verify the reported delta matches.
  start  <- TreeLength(trees, char)
  charQm <- char
  charQm[["t3"]] <- 5L
  expect_equal(unname(result[["t3"]]),
               unname(start - TreeLength(trees, charQm)))
})

test_that("LengthAdded() qm scalar-unwrap: ≥2 fully-ambiguous contrast rows", {
  # Regression for the analogous qm fix (commit e8b318c3): when ≥2 rows have
  # rowSums(cont) == ncol(cont), the pre-fix code left qm as a vector.
  # Assigning it to charQm[[leaf]] produced a wildly wrong TreeLength() result
  # and large negative deltas for tips coded with the "?" (fully ambiguous) token.
  skip_if_not_installed("phangorn")
  nTips     <- 8L
  tipLabels <- paste0("t", seq_len(nTips))
  # t4 and t5 are coded "?"; t1 and t6 are inapplicable; rest are applicable.
  tipCodes  <- c("-", "0", "1", "?", "?", "-", "0", "1")
  names(tipCodes) <- tipLabels

  # "?" and "also?" both have rowSums == ncol(cont): two qm candidates.
  levs <- c("-", "0", "1", "?", "also?")
  cont <- matrix(
    c(1, 0, 0,
      0, 1, 0,
      0, 0, 1,
      1, 1, 1,   # "?" — used by t4, t5; first qm candidate
      1, 1, 1),  # "also?" — unused; second candidate; pre-fix qm was c(4, 5)
    nrow = 5L, ncol = 3L, byrow = TRUE,
    dimnames = list(levs, c("-", "0", "1"))
  )
  char <- phangorn::phyDat(
    setNames(as.list(tipCodes), tipLabels),
    type = "USER", levels = levs, contrast = cont
  )

  set.seed(7L)
  trees <- c(TreeTools::RandomTree(char, root = TRUE))

  result <- LengthAdded(trees, char)

  # All deltas non-negative (violated pre-fix when qm was a vector).
  expect_true(all(result >= 0))

  # Tips coded "?" are already fully ambiguous: setting them to qm (= "?")
  # leaves tree length unchanged, so their deltas must be exactly 0.
  # Independently verify against manual TreeLength() calls.
  start <- TreeLength(trees, char)
  for (tip in c("t4", "t5")) {
    charQm <- char
    charQm[[tip]] <- 4L   # row 4 = "?" = qm[[1L]] after scalar-unwrap
    expect_equal(unname(result[[tip]]),
                 unname(start - TreeLength(trees, charQm)))
    expect_equal(unname(result[[tip]]), 0)
  }
})
