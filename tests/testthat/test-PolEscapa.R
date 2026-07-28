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

test_that("LengthAdded() qm-empty: no fully-ambiguous contrast row", {
  # Regression: when a character has a {-, state} (partial-inapplicable)
  # ambiguity but NO fully ambiguous ("?") contrast row, `qm` was integer(0).
  # A leaf whose starting token is inapplicable then hit
  # `charQm[[leaf]] <- qm`, assigning integer(0) and silently corrupting the
  # phyDat (dropping an element) — surfacing as a recycling warning and a
  # wrong instability score.  The fix appends an all-ones fallback row.
  skip_if_not_installed("phangorn")
  nTips     <- 6L
  tipLabels <- paste0("t", seq_len(nTips))
  # t3 inapplicable ("-"); t6 has the {-,0} partial ambiguity (app AND inapp).
  tipCodes  <- c("0", "1", "-", "0", "1", "{-0}")
  names(tipCodes) <- tipLabels

  # No "?" (all-ones) row exists, so `qm` is empty pre-fix.
  levs <- c("-", "0", "1", "{-0}")
  cont <- matrix(
    c(1, 0, 0,   # "-"
      0, 1, 0,   # "0"
      0, 0, 1,   # "1"
      1, 1, 0),  # "{-0}" — applicable AND inapplicable
    nrow = 4L, ncol = 3L, byrow = TRUE,
    dimnames = list(levs, c("-", "0", "1"))
  )
  char <- phangorn::phyDat(
    setNames(as.list(tipCodes), tipLabels),
    type = "USER", levels = levs, contrast = cont
  )
  set.seed(11L)
  trees <- c(TreeTools::RandomTree(char, root = TRUE))

  # Must not warn (phyDat recycling) and must return one finite, non-negative
  # value per tip.
  expect_no_warning(result <- LengthAdded(trees, char, concavity = 10))
  expect_length(result, nTips)
  expect_false(anyNA(result))
  expect_true(all(result >= 0))

  # EW path likewise.
  expect_no_warning(resultEw <- LengthAdded(trees, char))
  expect_length(resultEw, nTips)
  expect_false(anyNA(resultEw))
})


# Profile parsimony: regressions for T-365 ------------------------------------
# `PrepareDataProfile()` replaces the contrast matrix wholesale with
# `rbind(diag(k), rep(1, k))` and renumbers every token to `1:k`, with `k + 1`
# denoting ambiguity.  `LengthAdded()` used to read `cont`, `qm`, `qmApp`,
# `app` and `inapp` from the character the *user* supplied, so those row and
# token indices named rows of a contrast that no longer existed.  Every
# character whose raw and prepared token spaces differed in shape -- i.e.
# everything but a plain {0, 1, ?} character -- errored out.

# A single character with one token per leaf.
ProfileChar <- function (tokens) {
  TreeTools::MatrixToPhyDat(matrix(
    tokens, ncol = 1,
    dimnames = list(paste0("t", seq_along(tokens)), NULL)
  ))
}

# Two trees that put the character in conflict, so an informative leaf has a
# non-zero instability and the tests cannot pass on an all-zero result.
ProfileTrees <- function (tipOrder) {
  c(TreeTools::PectinateTree(tipOrder), TreeTools::BalancedTree(tipOrder))
}

# Recompute the expected instability independently of `LengthAdded()`, in the
# token space of the *prepared* character: ambiguity is its last contrast row.
ProfileExpectation <- function (trees, char) {
  prepared <- suppressMessages(PrepareDataProfile(char))
  ambiguous <- nrow(attr(prepared, "contrast"))
  rooted <- RootTree(trees, 1)
  start <- TreeLength(rooted, prepared, "profile")
  vapply(names(char), function (leaf) {
    ambiguated <- prepared
    ambiguated[[leaf]] <- ambiguous
    mean(start - TreeLength(rooted, ambiguated, "profile"))
  }, double(1))
}

test_that("LengthAdded(concavity = 'profile') collapses singleton states", {
  # Tokens 0,0,0,1,1,1,2,2,3,?: the raw contrast is 5x4, but `3` occurs once,
  # so profile parsimony treats it as ambiguous and the prepared contrast is
  # 4x3.  The raw `?` token (5) then indexes past the prepared contrast:
  # pre-fix, "`tip_data` values must be in [1, nrow(contrast)] (4); found 5".
  char <- ProfileChar(c("0", "0", "0", "1", "1", "1", "2", "2", "3", "?"))
  trees <- ProfileTrees(c("t1", "t4", "t7", "t2", "t5", "t8", "t3", "t6",
                          "t9", "t10"))

  expect_no_error(added <- LengthAdded(trees, char, concavity = "profile"))
  expect_named(added, names(char))
  expect_false(anyNA(added))
  expect_true(all(added >= 0))

  # The singleton `3` and the `?` are ambiguous once prepared, so ambiguating
  # them cannot change a tree length.
  expect_equal(unname(added[c("t9", "t10")]), c(0, 0))
  # Leaves that do carry profile information must move the score.
  expect_gt(added[["t1"]], 0)

  expect_equal(added, ProfileExpectation(trees, char))
})

test_that("LengthAdded(concavity = 'profile') handles inapplicable tokens", {
  # Tokens 0,0,0,1,1,1,-,-: the raw contrast has three columns (-, 0, 1), the
  # prepared one two.  `qmApp` was empty, so the fallback wrote the raw
  # three-column contrast onto a two-level character: pre-fix,
  # "`levels` length (2) must equal ncol(contrast) (3)".
  char <- ProfileChar(c("0", "0", "0", "1", "1", "1", "-", "-"))
  trees <- ProfileTrees(c("t1", "t4", "t2", "t5", "t3", "t6", "t7", "t8"))

  messages <- testthat::capture_messages(
    added <- LengthAdded(trees, char, concavity = "profile")
  )
  expect_true(any(grepl("Inapplicable tokens treated as ambiguous", messages,
                        fixed = TRUE)))

  expect_named(added, names(char))
  expect_false(anyNA(added))
  expect_true(all(added >= 0))

  # Profile parsimony folds `-` into `?` before scoring, so an inapplicable
  # leaf is ambiguous either way and adds no length -- the same zero the
  # applicability-preserving equal-weights path reports for it.
  expect_equal(unname(added[c("t7", "t8")]), c(0, 0))
  expect_gt(added[["t1"]], 0)

  expect_equal(added, ProfileExpectation(trees, char))
})

test_that("LengthAdded(concavity = 'profile') handles a missing `?` token", {
  # With no fully ambiguous token in the raw character, both the `qmApp`
  # fallback and the `qm` fallback `rbind()` a row onto the contrast.  Pre-fix
  # that row was of raw width, widening the prepared contrast beyond its
  # levels.  Tokens 0,0,0,1,1,1,2,2,3 (three informative states) gave
  # "`levels` length (3) must equal ncol(contrast) (4)".
  char <- ProfileChar(c("0", "0", "0", "1", "1", "1", "2", "2", "3"))
  trees <- ProfileTrees(c("t1", "t4", "t7", "t2", "t5", "t8", "t3", "t6",
                          "t9"))

  expect_no_error(added <- LengthAdded(trees, char, concavity = "profile"))
  expect_named(added, names(char))
  expect_false(anyNA(added))
  expect_true(all(added >= 0))
  expect_equal(unname(added[["t9"]]), 0)
  expect_gt(added[["t1"]], 0)
  expect_equal(added, ProfileExpectation(trees, char))

  # Two informative states and two singletons: the prepared contrast is 3x2,
  # so the raw four-column fallback row gave
  # "`levels` length (2) must equal ncol(contrast) (4)".
  char2 <- ProfileChar(c("0", "0", "0", "1", "1", "1", "2", "3"))
  trees2 <- ProfileTrees(c("t1", "t4", "t2", "t5", "t3", "t6", "t7", "t8"))

  expect_no_error(added2 <- LengthAdded(trees2, char2, concavity = "profile"))
  expect_named(added2, names(char2))
  expect_false(anyNA(added2))
  expect_true(all(added2 >= 0))
  expect_equal(unname(added2[c("t7", "t8")]), c(0, 0))
  expect_gt(added2[["t1"]], 0)
  expect_equal(added2, ProfileExpectation(trees2, char2))
})

test_that("LengthAdded(concavity = 'profile') handles an uninformative character", {
  # State `1` is a singleton, so `maxInformative < 2` and `PrepareDataProfile()`
  # returns a zero-character phyDat (T-372).  `QMScore()` then indexed
  # `char[[leaf]]` expecting one token per taxon, which is now empty for every
  # leaf -- pre-fix: "argument is of length zero" inside `if (!app[startToken])`.
  # Ambiguating a leaf of an already-uninformative character cannot create
  # information, so every delta must be exactly zero; `ProfileExpectation()`
  # is not used here, as it assumes a retained character to overwrite.
  char <- ProfileChar(c("0", "0", "0", "0", "1"))
  trees <- ProfileTrees(paste0("t", 1:5))

  expect_no_error(added <- LengthAdded(trees, char, concavity = "profile"))
  expect_named(added, names(char))
  expect_equal(unname(added), rep(0, length(char)))
})
