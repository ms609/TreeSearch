library("TreeTools", quietly = TRUE)

test_that("_Concordance() handles null input", {
  expect_warning(expect_null(QuartetConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(PhylogeneticConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(ClusteringConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(
    MutualClusteringConcordance(BalancedTree(8), NULL)))
  expect_warning(expect_null(
    SharedPhylogeneticConcordance(BalancedTree(8), NULL)))
})

test_that("_Concordance() handles tip mismatch", {
  char <- MatrixToPhyDat(cbind(c(a = 0, b = 0, c = 0, d = 1, e = 1)))
  tree <- BalancedTree(5)
  expect_warning(expect_null(QuartetConcordance(tree, char)),
                 "No overlap between tree labels and dataset.")
  expect_warning(expect_null(ClusteringConcordance(tree, char)),
                 "Could not find 't1', .* in names.dataset.")
})

test_that("QuartetConcordance() works", {
  # This block locks the raw-quartet contract, so it names `unit` and
  # `chanceCorrect` explicitly rather than relying on the (NRQS) defaults.
  Naive <- function(...) {
    QuartetConcordance(..., unit = "quartet", chanceCorrect = FALSE)
  }
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, 0, 0, 1, 1, 1, 1,  0,
                      0, 1, 0, 1, 0, 1, 0, 1,  0,
                      0, 0, 0, 1, 0, 1, 1, 1,  0,
                      0, 0, 0, 0, 1, 1, 2, 2,  0,
                      0, 0, 1, 1, 2, 2, 3, 3,  0,
                      0, 1, 2, 3, 0, 1, 2, 3,  0), 9,
                    dimnames = list(paste0("t", 1:9), NULL))
  expect_error(QuartetConcordance(tree, mataset),
               "`dataset` must be a phyDat object")
  dat <- MatrixToPhyDat(mataset)
  expect_equal(unname(Naive(tree, dat[, 1])), rep(1, 5))
  # plot(tree); nodelabels();
  expect_equal(Naive(tree, dat[, 2]),
               c("10" = 1/9, "11" = 0, "12" = 0,
                 "13" = 1/9, "14" = 0, "15" = 0)[names(as.Splits(tree))])
  
  allQuartets <- combn(8, 4)
  for (charI in seq_len(ncol(mataset))) {
    qc <- Naive(tree, dat[, charI])
    for (splitI in seq_along(splits)) {
      split <- splits[[splitI]]
      logiSplit <- as.logical(split)
      case <- apply(allQuartets, 2, function (q) {
        qSplit <- logiSplit[q]
        qChar <- mataset[q, charI]
        if (identical(unique(table(qSplit)), 2L) &&
            identical(unique(table(qChar)), 2L)) {
          tbl <- table(qSplit, qChar)
          tab <- paste0(sort(tbl[tbl > 0]), collapse = "")
          switch(tab,
                 "1111" = FALSE,
                 "112" = NA,
                 "13" = NA,
                 "22" = TRUE,
                 "4" = NA,
                 stop(q, ": ", tab)
          )
        } else {
          NA
        }
      })
      expect_equal(sum(case, na.rm = TRUE) / sum(!is.na(case)),
                   unname(qc[as.character(names(split))]))
    }
  }
  
  expect_equal(Naive(tree, dat[, c(1:4, 6)]),
               c("10" = (36 + 2 + 9 + 12) / (36 + 18 + 18 + 12 + 6),
                 "11" = ( 6 + 0 + 6 +  2) / ( 6 +  9 +  6 +  2 + 1),
                 "12" = ( 6 + 0 + 0 +  2) / ( 6 +  9 +  9 +  2 + 1),
                 "13" = (36 + 2 + 9 + 12) / (36 + 18 + 18 + 12 + 6),
                 "14" = ( 6 + 0 + 0 +  7) / ( 6 +  9 +  9 +  7 + 1),
                 "15" = ( 6 + 0 + 6 +  7) / ( 6 +  9 +  6 +  7 + 1))[
                   names(as.Splits(tree))]
  )
})

test_that("QuartetConcordance() handles ambiguity", {
  library("TreeTools", quietly = TRUE)
  tree <- BalancedTree(12)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, "{01}", 0, 0, "{01}", 1, 1, "-", 1, 1, "-",
                      0, 1, "?", 0, 1, "?", 0, 1, "(01)", 0, 1, "(01)",
                      0, 0, "?", 0, 1, "(12)", 0, 1, "(12)", 1, 1, "(12)",
                      0, 0, "?", 0, 0, "?", 1, 1, "?", 2, 2, "?",
                      0, 0, "?", 0, 0, "?", 0, 0, "-", 0, 0, "-",
                      rep("?", 12),
                      0, 1, "?", 2, 3, "?", 0, 1, "-", 2, 3, "-"), 12,
                    dimnames = list(paste0("t", 1:12), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  expectation <- unname(QuartetConcordance(tree, dat)[
    c("14", "16", "18", "19", "21", "23")])
  expect_equal(
    unname(QuartetConcordance(DropTip(tree, paste0("t", 3 * 1:4)), dat)),
    expectation[!is.na(expectation)]
  )
  
  expectation <- unname(QuartetConcordance(tree, dat)[
    c("14", "15", "17", "19", "20", "22")])
  expect_equal(
    unname(QuartetConcordance(DropTip(tree, paste0("t", 3 * 1:4)), dat)),
    expectation[!is.na(expectation)]
  )
})

test_that("QuartetConcordance() handles incomplete data", {
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      rep("?", 8)), 8,
                    dimnames = list(paste0("t", 1:8), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  expect_equal(unname(QuartetConcordance(tree, dat)), rep(NA_real_, 5))
})

test_that("QuartetConcordance() handles non-integer data", {
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  mataset <- matrix(c("A", "A", "[AC]", "C", "C", "C", "T", "T", rep("?", 8)),
                    8, dimnames = list(paste0("t", 1:8), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  intSet <- matrix(c("1", "1", "[12]", "2", "2", "2", "4", "4", rep("?", 8)),
                    8, dimnames = list(paste0("t", 1:8), NULL))
  
  expect_equal(QuartetConcordance(tree, dat),
               QuartetConcordance(tree, MatrixToPhyDat(intSet)))
})

test_that("QuartetConcordance() unit = 'nrqs' locks the contract", {
  tree <- BalancedTree(8)

  # Character identical to the {t1..t4 | t5..t8} split: that split scores 1,
  # and it is the *only* split scoring 1 (nested/crossing get partial credit).
  identChar <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 1, 1, 1, 1), 8, dimnames = list(paste0("t", 1:8), NULL)))
  qt <- QuartetConcordance(tree, identChar, unit = "nrqs",
                           chanceCorrect = FALSE)
  expect_equal(max(qt, na.rm = TRUE), 1)
  expect_equal(sum(abs(qt - 1) < 1e-9, na.rm = TRUE), 1L)
  expect_true(all(qt >= 0 & qt <= 1, na.rm = TRUE))

  # `unit = "nrqs"` and `chanceCorrect = TRUE` are the defaults, so a bare call
  # is byte-identical to naming them.
  expect_identical(QuartetConcordance(tree, identChar),
                   QuartetConcordance(tree, identChar, unit = "nrqs",
                                      chanceCorrect = TRUE))
  # Coverage normalisation makes NRQS no laxer than quartet.
  expect_true(mean(qt, na.rm = TRUE) <=
                mean(QuartetConcordance(tree, identChar, unit = "quartet",
                                        chanceCorrect = FALSE),
                     na.rm = TRUE))

  # Uninformative characters (constant / autapomorphy) carry no NRQS -> NA.
  autap <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 0, 0, 0, 1), 8, dimnames = list(paste0("t", 1:8), NULL)))
  expect_equal(unname(QuartetConcordance(tree, autap, unit = "nrqs")),
               rep(NA_real_, 5))

  # `unit` is validated.
  expect_error(QuartetConcordance(tree, identChar, unit = "trits"),
               "should be one of")

  # `return` aliases mirror the quartet path.
  ByReturn <- function(r) {
    QuartetConcordance(tree, identChar, return = r, unit = "nrqs")
  }
  expect_equal(ByReturn("edge"), ByReturn("default"))
  cA <- ByReturn("char")
  expect_equal(cA, ByReturn("character"))
  expect_equal(cA, ByReturn("site"))
})

test_that("QuartetConcordance() unit = 'nrqs' gives nested partial credit", {
  # {t1..t5 | t6,t7,t8} split; state {t1,t2,t4} nested within side A but not a
  # clade, giving cells (p,q,r,s) = (3,0,2,3) -> A/Wk = 4/8 = 0.5 exactly.
  tree <- ape::read.tree(text = "(((((t1,t2),t3),t4),t5),(t6,(t7,t8)));")
  char <- MatrixToPhyDat(matrix(
    c(0, 0, 1, 0, 1, 1, 1, 1), 8, dimnames = list(paste0("t", 1:8), NULL)))
  qt <- QuartetConcordance(tree, char, unit = "nrqs", chanceCorrect = FALSE)

  # The split isolating {t6,t7,t8} scores exactly the nested value 0.5.
  sp <- as.Splits(tree)
  tips <- TipLabels(sp)
  member <- vapply(seq_along(sp), function(i) {
    side <- tips[as.logical(sp[[i]])]
    setequal(side, c("t6", "t7", "t8")) ||
      setequal(setdiff(tips, side), c("t6", "t7", "t8"))
  }, logical(1))
  expect_equal(unname(qt[member]), 0.5)
})

test_that("QuartetConcordance() unit = 'nrqs' supports multistate", {
  # Multistate no longer errors: NRQS sum over state-pairs (same currency).
  tree <- BalancedTree(8)
  ms <- MatrixToPhyDat(matrix(
    c(0, 0, 1, 1, 2, 2, 2, 0,      # 3 states
      0, 0, 0, 1, 1, 2, 3, 3), 8,  # 4 states
    dimnames = list(paste0("t", 1:8), NULL)))
  qt <- QuartetConcordance(tree, ms, unit = "nrqs", chanceCorrect = FALSE)
  expect_length(qt, 5L)
  expect_true(all(qt >= 0 & qt <= 1, na.rm = TRUE))
  expect_false(anyNA(qt))  # both characters are informative on this tree
})

test_that("QuartetConcordance() unit = 'nrqs' scores 1 iff split displayed", {
  # A multistate character need not be *identical* to a split to score 1: it
  # scores full marks for every split its own tree displays (each state block
  # wholly on one side; the split-orthogonal state-pair drops via M = 0).  This
  # generalises "only identical scores 1" from binary to multistate -- a cleanly
  # refining reproductive character fully supports the clade it refines.
  tree <- ape::read.tree(text = "((((t1,t2),t3),((t4,t5),t6)),((t7,t8),t9));")
  char <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 1, 1, 1, 2, 2, 2), 9,
    dimnames = list(paste0("t", 1:9), NULL)))
  qt <- QuartetConcordance(tree, char, unit = "nrqs", chanceCorrect = FALSE)
  sp <- as.Splits(tree)
  tips <- TipLabels(sp)
  atSplit <- function(members) {
    i <- which(vapply(seq_along(sp), function(j) {
      side <- tips[as.logical(sp[[j]])]
      setequal(side, members) || setequal(setdiff(tips, side), members)
    }, logical(1)))
    unname(qt[i])
  }
  # Splits the character displays -> exactly 1
  expect_equal(atSplit(c("t1", "t2", "t3")), 1)
  expect_equal(atSplit(c("t4", "t5", "t6")), 1)
  expect_equal(atSplit(c("t1", "t2", "t3", "t4", "t5", "t6")), 1)
  # A split cutting through a state block is only partially supported
  expect_true(atSplit(c("t1", "t2")) < 1)
})

test_that("QuartetConcordance() 'nrqs' return = 'char' has no ceiling of 1", {
  # Documented limitation, not an oversight: the character score is the
  # M-weighted mean across EVERY split, and a character agrees exactly with at
  # most one of them while being merely compatible with -- silent about -- the
  # rest.  Coverage scores silence as failure to cover, so even a character
  # identical to a split of the tree falls short of 1.  `unit = "quartet"` is
  # the per-character measure that does attain 1.
  tree <- BalancedTree(8)
  tips <- TipLabels(tree)
  sp <- as.Splits(tree)
  ByUnit <- function(u, r) {
    QuartetConcordance(tree, char, return = r, unit = u, chanceCorrect = FALSE)
  }
  for (k in seq_along(sp)) {
    char <- MatrixToPhyDat(matrix(
      as.integer(as.logical(sp[[k]])), length(tips),
      dimnames = list(tips, NULL)))

    # Per-split bound holds: the character scores 1 at the split it duplicates
    expect_equal(max(ByUnit("nrqs", "edge")), 1)

    # Naive currency reaches 1 on the character path ...
    expect_equal(ByUnit("quartet", "char")[[1]], 1)

    # ... the NRQS currency does not, and must not be read against 1
    nrqsChar <- ByUnit("nrqs", "char")[[1]]
    expect_gt(nrqsChar, 0)
    expect_lt(nrqsChar, 1)
  }
})

test_that("QuartetConcordance() unit = 'nrqs' chanceCorrect contract", {
  tree <- ape::read.tree(text = "((((t1,t2),t3),(t4,(t5,t6))),(t7,(t8,t9)));")
  ident <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 0, 0, 1, 1, 1), 9,
    dimnames = list(paste0("t", 1:9), NULL)))

  # chanceCorrect = TRUE is the default; re-zeroing can only lower the measure,
  # since a non-negative baseline is subtracted.
  expect_identical(
    QuartetConcordance(tree, ident, unit = "nrqs", chanceCorrect = TRUE),
    QuartetConcordance(tree, ident, unit = "nrqs"))
  expect_true(all(
    QuartetConcordance(tree, ident, unit = "nrqs", chanceCorrect = TRUE) <=
      QuartetConcordance(tree, ident, unit = "nrqs", chanceCorrect = FALSE),
    na.rm = TRUE))

  # Invalid chanceCorrect is rejected.
  for (bad in list(0, -3, "x")) {
    expect_error(QuartetConcordance(tree, ident, chanceCorrect = bad),
                 "positive integer")
  }

  # Chance correction also works for the raw quartet currency.
  qc <- QuartetConcordance(tree, ident, unit = "quartet", chanceCorrect = TRUE)
  expect_type(qc, "double")
  expect_false(anyNA(qc))

  # An identical (displayed) split still scores exactly 1 after re-zeroing:
  # .Rezero(1, z) == 1, so chance correction lifts the floor without moving the
  # ceiling.
  sp <- as.Splits(tree)
  tips <- TipLabels(sp)
  col789 <- which(vapply(seq_along(sp), function(j) {
    side <- tips[as.logical(sp[[j]])]
    setequal(side, c("t7", "t8", "t9")) ||
      setequal(setdiff(tips, side), c("t7", "t8", "t9"))
  }, logical(1)))
  corrected <- QuartetConcordance(tree, ident, unit = "nrqs",
                                  chanceCorrect = TRUE)
  expect_equal(unname(corrected[col789]), 1)
})

test_that("QuartetConcordance() unit = 'nrqs' exact baseline matches MC", {
  tree <- ape::read.tree(text = "((((t1,t2),t3),(t4,(t5,t6))),(t7,(t8,t9)));")
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 0, 0, 1, 1, 1,      # displays {t7,t8,t9}
      0, 0, 0, 0, 0, 0, 0, 1, 1,      # nested: {t8,t9}
      0, 0, 1, 1, 0, 0, 1, 1, 0,      # crossing (binary)
      0, 0, 0, 1, 1, 1, 2, 2, 2), 9,  # 3-state: exercises the random-wk path
    dimnames = list(paste0("t", 1:9), NULL)))

  # The exact hypergeometric baseline and the Monte-Carlo tip-shuffle baseline
  # target the same null, so they agree to within MC error -- for edge and char,
  # under both weightings, and crucially through the MULTISTATE pair path where a
  # pair's wk and M = min(w_c, w_k) are themselves random under the null.
  for (ret in c("edge", "char")) {
    for (w in c(TRUE, FALSE)) {
      exact <- QuartetConcordance(tree, dat, unit = "nrqs", return = ret,
                                  weight = w, chanceCorrect = TRUE)
      set.seed(1)
      mc <- QuartetConcordance(tree, dat, unit = "nrqs", return = ret,
                               weight = w, chanceCorrect = 8000)
      # Same cells resolve to NA under exact and MC (guards the cell-matching in
      # the weight = FALSE path), and the finite values agree within MC error.
      expect_identical(is.na(exact), is.na(mc))
      expect_lt(max(abs(exact - mc), na.rm = TRUE), 0.04)
    }
  }
})

test_that("QuartetConcordance() unit = 'quartet' exact baseline matches MC", {
  tree <- ape::read.tree(text = "((((t1,t2),t3),(t4,(t5,t6))),(t7,(t8,t9)));")
  dat <- MatrixToPhyDat(matrix(
    c(0, 0, 0, 0, 0, 0, 1, 1, 1,      # displays {t7,t8,t9}
      0, 0, 0, 0, 0, 0, 0, 1, 1,      # nested
      0, 0, 1, 1, 0, 0, 1, 1, 0,      # crossing (binary)
      0, 0, 0, 1, 1, 1, 2, 2, 2), 9,  # 3-state: multistate decisive path
    dimnames = list(paste0("t", 1:9), NULL)))

  # The raw-currency chance baseline (exact E[conc]/E[dec] from the
  # hypergeometric pmf) matches the Monte-Carlo tip-shuffle baseline within MC
  # error, across edge/char and both weightings, incl. the multistate path.
  for (ret in c("edge", "char")) {
    for (w in c(TRUE, FALSE)) {
      exact <- QuartetConcordance(tree, dat, unit = "quartet", return = ret,
                                  weight = w, chanceCorrect = TRUE)
      set.seed(1)
      mc <- QuartetConcordance(tree, dat, unit = "quartet", return = ret,
                               weight = w, chanceCorrect = 8000)
      expect_identical(is.na(exact), is.na(mc))
      expect_lt(max(abs(exact - mc), na.rm = TRUE), 0.04)
    }
  }

  # The published raw measure is recovered by naming both non-default settings,
  # and the default no longer returns it.
  for (ret in c("edge", "char")) for (w in c(TRUE, FALSE)) {
    raw <- QuartetConcordance(tree, dat, unit = "quartet", return = ret,
                              weight = w, chanceCorrect = FALSE)
    expect_true(all(raw >= 0 & raw <= 1, na.rm = TRUE))
    expect_false(identical(
      raw,
      QuartetConcordance(tree, dat, return = ret, weight = w)))
  }
})

test_that("QuartetConcordance() 'nrqs' correction: conflict below random", {
  tree <- ape::read.tree(text = "((((t1,t2),t3),(t4,(t5,t6))),(t7,(t8,t9)));")
  support <- MatrixToPhyDat(matrix(   # agrees with the tree: displays {t7,t8,t9}
    c(0, 0, 0, 0, 0, 0, 1, 1, 1), 9,
    dimnames = list(paste0("t", 1:9), NULL)))
  conflict <- MatrixToPhyDat(matrix(  # crosses the two halves of the tree
    c(0, 0, 1, 1, 0, 0, 1, 1, 0), 9,
    dimnames = list(paste0("t", 1:9), NULL)))

  sChar <- QuartetConcordance(tree, support, return = "char", unit = "nrqs",
                              chanceCorrect = TRUE)
  cNone <- QuartetConcordance(tree, conflict, return = "char", unit = "nrqs",
                              chanceCorrect = FALSE)
  cChar <- QuartetConcordance(tree, conflict, return = "char", unit = "nrqs",
                              chanceCorrect = TRUE)
  # A character that agrees with the tree scores above random expectation;
  # a crossing character is pulled below it, and below its uncorrected score.
  expect_gt(sChar, 0)
  expect_lt(cChar, cNone)
  expect_lt(cChar, 0)
})

test_that(".Rezero() works", {
  expect_equal(TreeSearch:::.Rezero(seq(0, 1, by = 0.1), 0.1), -1:9 / 9)
})

test_that("ConcordanceTable() marginSize top/right strips", {
  skip_if_not_installed("vdiffr")
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]][, 1:20]
  tree <- TreeSearch::referenceTree

  vdiffr::expect_doppelganger("conc-tbl-xx34", function() {
    expect_named(
      ConcordanceTable(tree, dataset, marginSize = c(NA, NA, 2, 2)),
      c("info", "relInfo", "quality", "col")
    )
  })

  vdiffr::expect_doppelganger("conc-tbl-all", function() {
    ConcordanceTable(tree, dataset, marginSize = c(2, 2, 2, 2))
  })
  vdiffr::expect_doppelganger("conc-tbl-2", function() {
    capture.output(ConcordanceTable(tree, dataset, marginSize = 2))
  })
})

test_that("ConcordanceTable() paintSize strips", {
  skip_if_not_installed("vdiffr")
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]][, 1:20]
  tree <- TreeSearch::referenceTree

  vdiffr::expect_doppelganger("conc-tbl-paint-scalar", function() {
    ConcordanceTable(tree, dataset, paintSize = 1)
  })
  vdiffr::expect_doppelganger("conc-tbl-paint-with-margin", function() {
    ConcordanceTable(tree, dataset, marginSize = 2, paintSize = 1)
  })
})

test_that("ClusteringConcordance() gives sensible values", {
  tree <- BalancedTree(8)
  splits <- as.Splits(tree)
  # None of these characters are informative
  mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      rep("?", 8)), 8,
                    dimnames = list(paste0("t", 1:8), NULL))
  dat <- MatrixToPhyDat(mataset)
  
  expect_equal(unname(ClusteringConcordance(tree, dat)), rep(NA_real_, 5))
  
  tree <- ape::read.tree(text = "((a, b, c, d, e), (f, g, h));")
  split <- as.Splits(tree)
  
  mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 1, 1, 1, # Matches split
                      0, 0, 0, 0, 1, 1, 1, 1, # Consistent but not identical
                      0, 0, 0, 1, 1, 1, 1, 1, # Consistent, more different
                      0, 0, 0, 0, 0, 0, 1, 1, # Consistent other way
                      0, 1, 0, 1, 0, 1, 0, 1, # Worst possible
                      0, 0, 0, 0, rep("?", 4), # No information
                      0, 0, 1, 1, rep("?", 4), # No relevant information
                      rep("?", 8)), 8,
                    dimnames = list(letters[1:8], NULL))
  dat <- MatrixToPhyDat(mataset)
  cc <- ClusteringConcordance(tree, dat, return = "all")[, "10", ]
  .Entropy <- function(...) {
    TreeDist::Entropy(c(...) / sum(...))
  }
  .NormExp <- function(a, b, ab) {
    TreeSearch:::.Rezero(
      (.Entropy(a) + .Entropy(b) - .Entropy(ab)) / .Entropy(a),
      TreeSearch:::.ExpectedMI(a, b) / .Entropy(a)
    )
  }
  expect_equal(cc["normalized", ],
               c(NA_real_, 1, 
                 .NormExp(c(3, 5), c(4, 4), c(1, 3, 4)),
                 .NormExp(c(3, 5), c(3, 5), c(2, 3, 3)),
                 .NormExp(c(2, 6), c(3, 5), c(2, 1, 5)),
                 .NormExp(c(3, 5), c(4, 4), c(2, 1, 2, 3)), 
                 NA, NA, NA))
  
  randomset <- matrix(sample(0:1, 8 * 1000, replace = TRUE), 8,
                      dimnames = list(letters[1:8], NULL))
  rat <- MatrixToPhyDat(randomset)
  expect_equal(ClusteringConcordance(tree, rat, chanceCorrect = TRUE),
               c("10" = 0),
               tolerance = 0.05)
})

test_that("ClusteringConcordance(return = 'char') Monte-Carlo handles ambiguity", {
  # Regression for T-330: characters whose ambiguous tokens drop different tips
  # give `charSplits` over heterogeneous tip sets. The Monte-Carlo
  # `chanceCorrect`
  # path scored a *list* of random trees against that list in one call, which
  # could not reconcile a common label set ("Old and new labels must match").
  tree <- ape::read.tree(text = "((a, b, c, d, e), (f, g, h));")
  mataset <- matrix(c(0, 0, 0, 0, 0, 0, 0, 1,
                      0, 0, 0, 0, 0, 1, 1, 1,
                      0, 0, 0, 0, rep("?", 4), # drops 4 tips
                      0, 0, 1, 1, rep("?", 4), # drops a different 4 tips
                      rep("?", 8)),            # all ambiguous
                    8, dimnames = list(letters[1:8], NULL))
  dat <- MatrixToPhyDat(mataset)

  set.seed(1)
  # Previously errored: "Old and new labels must match"
  cc <- ClusteringConcordance(tree, dat, return = "char", chanceCorrect = 10L)
  nChar <- length(attr(dat, "index"))
  expect_length(cc, nChar)
  expect_length(attr(cc, "mcse"), nChar)
  # Character 2 matches the only split perfectly, so normalization leaves it at 1
  expect_equal(unname(cc[2]), 1)
  expect_equal(unname(attr(cc, "mcse")[2]), 0)
  # Monte-Carlo score never exceeds the un-normalized score (subtracts a baseline)
  bare <- ClusteringConcordance(tree, dat, return = "char",
                                chanceCorrect = FALSE)
  finite <- is.finite(cc) & is.finite(bare)
  expect_true(all(cc[finite] <= bare[finite] + 1e-8))
})

test_that("ConcordantInformation() works", {
  data(congreveLamsdellMatrices)
  dat <- congreveLamsdellMatrices[[10]]
  tree <- TreeTools::NJTree(dat)
  
  expect_message(ci <- ConcordantInformation(tree, dat),
                 "dataset contains .* bits")
  expect_warning(suppressMessages(eval_val <- Evaluate(tree, dat)))
  expect_equal(eval_val, ci)
  expect_equal(TreeLength(tree, dat, concavity = "prof"),
               unname(ci["noise"]))
  expect_equal(Log2Unrooted(22), unname(ci["treeInformation"]))
  expect_equal(sum(apply(PhyDatToMatrix(dat), 2, CharacterInformation)),
               unname(ci["informationContent"]))
  
  dataset <- MatrixToPhyDat(cbind(setNames(c(rep(1, 11), 2:5), paste0("t", 1:15))))
  tree <- TreeTools::PectinateTree(length(dataset))
  # All non-1 states are singletons → no informative characters → zero concordance
  ci_empty <- suppressMessages(ConcordantInformation(tree, dataset))
  expect_equal(0, unname(ci_empty["signal"]))
  expect_equal(0, unname(ci_empty["noise"]))
  
  dataset <- MatrixToPhyDat(c(a = 1, b = 2, c = 1, d = 2, e = 3, f = 3))
  tree <- TreeTools::PectinateTree(dataset)
  # After T-107, 3-state chars with 6 tips are within the MaddisonSlatkin
  # feasibility threshold (k=3, max=15 tips), so no binary reduction occurs.
  # Signal/noise are computed via the full 3-state profile (no warning).
  ci <- suppressMessages(ConcordantInformation(tree, dataset))
  expect_equal(c(signal = 0.7835082), ci["signal"], tolerance = 1e-5)
  expect_equal(c(noise = 3.1233824), ci["noise"], tolerance = 1e-5)
  expect_equal(c(ignored = 0), ci["ignored"])
  
})

test_that("QACol() handles input", {
  expect_equal(is.na(QACol(c(0, 1, NA, NA, 0), c(0, 1, NA, 0, NA))),
               c(FALSE, FALSE, TRUE,
                 FALSE, # No data = black, quality NA by definition
                 TRUE))
  expect_equal(is.na(QCol(c(0, 1, NA, NA), NA)), c(FALSE, FALSE, TRUE, TRUE))
})
