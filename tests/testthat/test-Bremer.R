library("TreeTools", quietly = TRUE)

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

# Exact Bremer by exhaustive enumeration of every binary tree on the taxa
# (feasible for <= 7 tips): the minimum length among trees LACKING each
# reference clade, minus the global minimum length.  The gold-standard oracle
# for the converse-constraint engine.
allBinaryTrees <- function(tips) {
  trees <- list(ape::read.tree(
    text = paste0("(", tips[1], ",", tips[2], ",", tips[3], ");")))
  for (k in 4:length(tips)) {
    trees <- unlist(lapply(trees, function(tr) {
      AddTipEverywhere(tr, tips[k], includeRoot = FALSE)
    }), recursive = FALSE)
  }
  trees
}

oracleBremer <- function(reference, dataset, concavity = Inf) {
  tips <- names(dataset)
  trees <- allBinaryTrees(tips)
  trees <- structure(lapply(trees, RootTree, tips[1]), class = "multiPhylo")
  lengths <- TreeLength(trees, dataset, concavity = concavity)
  Lstar <- min(lengths)
  refFreq <- SplitFrequency(reference, trees)
  nSplits <- length(refFreq)
  disp <- vapply(seq_along(trees), function(j) {
    SplitFrequency(reference, trees[j])
  }, double(nSplits))
  dim(disp) <- c(nSplits, length(trees))
  out <- vapply(seq_len(nSplits), function(i) {
    lacking <- disp[i, ] < 0.5
    if (any(lacking)) min(lengths[lacking]) - Lstar else Inf
  }, double(1))
  setNames(out, names(refFreq))
}

# Enumeration oracle scoring under the inapplicable "bgs" kernel (the DEFAULT
# scoring mode).  Validates the converse-constraint engine on data with
# inapplicable (`-`) tokens -- the path that runs exact_verify_sweep.
oracleBremerNA <- function(reference, dataset) {
  tips <- names(dataset)
  trees <- structure(lapply(allBinaryTrees(tips), RootTree, tips[1]),
                     class = "multiPhylo")
  lengths <- TreeLength(trees, dataset, inapplicable = "bgs")
  Lstar <- min(lengths)
  refFreq <- SplitFrequency(reference, trees)
  nSplits <- length(refFreq)
  disp <- vapply(seq_along(trees), function(j) {
    SplitFrequency(reference, trees[j])
  }, double(nSplits))
  dim(disp) <- c(nSplits, length(trees))
  setNames(vapply(seq_len(nSplits), function(i) {
    lacking <- disp[i, ] < 0.5
    if (any(lacking)) min(lengths[lacking]) - Lstar else Inf
  }, double(1)), names(refFreq))
}

# --- method = "pool" (approximate engine) ---
# A multiPhylo of several MPTs is summarised by its strict consensus; support
# is calculated for the consensus's resolved clades only.

test_that("Bremer(method = 'pool') returns node-keyed non-negative support", {
  set.seed(3418)
  trees <- MaximizeParsimony(ds, maxReplicates = 6L, targetHits = 2L,
                             verbosity = 0L)
  refCons <- if (length(trees) == 1L) trees[[1]] else ape::consensus(trees, p = 1)

  decay <- Bremer(trees, ds, method = "pool", maxBremer = 6,
                  maxReplicates = 6L, targetHits = 2L, verbosity = 0L)

  expect_type(decay, "double")
  # One value per resolved clade of the consensus reference.
  expect_length(decay, NSplits(refCons))
  # Names are internal node numbers.
  expect_true(all(as.integer(names(decay)) > NTip(ds)))

  # Finite (uncensored) support values are non-negative.
  expect_true(all(decay[is.finite(decay)] >= 0))

  # Censored attribute marks exactly the infinite entries.
  cens <- attr(decay, "censored")
  expect_type(cens, "logical")
  expect_length(cens, length(decay))
  expect_equal(unname(cens), unname(!is.finite(decay)))
})

test_that("Bremer(format = 'character') yields a node.label-shaped vector", {
  set.seed(3418)
  trees <- MaximizeParsimony(ds, maxReplicates = 6L, targetHits = 2L,
                             verbosity = 0L)
  refCons <- if (length(trees) == 1L) trees[[1]] else ape::consensus(trees, p = 1)

  lab <- Bremer(trees, ds, method = "pool", maxBremer = 6, format = "character",
                maxReplicates = 6L, targetHits = 2L, verbosity = 0L)
  expect_type(lab, "character")
  expect_length(lab, refCons[["Nnode"]])
  # Assignable straight onto the reference tree.
  refCons[["node.label"]] <- lab
  expect_length(refCons[["node.label"]], refCons[["Nnode"]])
})

test_that("Bremer accepts a single phylo reference with an explicit L*", {
  set.seed(3418)
  trees <- MaximizeParsimony(ds, maxReplicates = 6L, targetHits = 2L,
                             verbosity = 0L)
  Lstar <- attr(trees, "score")
  decay <- Bremer(trees[[1]], ds, method = "pool", maxBremer = 6,
                  optimalScore = Lstar,
                  maxReplicates = 6L, targetHits = 2L, verbosity = 0L)
  expect_type(decay, "double")
  expect_length(decay, NSplits(trees[[1]]))
  expect_true(all(decay[is.finite(decay)] >= 0))
  expect_setequal(names(decay), rownames(as.matrix(as.Splits(trees[[1]]))))
})

test_that("Bremer errors on a bad reference type", {
  expect_error(Bremer("not a tree", ds, method = "pool"),
               "must be a `phylo` or `multiPhylo`")
})

# --- method = "constraint" (rigorous engine): validated against brute force ---

test_that("the oracle enumerates every binary tree (guards the gate itself)", {
  # (2n - 5)!! unrooted binary trees; a gappy enumerator would make the oracle
  # itself an over-estimate that could rubber-stamp an over-estimating search.
  doubleFactorial <- function(m) prod(seq(m, 1, by = -2))
  expect_equal(length(allBinaryTrees(LETTERS[1:6])), doubleFactorial(2 * 6 - 5))
  expect_equal(length(allBinaryTrees(c(LETTERS[1:6], "out"))),
               doubleFactorial(2 * 7 - 5))
})

test_that("constraint Bremer matches enumeration when L* is computed internally", {
  # Exercises the full-strength L* search path (no optimalScore supplied).
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  con <- Bremer(ref, dat, method = "constraint",   # no optimalScore -> internal L*
                maxReplicates = 25L, verbosity = 0L)
  orc <- oracleBremer(ref, dat)
  expect_equal(as.numeric(con[names(orc)]), as.numeric(orc), tolerance = 1e-6)
})

test_that("constraint Bremer matches exact enumeration (graded pectinate)", {
  dat <- StringToPhyDat(
    "1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  con <- Bremer(ref, dat, method = "constraint",
                optimalScore = attr(mpts, "score"),
                maxReplicates = 25L, verbosity = 0L)
  orc <- oracleBremer(ref, dat)
  expect_equal(as.numeric(con[names(orc)]), as.numeric(orc), tolerance = 1e-6)
})

test_that("constraint Bremer matches enumeration in the adversarial case", {
  # {A,B,C} is supported by three characters (decay 3); the optimum DISPLAYS it
  # and the nearest tree lacking it is three steps longer -- the case a stalling
  # post-hoc filter would over-estimate.
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  con <- Bremer(ref, dat, method = "constraint",
                optimalScore = attr(mpts, "score"),
                maxReplicates = 25L, verbosity = 0L)
  orc <- oracleBremer(ref, dat)
  expect_equal(as.numeric(con[names(orc)]), as.numeric(orc), tolerance = 1e-6)
  expect_gte(max(con), 3)
})

test_that("constraint Bremer matches enumeration under implied weights", {
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, concavity = 3, maxReplicates = 8L,
                            verbosity = 0L)
  ref <- mpts[[1]]
  con <- Bremer(ref, dat, method = "constraint", concavity = 3,
                optimalScore = attr(mpts, "score"),
                maxReplicates = 25L, verbosity = 0L)
  orc <- oracleBremer(ref, dat, concavity = 3)
  expect_equal(as.numeric(con[names(orc)]), as.numeric(orc), tolerance = 1e-6)
  # Fractional decay under implied weights.
  expect_false(all(con == round(con)))
})

test_that("constraint is the default and never looser than an uncensored pool", {
  # The matrix leaves parts of the tree unresolved, so a fully-resolved reference
  # has decay-0 clades that OTHER equally-optimal trees lack.  A best-score pool
  # therefore breaks them, making the cross-engine comparison NON-vacuous (guards
  # the earlier bug where a consensus reference left every shared clade censored
  # and the key assertion compared empty vectors).
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(42)
  mptsFull <- MaximizeParsimony(dat, collapse = FALSE, maxReplicates = 12L,
                                verbosity = 0L)
  ref <- mptsFull[[1]]                       # a fully-resolved (binary) MPT
  Lstar <- attr(mptsFull, "score")

  con <- Bremer(ref, dat, optimalScore = Lstar, maxReplicates = 20L,
                verbosity = 0L)              # default method (constraint)
  expect_type(con, "double")
  expect_true(all(con >= 0))
  expect_true(all(is.finite(con)))

  pool <- Bremer(ref, dat, method = "pool", optimalScore = Lstar, maxBremer = 10,
                 maxReplicates = 50L, verbosity = 0L)
  shared <- intersect(names(con), names(pool))
  uncensored <- shared[is.finite(pool[shared])]
  expect_gt(length(uncensored), 0L)   # the pool broke at least one clade
  # Both are upper bounds on the true decay; the directed constraint search is
  # never looser than the incidental pool where the pool is uncensored.
  expect_true(all(con[uncensored] <= pool[uncensored] + 1e-6))
})

test_that("Bremer warns (does not error) on an inconsistent optimalScore, then proceeds", {
  dat <- StringToPhyDat("1100000 1110000 1111000 1111100 1001000",
                        1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  # A supplied L* inconsistent with the reference length under these scoring
  # arguments is WARNED, not errored: Bremer trusts the user's scoring choice and
  # proceeds (they may deliberately want a different analysis).  A value BELOW the
  # reference length exercises the warning without triggering the adopt-shorter-L*
  # path (which would add a second, unrelated warning).
  expect_warning(
    con <- Bremer(ref, dat, method = "constraint",
                  optimalScore = attr(mpts, "score") - 3,
                  maxReplicates = 20L, verbosity = 0L),
    "differs from the reference")
  expect_type(con, "double")
})

test_that("constraint Bremer matches enumeration on inapplicable (bgs) data", {
  # The default scoring mode (inapplicable = "bgs") previously had no ground-
  # truth oracle.  Validate the negative-constraint machinery + exact_verify_sweep
  # against exact enumeration on a 6-taxon matrix with inapplicable (`-`) tokens.
  dat <- StringToPhyDat("111000 111000 111000 110000 000110 -00-11 -00-11",
                        1:6, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:5], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 10L, verbosity = 0L)
  ref <- mpts[[1]]
  con <- Bremer(ref, dat, method = "constraint",
                optimalScore = attr(mpts, "score"),
                maxReplicates = 30L, verbosity = 0L)
  orc <- oracleBremerNA(ref, dat)
  expect_equal(as.numeric(con[names(orc)]), as.numeric(orc), tolerance = 1e-6)
})

test_that("Bremer warns on a likely scoring-units mismatch (IW trees vs EW default)", {
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, concavity = 3, maxReplicates = 8L, verbosity = 0L)
  # Trees found under implied weights; the default equal-weights Bremer would
  # subtract an IW optimum from EW lengths.  The recorded scoring signature makes
  # this an EXACT criterion mismatch (concavity 3 vs Inf), so it WARNS (accepting
  # and proceeding) rather than erroring or silently mis-scoring (BR-1/B-7).
  expect_warning(Bremer(mpts, dat, maxReplicates = 8L, verbosity = 0L),
                 "found under")
  # Passing the matching concavity is accepted without the mismatch warning.
  ok <- suppressWarnings(
    Bremer(mpts, dat, concavity = 3, maxReplicates = 12L, verbosity = 0L))
  expect_type(ok, "double")
})

test_that("MaximizeParsimony records a scoring signature Bremer trusts (B-7 provenance)", {
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)  # EW default
  sig <- attr(mpts, "scoring")
  expect_false(is.null(sig))
  expect_true(is.infinite(sig[["concavity"]]))       # equal weights
  expect_identical(sig[["inapplicable"]], "bgs")     # default kernel

  # Matching scoring arguments -> the exact provenance check passes -> no warning
  # (and, being all defaults here, this is the canonical Bremer() call).
  expect_warning(
    con <- Bremer(mpts, dat, maxReplicates = 12L, verbosity = 0L), NA)
  expect_true(all(is.finite(con)))

  # An IW analysis under matching concavity also matches its own signature.
  set.seed(1)
  mptsIW <- MaximizeParsimony(dat, concavity = 5, maxReplicates = 8L,
                              verbosity = 0L)
  expect_equal(attr(mptsIW, "scoring")[["concavity"]], 5)
  expect_warning(
    Bremer(mptsIW, dat, concavity = 5, maxReplicates = 8L, verbosity = 0L), NA)
})

test_that("the scoring-mismatch warning catches a small gap the 5% band missed (B-7)", {
  dat <- StringToPhyDat(
    "1110000 1110000 1110000 1100000 0001100 0001100 0000110 1111111",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  Lstar <- attr(mpts, "score")
  # A gap of 0.3 is < 5% of the reference length, so the old one-sided
  # 5%-of-length heuristic would NOT have warned (the P1 blind spot two criteria
  # placing L* close together); the tight two-sided tolerance now does.  Below
  # the reference length, so the adopt-shorter-L* path does not add a 2nd warning.
  expect_warning(
    Bremer(ref, dat, method = "constraint", optimalScore = Lstar - 0.3,
           maxReplicates = 15L, verbosity = 0L),
    "differs from the reference")
})

test_that("Bremer ignores (with a warning) drift/anneal args in a converse search", {
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  # driftCycles is a documented forwarded arg but unsound in a converse search;
  # it must be ignored (with a warning), not silently empty the pool (BR-2).
  expect_warning(
    con <- Bremer(mpts, dat, driftCycles = 5L, annealCycles = 3L,
                  maxReplicates = 12L, verbosity = 0L),
    "always disabled")
  expect_true(all(is.finite(con)))
})

# --- method = "constraint" on >= 12 tips: exercises the SECTORIAL search path ---
# rss_search / xss_search engage once a tree has >= 12 tips.  Their reduced-
# dataset sector solve is constraint-blind, so a reinsertion could rebuild a
# forbidden clade; the fix rejects such a reinsertion post-hoc (a reach guard --
# the pool backstop already guaranteed soundness).  This whole path is
# structurally invisible to the <= 7-tip enumeration oracle above (sectors need
# >= 12 tips), so it needs its own coverage.

# A graded pectinate matrix on `nIn` ingroup taxa + one outgroup: nested clades
# {t1..tk} (k = 2..nIn-1), each duplicated for graded support.
sectorFixture <- function(nIn = 13) {
  tips <- c(paste0("t", seq_len(nIn)), "out")
  nTip <- length(tips)
  mkChar <- function(k) paste0(c(rep("1", k), rep("0", nTip - k)), collapse = "")
  chars <- unlist(lapply(rep(2:(nIn - 1L), each = 2L), mkChar))
  dat <- StringToPhyDat(paste(chars, collapse = " "), seq_len(nTip),
                        byTaxon = FALSE)
  names(dat) <- tips
  dat
}

test_that("converse search on >= 12 tips (sectors engaged) returns clade-free trees", {
  skip_on_cran()
  dat <- sectorFixture(13)            # 14 tips -> sectorial search engages
  expect_gte(length(dat), 12L)
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 6L, verbosity = 0L)
  ref <- mpts[[1]]
  splits <- as.Splits(ref, tipLabels = names(dat))
  splitNames <- rownames(as.matrix(splits))
  expect_gt(length(splitNames), 0L)

  # Directly drive the negative-constraint sector path for a spread of clades and
  # confirm each returned tree genuinely LACKS its forbidden clade.  The pool
  # backstop guarantees this regardless of the sector fix, so it is a robust
  # soundness assertion; the fix's REACH benefit (fewer stranded replicates) is a
  # Hamilton-scale measurement, not asserted here.
  probe <- unique(round(seq(1, length(splitNames),
                            length.out = min(4L, length(splitNames)))))
  for (i in probe) {
    set.seed(100L + i)
    res <- MaximizeParsimony(dat, collapse = TRUE, nThreads = 1L,
                             .negativeConstraint = splits[[i]],
                             driftCycles = 0L, sectorGoDrift = 0L,
                             sectorDriftCycles = 0L, annealCycles = 0L,
                             maxReplicates = 12L, verbosity = 0L)
    tr <- if (inherits(res, "phylo")) res else res[[1L]]
    disp <- SplitFrequency(ref, structure(list(tr), class = "multiPhylo"))
    expect_lt(unname(disp[splitNames[i]]), 0.5)   # returned tree lacks the clade
  }
})

test_that("constraint Bremer on >= 12 tips is non-negative and <= pool", {
  skip_on_cran()
  dat <- sectorFixture(13)
  set.seed(2)
  mpts <- MaximizeParsimony(dat, maxReplicates = 6L, verbosity = 0L)
  ref <- mpts[[1]]
  Lstar <- attr(mpts, "score")

  con <- Bremer(ref, dat, method = "constraint", optimalScore = Lstar,
                maxReplicates = 10L, verbosity = 0L)
  expect_type(con, "double")
  expect_true(any(is.finite(con)))               # the sector path did not strand
  expect_true(all(con[is.finite(con)] >= 0))     # correct L* (no negative decay)

  pool <- Bremer(ref, dat, method = "pool", optimalScore = Lstar, maxBremer = 12,
                 maxReplicates = 30L, verbosity = 0L)
  shared <- intersect(names(con), names(pool))
  uncensored <- shared[is.finite(pool[shared]) & is.finite(con[shared])]
  # Both are upper bounds on the true decay; the directed constraint search is
  # never looser than the incidental pool where both are finite.
  expect_true(all(con[uncensored] <= pool[uncensored] + 1e-6))
})

# --- input validation / format consistency (red-team round 2) ---

test_that("Bremer rejects a non-finite or non-scalar optimalScore (B-8)", {
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  # NA_real_ slips past is.null() and previously crashed the one-sided scoring
  # guard with "missing value where TRUE/FALSE needed"; now a clear error.
  expect_error(
    Bremer(mpts, dat, optimalScore = NA_real_, maxReplicates = 6L,
           verbosity = 0L),
    "single finite number")
  expect_error(
    Bremer(mpts, dat, optimalScore = c(1, 2), maxReplicates = 6L,
           verbosity = 0L),
    "single finite number")
})

test_that("Bremer on a fully-unresolved reference warns and returns nothing (TA-9)", {
  dat <- StringToPhyDat("111000 110000 001100", 1:6, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:5], "out")
  star <- ape::stree(length(dat), type = "star", tip.label = names(dat))
  # No resolved internal clades -> nothing to calculate; early return, no search.
  expect_warning(
    res <- Bremer(star, dat, verbosity = 0L),
    "no resolved internal clades")
  expect_type(res, "double")
  expect_length(res, 0L)
})

test_that("Bremer(format='character') always carries a censored attribute (B-9)", {
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, collapse = FALSE, maxReplicates = 8L,
                            verbosity = 0L)
  ref <- mpts[[1]]
  Lstar <- attr(mpts, "score")
  # method = "constraint" never censors, so this is the "nothing censored" path
  # where the character format previously dropped the attribute entirely while
  # the numeric format kept an all-FALSE vector.
  lab <- Bremer(ref, dat, optimalScore = Lstar, format = "character",
                maxReplicates = 12L, verbosity = 0L)
  cens <- attr(lab, "censored")
  expect_false(is.null(cens))                 # attribute is present...
  expect_type(cens, "logical")
  expect_length(cens, ref[["Nnode"]])         # ...node.label-shaped...
  expect_false(any(cens))                     # ...and correctly all-FALSE.
  # The numeric format for the same call also carries a censored attribute.
  num <- Bremer(ref, dat, optimalScore = Lstar,
                maxReplicates = 12L, verbosity = 0L)
  expect_false(is.null(attr(num, "censored")))
})

# --- engine-failure handling, via a test seam (no real search needed) ---
# These deterministically drive the NA / clade-displaying paths that the real
# engine (C++ guard + pool backstop) is designed never to reach, so they cannot
# be provoked by an ordinary search.

test_that(".BremerProcessResult maps engine sentinels to NA (BR-6/TA-4)", {
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, collapse = FALSE, maxReplicates = 8L,
                            verbosity = 0L)
  ref <- mpts[[1]]
  sn <- rownames(as.matrix(as.Splits(ref, tipLabels = names(dat))))[[1]]

  # "No tree lacking the clade found" sentinel (negative / non-finite) -> NA.
  expect_identical(
    TreeSearch:::.BremerProcessResult(list(score = -1, tree = NULL), ref, sn),
    NA_real_)
  expect_identical(
    TreeSearch:::.BremerProcessResult(list(score = NA_real_, tree = NULL), ref, sn),
    NA_real_)

  # A returned tree that STILL displays the clade (the reference does) -> the
  # BR-6 defence-in-depth check warns and returns NA rather than a deflated decay.
  expect_warning(
    r <- TreeSearch:::.BremerProcessResult(list(score = 5, tree = ref), ref, sn),
    "still displays it")
  expect_identical(r, NA_real_)

  # A genuinely clade-free tree (a star displays no clade) -> score passes through.
  star <- ape::stree(length(dat), type = "star", tip.label = names(dat))
  expect_identical(
    TreeSearch:::.BremerProcessResult(list(score = 5, tree = star), ref, sn), 5)
})

test_that(".BremerConstraint emits the aggregate NA warning (VF-4/TA-3)", {
  dat <- StringToPhyDat("1111000 1111000 1100000", 1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, collapse = FALSE, maxReplicates = 8L,
                            verbosity = 0L)
  ref <- TreeSearch:::.BremerReference(mpts, dat,
                                       optimalScore = attr(mpts, "score"))
  expect_gte(length(ref$splitNames), 2L)
  sa <- list(concavity = Inf, extended_iw = TRUE, xpiwe_r = 0.5, xpiwe_max_f = 5,
             hierarchy = NULL, inapplicable = "bgs", hsj_alpha = 1.0)

  # Inject an engine stub returning the "no tree found" sentinel for every clade,
  # so every converse search NAs and the aggregate warning fires.  Lstar is
  # supplied via the reference, so no real search runs.
  stubNA <- function(negSplit) list(score = -1, tree = NULL)
  expect_warning(
    res <- TreeSearch:::.BremerConstraint(ref, dat, sa, Inf, list(),
                                          cl = NULL, .runConverse = stubNA),
    "returned NA")
  expect_length(res$bremer, length(ref$splitNames))
  expect_true(all(is.na(res$bremer)))
})
