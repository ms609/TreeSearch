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

test_that("Bremer errors when optimalScore exceeds the reference tree length", {
  dat <- StringToPhyDat("1100000 1110000 1111000 1111100 1001000",
                        1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  # A claimed optimum longer than the reference tree's own length is impossible
  # under these scoring arguments; the scoring guard catches it before searching
  # (a stale or mis-scaled L*) rather than reporting a bogus decay.
  expect_error(
    Bremer(ref, dat, method = "constraint",
           optimalScore = attr(mpts, "score") + 5,
           maxReplicates = 20L, verbosity = 0L),
    "exceeds an achievable length")
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
  # subtract an IW optimum from EW lengths.  Exact detection is infeasible (an
  # arbitrary resolution only upper-bounds L*, and computing L* is NP-hard), so
  # the guard WARNS on the large gap rather than erroring or silently mis-scoring
  # (BR-1).
  expect_warning(Bremer(mpts, dat, maxReplicates = 8L, verbosity = 0L),
                 "differs substantially")
  # Passing the matching concavity is accepted without the mismatch warning.
  ok <- suppressWarnings(
    Bremer(mpts, dat, concavity = 3, maxReplicates = 12L, verbosity = 0L))
  expect_type(ok, "double")
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
