library("TreeTools", quietly = TRUE)

# Exact Bremer by exhaustive enumeration of every binary tree (<= 7 tips) -- the
# same gold-standard oracle used by test-Bremer.R.  The parallel fan-out is
# validated to hit exactly this, seed- and schedule-independently.
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
  trees <- structure(lapply(allBinaryTrees(tips), RootTree, tips[1]),
                     class = "multiPhylo")
  lengths <- TreeLength(trees, dataset, concavity = concavity)
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

# A per-clade task mirroring .BremerConstraint()'s processConverse: an EW,
# serial converse search (worse-accepting phases disabled) returning the raw
# optimal score of the tree forced to lack clade `i`, with the engine's
# "no tree found" sentinel mapped to NA.  `TreeSearch::` qualification so the
# closure resolves the same way inside a cluster worker.
makeProcessFn <- function(dataset, splitList) {
  function(i) {
    s <- attr(TreeSearch::MaximizeParsimony(
      dataset, collapse = TRUE, nThreads = 1L,
      .negativeConstraint = splitList[[i]],
      driftCycles = 0L, sectorGoDrift = 0L, sectorDriftCycles = 0L,
      annealCycles = 0L,
      maxReplicates = 25L, verbosity = 0L), "score")
    if (!is.finite(s) || s < 0) NA_real_ else s
  }
}

# Shared fixture: a graded pectinate 7-taxon matrix with clades of varying
# support, so the fan-out is exercised over several distinct decay values.
bremerFixture <- function() {
  dat <- StringToPhyDat(
    "1100000 1110000 1111000 1111100 1100000 1110000 1111000 1111100 1001000",
    1:7, byTaxon = FALSE)
  names(dat) <- c(LETTERS[1:6], "out")
  set.seed(1)
  mpts <- MaximizeParsimony(dat, maxReplicates = 8L, verbosity = 0L)
  ref <- mpts[[1]]
  splits <- as.Splits(ref, tipLabels = names(dat))
  splitNames <- rownames(as.matrix(splits))
  list(dat = dat, ref = ref, Lstar = attr(mpts, "score"),
       splits = splits, splitNames = splitNames,
       splitList = lapply(seq_along(splitNames), function(i) splits[[i]]))
}

test_that(".BremerConverseScores (serial fan-out) matches the enumeration oracle", {
  skip_on_cran()
  fx <- bremerFixture()
  n <- length(fx$splitList)

  set.seed(7)
  raw <- TreeSearch:::.BremerConverseScores(
    n, makeProcessFn(fx$dat, fx$splitList), cl = NULL)
  bremer <- setNames(raw - fx$Lstar, fx$splitNames)

  orc <- oracleBremer(fx$ref, fx$dat)
  expect_length(bremer, length(orc))
  expect_equal(as.numeric(bremer[names(orc)]), as.numeric(orc), tolerance = 1e-6)
})

test_that(".BremerConverseScores is reproducible under set.seed (per-clade seeding)", {
  skip_on_cran()
  fx <- bremerFixture()
  n <- length(fx$splitList)
  pf <- makeProcessFn(fx$dat, fx$splitList)

  set.seed(7)
  a <- TreeSearch:::.BremerConverseScores(n, pf, cl = NULL)
  set.seed(7)
  b <- TreeSearch:::.BremerConverseScores(n, pf, cl = NULL)
  # Serial MaximizeParsimony is reproducible under set.seed, and the fan-out
  # draws the per-clade seeds in a fixed order, so two runs are bit-identical.
  expect_identical(a, b)

  # Explicit per-clade seeds bypass the RNG draw entirely.
  seeds <- c(11L, 22L, 33L, 44L, 55L)[seq_len(n)]
  c1 <- TreeSearch:::.BremerConverseScores(n, pf, cl = NULL, seeds = seeds)
  c2 <- TreeSearch:::.BremerConverseScores(n, pf, cl = NULL, seeds = seeds)
  expect_identical(c1, c2)
})

test_that(".BremerConverseScores handles an empty clade list", {
  expect_identical(
    TreeSearch:::.BremerConverseScores(0L, function(i) 0, cl = NULL),
    numeric(0))
})

test_that(".BremerConverseScores rejects a non-cluster cl", {
  expect_error(
    TreeSearch:::.BremerConverseScores(3L, function(i) 0, cl = 4L),
    "must be a `parallel` cluster")
})

test_that(".BremerConverseScores (PSOCK cluster) matches serial and the oracle", {
  skip_on_cran()
  skip_on_ci()  # PSOCK worker startup + load_all is slow/flaky in CI
  if (!requireNamespace("parallel", quietly = TRUE) ||
      !requireNamespace("pkgload", quietly = TRUE)) {
    skip("parallel / pkgload unavailable")
  }
  pkgPath <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
  if (is.null(pkgPath)) skip("cannot locate package source for worker load_all")

  cl <- tryCatch(parallel::makePSOCKcluster(2L), error = function(e) NULL)
  if (is.null(cl)) skip("cannot create PSOCK cluster")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  loaded <- tryCatch({
    parallel::clusterCall(cl, function(p) {
      suppressMessages(pkgload::load_all(p, quiet = TRUE, helpers = FALSE,
                                         attach_testthat = FALSE))
      TRUE
    }, pkgPath)
  }, error = function(e) NULL)
  if (is.null(loaded) || !all(vapply(loaded, isTRUE, logical(1)))) {
    skip("cannot load TreeSearch on cluster workers")
  }

  fx <- bremerFixture()
  n <- length(fx$splitList)
  pf <- makeProcessFn(fx$dat, fx$splitList)
  seeds <- seq_len(n) + 100L  # fixed seeds -> deterministic regardless of worker

  par <- TreeSearch:::.BremerConverseScores(n, pf, cl = cl, seeds = seeds)
  ser <- TreeSearch:::.BremerConverseScores(n, pf, cl = NULL, seeds = seeds)
  # Same per-clade seeds -> identical searches whether run in parallel or serial.
  expect_equal(par, ser, tolerance = 1e-6)

  bremer <- setNames(par - fx$Lstar, fx$splitNames)
  orc <- oracleBremer(fx$ref, fx$dat)
  expect_equal(as.numeric(bremer[names(orc)]), as.numeric(orc), tolerance = 1e-6)
})

test_that("public Bremer(cl=) matches the serial path on saturating data", {
  # Exercises the WHOLE public dispatch (Bremer -> .BremerConstraint -> fan-out),
  # which the .BremerConverseScores tests above do not.  On 7 tips the converse
  # searches saturate, so the serial and parallel paths -- despite their
  # different per-clade RNG streams -- reach the same optima and must agree.
  skip_on_cran()
  skip_on_ci()  # PSOCK worker startup + load_all is slow/flaky in CI
  if (!requireNamespace("parallel", quietly = TRUE) ||
      !requireNamespace("pkgload", quietly = TRUE)) {
    skip("parallel / pkgload unavailable")
  }
  pkgPath <- tryCatch(pkgload::pkg_path(), error = function(e) NULL)
  if (is.null(pkgPath)) skip("cannot locate package source for worker load_all")
  cl <- tryCatch(parallel::makePSOCKcluster(2L), error = function(e) NULL)
  if (is.null(cl)) skip("cannot create PSOCK cluster")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  loaded <- tryCatch(
    parallel::clusterCall(cl, function(p) {
      suppressMessages(pkgload::load_all(p, quiet = TRUE, helpers = FALSE,
                                         attach_testthat = FALSE))
      TRUE
    }, pkgPath), error = function(e) NULL)
  if (is.null(loaded) || !all(vapply(loaded, isTRUE, logical(1)))) {
    skip("cannot load TreeSearch on cluster workers")
  }

  fx <- bremerFixture()
  set.seed(1)
  mpts <- MaximizeParsimony(fx$dat, maxReplicates = 8L, verbosity = 0L)
  set.seed(5)
  ser <- Bremer(mpts, fx$dat, maxReplicates = 20L, verbosity = 0L)
  set.seed(5)
  par <- Bremer(mpts, fx$dat, cl = cl, maxReplicates = 20L, verbosity = 0L)
  expect_equal(unname(par), unname(ser), tolerance = 1e-6)
  expect_identical(names(par), names(ser))
})
