library("TreeTools", quietly = TRUE)

test_that("Profile handles multi-state characters", {
  # 3-state char with 6 tips: now natively supported by MaddisonSlatkin
  # (feasible for k=3, n=6, threshold=15). No warning or binary reduction.
  dataset <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 3, f = 3))
  pd <- PrepareDataProfile(dataset)
  expect_equal(3L, length(attr(pd, "levels")))
  expect_s3_class(pd, "phyDat")
})

test_that("Constraints work", {
  # Morphy() emits informational cli messages via message()
  # ("Initialized N distinct constraints", "Modifying tree to match
  # constraint", ...) that are not under test here.  Wrap the whole
  # body in suppressMessages() to keep test output clean; warnings and
  # errors are still raised normally.  Verbose paths are exercised in
  # separate tests below.
  suppressMessages({
  constraint <- MatrixToPhyDat(c(a = 1, b = 1, c = 0, d = 0, e = 0, f = 0))
  characters <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  set.seed(0)
  # Morphy() defaults to verbosity = 3 and prints a banner / score lines to
  # stdout. The verbose output is not under test here, so silence it with
  # verbosity = 0L. (Verbose paths are exercised separately below.)
  ewResults <- Morphy(characters,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint,
                                 verbosity = 0L)
  expect_equal_tree(PectinateTree(letters[1:6]), ewResults[[1]])
  expect_equal(c(seed = 0, start = 1, final = 0),
               attr(ewResults, "firstHit"))
  expect_equal(names(ewResults), "start_1")
  expect_equal_tree(PectinateTree(letters[1:6]),
               Morphy(characters, concavity = "p",
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint,
                                 verbosity = 0L)[[1]])
  expect_equal_tree(PectinateTree(letters[1:6]),
               Morphy(characters, concavity = 10,
                                 PectinateTree(c("a", "b", "f", "d", "e", "c")),
                                 ratchIter = 0, constraint = constraint,
                                 verbosity = 0L)[[1]])
  # Start tree not consistent with constraint
  dataset <- characters
  tree <- PectinateTree(c("a", "c", "f", "d", "e", "b"))
  expect_equal_tree(PectinateTree(letters[1:6]),
               Morphy(characters,
                      PectinateTree(c("a", "c", "f", "d", "e", "b")),
                                 ratchIter = 0, constraint = constraint,
                                 verbosity = 0L)[[1]])


  dataset <- MatrixToPhyDat(matrix(c(0, 0, 1, 1, 1, 1, 1,
                                     1, 1, 1, 1, 0, 0, 0), ncol = 2,
                                   dimnames = list(letters[1:7], NULL)))
  constraint <- MatrixToPhyDat(matrix(c(0, 0, 1, "?", 1, 1,
                                        1, 1, 1,   1, 0, 0), ncol = 2,
                                      dimnames = list(letters[1:6], NULL)))
  # T-039 fixed: column-major indexing in build_constraint + Wagner guards
  cons <- consensus(Morphy(dataset, constraint = constraint,
                           verbosity = 0L),
                    rooted = TRUE)
  # Avoid %in%.Splits — S3 dispatch breaks in testthat's cloned namespace
  # (test_check / R CMD check). Compare bipartitions as plain logical vectors.
  split_in_splits <- function(sp, table) {
    tips <- attr(table, "tip.label")
    sp <- as.Splits(sp, tipLabels = tips)
    s <- as.logical(as.logical(sp))  # flatten to plain vector
    tab <- as.logical(table)
    if (!is.matrix(tab)) tab <- matrix(tab, nrow = 1)
    any(apply(tab, 1, function(r) all(s == r) || all(s == !r)))
  }
  expect_true(split_in_splits(
    as.Splits(as.logical(c(0, 0, 1, 1, 1)), letters[c(1:3, 5:6)]),
    as.Splits(DropTip(cons, c("d", "g")))))
  
  expect_true(split_in_splits(
    as.Splits(as.logical(c(0, 0, 0, 0, 1, 1)), letters[1:6]),
    as.Splits(DropTip(cons, "g"))))

  })  # end suppressMessages
})

test_that("Inconsistent constraints fail", {
  constraint <- MatrixToPhyDat(matrix(
    c(0, 1, 1, 1, 0, 0,
      1, 1, 1, 0, 0, 0), ncol = 2,
    dimnames = list(letters[1:6], NULL)))
  # Morphy() may emit cli messages before the error fires; suppress so
  # they do not leak into testthat output.
  expect_error(
    suppressMessages(
      Morphy(constraint,
             PectinateTree(c("a", "b", "f", "d", "e", "c")),
             ratchIter = 0, constraint = constraint,
             verbosity = 0L)
    )
  )
})

test_that("Morphy() times out", {
  # Do not run on CRAN: Writing R Extensions discourages testing timings
  skip_if(Sys.getenv("GITHUB_PAT") == "") # Run only on GH Actions
  
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[42]]
  startTime <- Sys.time()
  # Discard verbose progress output — the test is about wall-clock timing.
  invisible(capture.output(
    Morphy(dataset, ratchIter = 10000, tbrIter = 1, maxHits = 1,
           maxTime = 0)
  ))
  expect_gt(as.difftime(5, units = "secs"), Sys.time() - startTime)
})

test_that("Seed trees retained", {
  tree1 <- read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  tree2 <- read.tree(text = "(a, (b, (c, (f, (e, d)))));")
  badTree <- read.tree(text = "(f, (b, (c, (a, (e, d)))));")
  dat <- StringToPhyDat("110000 110000 111000 111000 111100 111001",
                        letters[1:6], byTaxon = FALSE)
  # verbosity = 4 deliberately exercises the most verbose printing path,
  # which mixes cli messages (stderr) and Rprintf output (stdout).
  # Capture both streams so nothing leaks into testthat output, and
  # assert that the verbose marker is present so a future change that
  # silently breaks the print path would fail this test.
  msg_lines <- character()
  stdout_lines <- capture.output(
    msg_lines <- capture.output(
      results <- Morphy(dataset = dat,
                        tree = c(tree1, tree2, badTree),
                        ratchIter = 0, verbosity = 4),
      type = "message"
    )
  )
  all_lines <- c(stdout_lines, msg_lines)
  expect_true(any(grepl("TREE SEARCH|Score|Initial score|Starting search",
                        all_lines)))
  expect_equal(attr(results, "firstHit"),
               c(seed = 2, start = 0, final = 0))
})

test_that("Mismatched tree/dataset handled with warnings", {
  treeAf <- read.tree(text = "(a, (b, (c, (d, (e, f)))));")
  treeBg <- read.tree(text = "(g, (b, (c, (d, (e, f)))));")
  datAf <- StringToPhyDat("110000 110000 111100 111000",
                              letters[1:6], byTaxon = FALSE)
  datAe <- StringToPhyDat("11000 11000 11110 11100",
                              letters[1:5], byTaxon = FALSE)
  datAg <- StringToPhyDat("1100000 1100000 1111000 1110000",
                              letters[1:7], byTaxon = FALSE)
  
  # QP emits cli messages ("Ignoring taxa...", "Initialized N distinct
  # constraints") in addition to the R-level warning the tests check
  # for.  Wrap each call in suppressMessages() to silence the
  # informational cli output while preserving warning capture.
  QP <- function (...) suppressMessages(
    Morphy(..., ratchIter = 0, maxHits = 1, verbosity = 0)
  )

  # Some calls emit multiple R warnings (one per mismatch type: tree-only,
  # dataset-only, constraint-only).  expect_warning() only captures the first;
  # additional warnings propagate as "unexpected" in edition 3.  Use a helper
  # that muffles every warning while asserting at least one was raised.
  check_warns <- function(expr) {
    warns <- 0L
    val <- withCallingHandlers(expr, warning = function(w) {
      warns <<- warns + 1L
      invokeRestart("muffleWarning")
    })
    expect_gt(warns, 0L)
    val
  }

  r1 <- check_warns(QP(datAf, treeBg));                        expect_equal(5, unname(NTip(r1)))
  r2 <- check_warns(QP(datAe, treeAf));                        expect_equal(5, unname(NTip(r2)))
  r3 <- check_warns(QP(datAg, treeAf));                        expect_equal(6, unname(NTip(r3)))
  r4 <- check_warns(QP(datAf, treeBg, constraint = datAe));    expect_equal(5, unname(NTip(r4)))
  expect_equal(6, unname(NTip(QP(datAf, treeAf, constraint = datAe))))
  r5 <- check_warns(QP(datAf, treeAf, constraint = datAg));    expect_equal(6, unname(NTip(r5)))
})

test_that("Root retained if not 1", {
  tr <- RootTree(BalancedTree(8), "t5")
  dataset <- StringToPhyDat("11000000 11100000 11110000 11111000",
                            paste0("t", 1:8), byTaxon = FALSE)
  
  mpt <- Morphy(dataset, tr, verbosity = 0L)
  expect_equal(5, mpt[[1]]$edge[14, 2])
})

test_that("Resample() fails and works", {
  expect_error(Resample(0))
  dataset <- MatrixToPhyDat(rbind(
    a = c(0, 0, 0, 0, 0, 0),
    b = c(0, 0, 0, 0, 0, 0),
    c = c(1, 1, 0, 0, 0, 1),
    d = c(1, 1, 0, 0, 1, 0),
    e = c(1, 1, 1, 1, 1, 1),
    f = c(1, 1, 1, 1, 1, 1)))
  
  expect_error(Resample(dataset, method = "ERROR"))
  expect_error(Resample(dataset, proportion = 0))
  expect_error(Resample(dataset, proportion = 6 / 7))

  nRep <- 42L # Arbitrary number to balance runtime vs false +ves & -ves
  bal <- as.Splits(BalancedTree(dataset))
  set.seed(6034) # Fix seed: stochastic test has ~13% failure rate without it
  
  jackTrees <- replicate(nRep, Resample(dataset, NJTree(dataset), verbosity = 0L),
                         simplify = FALSE)
  jackSplits <- as.Splits(unlist(jackTrees, recursive = FALSE))
  jackSupport <- rowSums(
    # TODO replace :::.in.Splits with exported %in%
    # %in% works when testing file but not entire package
    # See https://github.com/r-lib/testthat/issues/1661
    vapply(jackSplits, function(sp) TreeTools:::.in.Splits(bal, sp), logical(3))
  )
  
  # Stochastic test — tolerance allows for sampling variability
  expect_equal(jackSupport, tolerance = 0.3,
               c("8" = 1/2, "9" = 1, "10" = 1/2, "11" = 0)[names(bal)] *
                 sum(vapply(jackTrees, length, 1L)))
  
  set.seed(4817) # Separate seed so jackknife changes don't shift bootstrap RNG
  bootTrees <- replicate(nRep, Resample(dataset, method = "bootstrap",
                                        verbosity = 0),
                         simplify = FALSE)
  bootSupport <- rowSums(vapply(
    unlist(bootTrees, recursive = FALSE),
    # TODO replace :::.in.Splits with exported %in%
    # %in% works when testing file but not entire package
    # See https://github.com/r-lib/testthat/issues/1661
    function(tr) TreeTools:::.in.Splits(bal, as.Splits(tr)),
    logical(3)
  ))
  expect_equal(bootSupport, tolerance = 0.3,
               c("8" = 1/2, "9" = 1, "10" = 1/2, "11" = 0)[names(bal)] * 
                 sum(vapply(bootTrees, length, 1L)))
    
})

test_that("TreeSearch:::.CombineResults() handles duplicates", {
  x <- structure(
    array(c(
      rep(1L, 8),
      rep(2L, 8),
      rep(3L, 8),
      rep(2L, 8),
      rep(1L, 8)
      ),
      dim = c(4, 2, 5)),
    firstHit = c(start = 5, test = 0, end = 0)
  )
  y <- array(c(rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8),
               rep(4L, 8),
               rep(1L, 8)),
          dim = c(4, 2, 5)
          )
  expect_warning(TreeSearch:::.CombineResults(x, y, stage = "test"))
  uX <- structure(unique(x, MARGIN = 3L),
                  firstHit = c(start = 3, test = 0, end = 0))
  expect_equal(attr(TreeSearch:::.CombineResults(uX, y, stage = "test"), "firstHit"),
               c(start = 3, test = 1, end = 0))
               
})
