library("TreeTools", quietly = TRUE)

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

# --- Input validation ---

test_that("MaximizeParsimony stops with message when dataset is NULL", {
  expect_error(
    MaximizeParsimony(NULL, maxReplicates = 1L, targetHits = 1L,
                      verbosity = 0L),
    "`dataset` cannot be NULL."
  )
})

test_that("MaximizeParsimony rejects maxReplicates < 1 (T-341)", {
  # maxReplicates = 0 runs the search loop zero times, leaving the pool
  # empty and best_score at the C++ sentinel of -1; without this guard, the
  # empty-pool fallback silently returned the random starting tree tagged
  # with that bogus score instead of erroring.
  expect_error(
    MaximizeParsimony(ds, maxReplicates = 0L, targetHits = 1L,
                      verbosity = 0L),
    "`maxReplicates` must be"
  )
  expect_error(
    MaximizeParsimony(ds, maxReplicates = -1L, targetHits = 1L,
                      verbosity = 0L),
    "`maxReplicates` must be"
  )
  expect_error(
    MaximizeParsimony(ds, maxReplicates = NA_integer_, targetHits = 1L,
                      verbosity = 0L),
    "`maxReplicates` must be"
  )
})

test_that("replicate-adequacy warning uses unscaled character count (T-342)", {
  # `weight` is the .ScaleWeight()-integerised value (up to ~1260x for
  # fractional weights); the printed `nChars` must reflect the true number
  # of characters, not that internal scale factor. The warning only fires
  # for nTip >= 30, so use a synthetic dataset large enough to trigger it.
  set.seed(1)
  dat <- TreeTools::MatrixToPhyDat(matrix(
    sample(0:1, 30 * 10, replace = TRUE), nrow = 30,
    dimnames = list(paste0("t", 1:30), NULL)))
  attr(dat, "weight") <- rep(0.5, attr(dat, "nr"))
  nCharsTrue <- sum(attr(dat, "weight"))
  expect_warning(
    MaximizeParsimony(dat, maxReplicates = 1L, targetHits = 1L,
                      verbosity = 1L),
    paste0(nCharsTrue, " characters")
  )
})

# --- Strategy presets ---

test_that("the sprint rung runs and returns a valid result", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(3418)
  result <- MaximizeParsimony(ds, effort = -9L,
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_true(attr(result, "score") > 0)
  expect_equal(NTip(result[[1]]), NTip(ds))
})

test_that("candidates_evaluated attribute is reported for serial search", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  # Diagnostic counter (TNT "rearrangements examined" analogue): a positive,
  # finite scalar for a single-threaded search. See MaximizeParsimony @return.
  set.seed(3418)
  result <- MaximizeParsimony(ds, effort = -9L,
                               maxReplicates = 2L, targetHits = 1L,
                               nThreads = 1L, verbosity = 0L)
  ce <- attr(result, "candidates_evaluated")
  expect_type(ce, "double")
  expect_length(ce, 1L)
  expect_true(is.finite(ce) && ce > 0)
})

test_that("the default rung runs and returns a valid result", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(5726)
  result <- MaximizeParsimony(ds, .rung = "default",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("the thorough rung runs and returns a valid result", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(8103)
  result <- MaximizeParsimony(ds, .rung = "thorough",
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that(".AutoRung selects on size and signal density", {
  AS <- TreeSearch:::.AutoRung

  # Small datasets: always sprint
  expect_equal(AS(20, 100), 1L)
  expect_equal(AS(30, 500), 1L)

  # Few chars (< 100 patterns) -> flat landscape -> always default
  expect_equal(AS(50, 25),   2L)  # small, very few chars
  expect_equal(AS(65, 80),   2L)  # large enough tip count, but nChar < 100
  expect_equal(AS(200, 99),  2L)  # large, but still nChar < 100

  # Mid-size (31-64 tips) with enough chars -> default (not large enough)
  expect_equal(AS(60, 300), 2L)  # nChar >= 100 but nTip < 65
  expect_equal(AS(64, 200), 2L)  # nChar >= 100 but nTip < 65

  # Large (>= 65 tips) with enough chars -> thorough
  # Signal density does NOT gate thorough: more chars = more benefit (T-068 benchmark)
  expect_equal(AS(65, 100),   3L)  # boundary case: 65 tips, 100 chars
  expect_equal(AS(74, 200),   3L)  # 74 tips, ratio 2.7
  expect_equal(AS(75, 250),   3L)  # ratio 3.3
  expect_equal(AS(100, 200),  3L)  # ratio 2.0
  expect_equal(AS(75, 400),   3L)  # ratio 5.3 — high ratio still benefits
  expect_equal(AS(119, 2800), 3L)  # just below large threshold
  expect_equal(AS(125, 2800), 4L)    # >= 120 tips -> large
  expect_equal(AS(200, 100),  4L)    # >= 120 tips -> large
  expect_equal(AS(200, 1200), 4L)    # >= 120 tips -> large
})

test_that("effort = 0 selects the rung from dataset size", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(2944)
  # Vinther2008 has 23 tips -> should auto-select "sprint"
  result <- MaximizeParsimony(ds, effort = 0L,
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("the no-preset escape uses raw parameter defaults", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(6017)
  result <- MaximizeParsimony(ds, .rung = "none",
                               maxReplicates = 2L, targetHits = 1L,
                               ratchetCycles = 1L, driftCycles = 0L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("explicit params override the rung preset", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(1589)
  # Sprint has driftCycles=0; override to 1
  result <- MaximizeParsimony(ds, effort = -9L,
                               driftCycles = 1L,
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("`effort` rejects a non-integer offset", {
  # `effort` is an OFFSET, so there is no vocabulary of names to mistype and no
  # "unknown strategy" fallback to warn about.  What can go wrong instead is a
  # value that is not a whole number, which must error rather than silently
  # truncate -- a search quietly run at the wrong rung is worse than a stop.
  expect_error(MaximizeParsimony(ds, effort = 1.5, maxReplicates = 1L,
                                 verbosity = 0L),
               "whole number")
  expect_error(MaximizeParsimony(ds, effort = "thorough", maxReplicates = 1L,
                                 verbosity = 0L),
               "whole number")
  expect_error(MaximizeParsimony(ds, effort = c(1L, 2L), maxReplicates = 1L,
                                 verbosity = 0L),
               "whole number")
})

test_that("`effort` clamps at both ends of the ladder", {
  AR <- TreeSearch:::.AutoRung
  ER <- TreeSearch:::.EffortRung
  # Bottom clamp is load-bearing: it is what lets a caller write a large
  # negative offset and reliably get `sprint` whatever the dataset.
  expect_equal(ER(AR(200L, 200L), -99L, 0L), 1L)
  expect_equal(ER(AR(20L, 100L), -1L, 0L), 1L)
  # The only top limit is representability, and it announces itself rather than
  # silently pretending a bigger number meant something.  There is deliberately
  # no policy ceiling below it: extra replicates cost wall but cannot cost
  # reach, so refusing to go further would just obstruct the request.
  expect_equal(ER(1L, 99L, 0L), TreeSearch:::.effortMaxRung)
  expect_message(ER(1L, 99L, 1L), "clamped to rung")
  # The ceiling must be exactly where the budget stops being an integer: one
  # rung lower would be arbitrary, one higher would overflow.
  RS <- TreeSearch:::.RungSpec
  expect_true(RS(TreeSearch:::.effortMaxRung)[["maxReplicates"]] > 0L)
  expect_true(is.na(suppressWarnings(
    as.integer(500 * 2^(TreeSearch:::.effortMaxRung + 1L - 4L)))))
})

test_that("both budget knobs double, so a notch is the same size either way", {
  # Mixed rates would make one notch 2x the work on hard datasets (where
  # maxReplicates binds) but only (k+1)/k on easy ones (where targetHits does),
  # so notches would shrink as you climb on the easy population.
  RS <- TreeSearch:::.RungSpec
  for (r in 5:9) {
    expect_equal(RS(r)[["maxReplicates"]] / RS(r - 1L)[["maxReplicates"]], 2)
    expect_equal(RS(r)[["hitMultiplier"]] / RS(r - 1L)[["hitMultiplier"]], 2)
  }
  # Rungs 1-4 leave both alone: they differ in provisioning, not budget.
  for (r in 1:4) expect_equal(RS(r)[["hitMultiplier"]], 1L)
})

test_that("effort = 0 reproduces the automatic choice on every size band", {
  # The whole point of an offset rather than an absolute level: the default
  # must be byte-identical to what the package chose before `effort` existed.
  AR <- TreeSearch:::.AutoRung
  ER <- TreeSearch:::.EffortRung
  for (nTip in c(20L, 50L, 70L, 200L)) {
    for (nChar in c(50L, 300L)) {
      expect_equal(ER(AR(nTip, nChar), 0L, 0L), AR(nTip, nChar))
    }
  }
})

# --- maxSeconds timeout ---

test_that("maxSeconds stops search before maxReplicates", {
  set.seed(7392)
  # Use a near-zero timeout so even a small dataset triggers it reliably.
  # The timeout is checked between replicates, so the first may complete

  # but subsequent ones won't start.
  result <- MaximizeParsimony(ds, maxReplicates = 1000L, targetHits = 1000L,
                               maxSeconds = 0.001, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_lt(attr(result, "replicates"), 1000L)
})

test_that("maxSeconds = 0 means no timeout", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(8456)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               maxSeconds = 0, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_false(attr(result, "timed_out"))
})

test_that("verbosity = 1 prints 'Search complete' summary to console", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(3071)
  # MaximizeParsimony emits two streams at verbosity = 1: cli messages
  # via message() ("Strategy: ...", "Search complete: ...") and C++
  # Rprintf progress via stdout ("Replicate N/M", "Converged: ...").
  # Capture both so they don't leak into testthat output, then assert
  # that the expected lines were produced.
  msg_lines <- character()
  stdout_lines <- capture.output(
    msg_lines <- capture.output(
      MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                        verbosity = 1L),
      type = "message"
    )
  )
  expect_true(any(grepl("Search complete", msg_lines)))
  expect_true(any(grepl("Replicate", stdout_lines)))
  expect_true(any(grepl("Converged|score", stdout_lines)))
})

# --- nThreads ---

test_that("nThreads = 1 (serial) runs correctly", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(5193)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               nThreads = 1L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
})

test_that("nThreads = 2 (parallel) runs correctly", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(6274)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               nThreads = 2L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  # Score should be reasonable (not garbage from parallel corruption)
  expect_true(attr(result, "score") < 200)
})

# --- User-supplied starting tree (warm-start) ---

test_that("user tree is used as warm start", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  # Build a known tree
  set.seed(9847)
  user_tree <- RandomTree(ds, root = TRUE)
  user_tree <- Preorder(user_tree)

  result <- MaximizeParsimony(ds, tree = user_tree,
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  result_score <- attr(result, "score")

  # Result should be at least as good as the input tree
  input_score <- TreeLength(user_tree, ds)
  expect_true(result_score <= input_score)
})

test_that("multiPhylo warm starts survive tip renumbering and polytomies", {
  skip_on_cran()
  # Trees whose tip order differs from the dataset, one of them unresolved:
  # every start must be normalized (renumbered, resolved, rerooted) before it
  # reaches the engine, not just the first.
  set.seed(1002)
  shuffled <- RandomTree(ds, root = TRUE)
  shuffled <- KeepTip(shuffled, rev(TipLabels(shuffled)))
  polytomous <- CollapseNode(Preorder(RandomTree(ds, root = TRUE)),
                             NTip(ds) + 3L)
  trees <- structure(list(shuffled, polytomous), class = "multiPhylo")

  set.seed(4004)
  result <- MaximizeParsimony(ds, tree = trees, maxReplicates = 2L,
                              targetHits = 99L, verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_equal(attr(result, "score"), TreeLength(result[[1]], ds))
  expect_setequal(TipLabels(result[[1]]), names(ds))
})

test_that("multiPhylo warm starts are validated", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(6003)
  t1 <- Preorder(RandomTree(ds, root = TRUE))
  t2 <- Preorder(RandomTree(ds, root = TRUE))
  t2[["tip.label"]][[1]] <- "not_a_taxon"

  expect_error(
    MaximizeParsimony(ds, tree = structure(list(), class = "multiPhylo"),
                      maxReplicates = 1L, verbosity = 0L),
    "contains no trees"
  )
  expect_error(
    MaximizeParsimony(ds, tree = structure(list(t1, t2), class = "multiPhylo"),
                      maxReplicates = 1L, verbosity = 0L),
    "same tip labels"
  )
  expect_error(
    MaximizeParsimony(ds, tree = structure(list(t1, "not a tree"),
                                           class = "multiPhylo"),
                      maxReplicates = 1L, verbosity = 0L),
    "class 'phylo'"
  )

  # A structurally invalid `phylo` must be rejected in R.  ape::unroot()
  # manufactures one from any TreeTools `order = "preorder"` tree, and passing
  # it on segfaults inside TreeTools' rooting code -- unrecoverable, so the
  # error has to come first.  The index tells the user which tree is at fault.
  broken <- ape::unroot(t1)
  expect_error(
    MaximizeParsimony(ds, tree = broken, maxReplicates = 1L, verbosity = 0L),
    "`tree` is not a valid tree"
  )
  expect_error(
    MaximizeParsimony(ds, tree = structure(list(t1, broken),
                                           class = "multiPhylo"),
                      maxReplicates = 1L, verbosity = 0L),
    "`tree\\[\\[2\\]\\]` is not a valid tree"
  )
  # A genuinely unrooted tree is fine; only the corrupted object is refused.
  expect_silent(
    MaximizeParsimony(ds, tree = RandomTree(ds, root = FALSE),
                      maxReplicates = 1L, targetHits = 1L, verbosity = 0L)
  )

  # A validly unrooted tree from outside TreeTools (e.g. ape::rtree() then
  # ape::unroot(), unlike `broken` above) has nrow(edge) == 2 * n - 3, which
  # fails the "is this already bifurcating?" test the same way a genuine
  # polytomy would, because that test conflates "needs resolving" with "needs
  # rooting".  Passing it to MakeTreeBinary() misreads the unrooted root's
  # legitimate degree-3 trifurcation as a polytomy, corrupting the tree and
  # previously surfacing downstream as "argument is of length zero".
  set.seed(9)
  apeUnrooted <- ape::unroot(ape::rtree(NTip(ds), tip.label = names(ds)))
  expect_silent(
    MaximizeParsimony(ds, tree = apeUnrooted,
                      maxReplicates = 1L, targetHits = 1L, verbosity = 0L)
  )

  # Unused pool members are reported against the replicates actually run,
  # which targetHits can cut short well below maxReplicates.
  expect_warning(
    MaximizeParsimony(ds, tree = structure(list(t1, t1, t1),
                                           class = "multiPhylo"),
                      maxReplicates = 9L, targetHits = 1L, verbosity = 0L),
    "of the 3 trees supplied"
  )
})

# --- timings attribute ---

test_that("timings attribute is returned", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(2689)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  timings <- attr(result, "timings")
  expect_false(is.null(timings))
  expect_true(is.numeric(timings))
  expect_true(all(timings >= 0))
  expect_true("wagner_ms" %in% names(timings))
  expect_true("ratchet_ms" %in% names(timings))
})

# --- IW with strategy ---

test_that("IW mode works with effort rungs", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(4012)
  result <- MaximizeParsimony(ds, concavity = 10, effort = -9L,
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  cpp_score <- attr(result, "score")
  tl_score <- TreeLength(result[[1]], ds, concavity = 10)
  expect_equal(cpp_score, tl_score, tolerance = 0.01)
})

# --- T-340 regression: `concavity` normalization ---

test_that("concavity as a numeric-coercible string behaves like the number", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  # A bare "10" must not silently drop into IW mode with unpopulated
  # min_steps (uncorrected homoplasy) -- it should match concavity = 10
  # exactly, seed-for-seed.
  set.seed(4012)
  asString <- MaximizeParsimony(ds, concavity = "10", effort = -9L,
                                 maxReplicates = 2L, targetHits = 1L,
                                 verbosity = 0L)
  set.seed(4012)
  asNumber <- MaximizeParsimony(ds, concavity = 10, effort = -9L,
                                 maxReplicates = 2L, targetHits = 1L,
                                 verbosity = 0L)
  expect_equal(attr(asString, "score"), attr(asNumber, "score"))
  expect_equal(asString[[1]][["edge"]], asNumber[[1]][["edge"]])
})

test_that("concavity = 'Profile'/'prof' route to profile mode like 'profile'", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(4012)
  canonical <- MaximizeParsimony(ds, concavity = "profile", effort = -9L,
                                  maxReplicates = 2L, targetHits = 1L,
                                  verbosity = 0L)
  for (spelling in c("Profile", "prof")) {
    set.seed(4012)
    result <- MaximizeParsimony(ds, concavity = spelling, effort = -9L,
                                 maxReplicates = 2L, targetHits = 1L,
                                 verbosity = 0L)
    expect_equal(attr(result, "score"), attr(canonical, "score"))
  }
})

test_that("an invalid concavity string errors cleanly instead of silently using EW", {
  expect_error(
    MaximizeParsimony(ds, concavity = "banana", maxReplicates = 1L,
                       targetHits = 1L, verbosity = 0L),
    "`concavity` must be a single positive number"
  )
})

# --- Output tree validity ---

test_that("output trees have valid preorder numbering", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(8734)
  result <- MaximizeParsimony(ds, maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  for (tree in result) {
    # Trees should have correct number of tips and edges
    expect_equal(NTip(tree), NTip(ds))
    expect_equal(nrow(tree$edge), 2 * (NTip(ds) - 1))
    # Tip labels should match
    expect_true(all(TipLabels(tree) %in% names(ds)))
  }
})

# --- T-039 regression: constraint on small fully-resolved trees ---

test_that("Fully-resolving constraint on 5-tip tree does not crash", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(4172)
  result <- MaximizeParsimony(ds5, constraint = cons,
                               maxReplicates = 1L, targetHits = 1L,
                               verbosity = 0L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_equal(NTip(result[[1]]), 5L)
})

test_that("AdditionTree with fully-resolving constraint works", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  set.seed(6091)
  wt <- AdditionTree(ds5, constraint = cons)
  expect_s3_class(wt, "phylo")
  expect_equal(NTip(wt), 5L)
  expect_equal(nrow(wt$edge), 8L)
})

test_that("Constrained Wagner tree works with multiple seeds", {
  ds5 <- phangorn::phyDat(
    matrix(c("0","0","0","1","1","0","1","0","1","0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1"))
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")
  at <- attributes(ds5)
  consArgs <- TreeSearch:::.PrepareConstraint(cons, ds5)

  for (s in c(1, 6, 42)) {
    set.seed(s)
    result <- do.call(TreeSearch:::ts_random_wagner_tree, c(
      list(contrast = at$contrast,
           tip_data = matrix(unlist(ds5, use.names = FALSE), nrow = 5, byrow = TRUE),
           weight = at$weight, levels = at$levels),
      consArgs))
    expect_true(is.finite(result$score), info = paste("seed", s))
    expect_equal(nrow(result$edge), 8L, info = paste("seed", s))
  }
})

# --- Intra-replicate fusing (T-258) ---

test_that("intraFuse runs without error", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  set.seed(8517)
  result <- MaximizeParsimony(ds, effort = -9L,
                              maxReplicates = 5L, targetHits = 2L,
                              maxSeconds = 3, intraFuse = TRUE,
                              verbosity = 0L, nThreads = 1L)
  expect_s3_class(result, "multiPhylo")
  expect_true(is.finite(attr(result, "score")))
  expect_lte(attr(result, "score"), 100)  # should find reasonable score
})

test_that("intraFuse with dataset size change does not crash", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  ds_large <- inapplicable.phyData[["Agnarsson2004"]]  # 62 tips
  ds_small <- inapplicable.phyData[["Vinther2008"]]     # 23 tips

  # Run on larger dataset first with intra-fuse
  set.seed(9014)
  r1 <- MaximizeParsimony(ds_large, effort = -9L,
                          maxReplicates = 3L, targetHits = 2L,
                          maxSeconds = 3, intraFuse = TRUE,
                          verbosity = 0L, nThreads = 1L)
  expect_true(is.finite(attr(r1, "score")))

  # Then run on smaller dataset with intra-fuse (regression test for segfault)
  set.seed(9015)
  r2 <- MaximizeParsimony(ds_small, effort = -9L,
                          maxReplicates = 3L, targetHits = 2L,
                          maxSeconds = 3, intraFuse = TRUE,
                          verbosity = 0L, nThreads = 1L)
  expect_true(is.finite(attr(r2, "score")))
})

# --- collapse = TRUE (TNT-style zero-length-branch contraction) ---

# A soft polytomy: (b1,b2) & (b3,b4) supported cherries plus a fan of identical
# taxa that has no internal support.  Every resolution of the fan is an MPT, so
# collapse = FALSE returns many fully-resolved trees while collapse = TRUE must
# contract the fan to a single polytomy -> exactly one distinct collapsed tree.
.SoftPolytomyData <- function(nFan = 5L) {
  backbone <- c("b1", "b2", "b3", "b4")
  fan <- paste0("f", seq_len(nFan))
  taxa <- c(backbone, fan)
  mat <- matrix("0", nrow = length(taxa), ncol = 0,
                dimnames = list(taxa, NULL))
  addChar <- function(m, ones) {
    cbind(m, ifelse(rownames(m) %in% ones, "1", "0"))
  }
  for (i in 1:3) mat <- addChar(mat, fan)            # fan synapomorphy x3
  for (i in 1:2) mat <- addChar(mat, c("b1", "b2"))  # backbone cherry x2
  for (i in 1:2) mat <- addChar(mat, c("b3", "b4"))  # backbone cherry x2
  phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))
}

test_that("collapse = TRUE contracts a soft polytomy to one collapsed tree", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  phy <- .SoftPolytomyData(5L)
  set.seed(1L)
  resolved <- MaximizeParsimony(
    phy, concavity = Inf, maxReplicates = 12L, .rung = "thorough",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = FALSE)
  set.seed(1L)
  collapsed <- MaximizeParsimony(
    phy, concavity = Inf, maxReplicates = 12L, .rung = "thorough",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = TRUE)

  # Many resolved variants of one collapsed topology.
  expect_gt(length(resolved), 1L)
  # All collapse to a single distinct topology.
  expect_equal(length(collapsed), 1L)
  expect_equal(attr(collapsed, "n_topologies"), 1L)
  # The single tree is a genuine polytomy (fewer internal nodes than binary).
  expect_lt(collapsed[[1]][["Nnode"]], NTip(phy) - 1L)
  # The supported backbone cherries survive collapse (the fan, not the
  # backbone, is contracted).  Root on a fan tip so neither cherry sits at the
  # root, then each cherry is the exclusive descendant set of its MRCA.
  rooted <- TreeTools::RootTree(collapsed[[1]], "f1")
  isClade <- function(members) {
    mrca <- ape::getMRCA(rooted, members)
    !is.null(mrca) &&
      setequal(rooted[["tip.label"]][phangorn::Descendants(rooted, mrca,
                                                           "tips")[[1]]],
               members)
  }
  expect_true(isClade(c("b1", "b2")))
  expect_true(isClade(c("b3", "b4")))
})

test_that("collapse = TRUE is a no-op when no branch is unsupported", {
  skip_on_cran() # mixed-tier exception, see tests/testing-strategy.md
  # Vinther2008 MPTs are fully resolved (no zero-length branches): collapse must
  # leave every tree binary and topologically unchanged.
  set.seed(3418)
  resolved <- MaximizeParsimony(ds, effort = -9L,
                                maxReplicates = 3L, targetHits = 1L,
                                verbosity = 0L, collapse = FALSE)
  set.seed(3418)
  collapsed <- MaximizeParsimony(ds, effort = -9L,
                                 maxReplicates = 3L, targetHits = 1L,
                                 verbosity = 0L, collapse = TRUE)
  nTip <- NTip(ds)
  # No spurious collapse: every returned tree stays fully binary.
  expect_true(all(vapply(collapsed, function(t) t[["Nnode"]] == nTip - 1L,
                         logical(1))))
  # Same set of topologies (each collapsed tree matches a resolved tree, RF 0).
  expect_true(all(apply(
    as.matrix(TreeDist::RobinsonFoulds(collapsed, resolved)), 1L, min) == 0))
})

test_that("collapse = TRUE shows an enforced clade and collapses the rest", {
  # "Show the enforced clade": a constraint encodes external evidence for a
  # grouping the matrix does not support, so it must stay visible even though it
  # sits on a zero-length branch.  Forcing (f1, f2) inside the otherwise
  # unsupported fan gives it a zero-length branch; collapse must keep it while
  # still contracting the unconstrained fan resolutions (the philosophy that
  # makes a constrained collapse meaningful, vs the old skip-under-constraint).
  phy <- .SoftPolytomyData(5L)
  constraint <- ape::read.tree(
    text = "((f1, f2), (b1, b2, b3, b4, f3, f4, f5));")
  set.seed(1L)
  collapsed <- MaximizeParsimony(
    phy, constraint = constraint, concavity = Inf, maxReplicates = 12L,
    .rung = "thorough",
    control = SearchControl(poolMaxSize = 2000L, tbrMaxHits = 200L),
    verbosity = 0L, collapse = TRUE)

  isClade <- function(tr, members) {
    tr <- TreeTools::RootTree(tr, "b1")
    mrca <- ape::getMRCA(tr, members)
    !is.null(mrca) &&
      setequal(tr[["tip.label"]][phangorn::Descendants(tr, mrca, "tips")[[1]]],
               members)
  }
  # The enforced (f1, f2) clade is shown in every returned tree ...
  expect_true(all(vapply(collapsed, isClade, logical(1),
                         members = c("f1", "f2"))))
  # ... yet the rest of the fan is still contracted (tree is a polytomy), so the
  # collapse acted everywhere except on the protected constraint split.
  expect_true(all(vapply(collapsed, function(t)
    t[["Nnode"]] < NTip(phy) - 1L, logical(1))))
})

test_that("collapse = TRUE contracts a fully-uninformative dataset to a star (T-331)", {
  # When every character is parsimony-uninformative (all-constant here), every
  # binary resolution ties at the same (zero) score, so the collapsed answer is
  # a single star: one polytomy, one topology -- never the inflated binary
  # count that a total_words == 0 no-op would produce.
  taxa <- paste0("t", 1:5)
  mat <- matrix("0", nrow = length(taxa), ncol = 2,
                dimnames = list(taxa, NULL))
  phy <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))

  collapsed <- MaximizeParsimony(phy, verbosity = 0L, collapse = TRUE)

  expect_equal(length(collapsed), 1L)
  expect_equal(attr(collapsed, "n_topologies"), 1L)
  expect_equal(collapsed[[1]][["Nnode"]], 1L)
})
