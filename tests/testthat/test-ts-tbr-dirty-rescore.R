# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# T-304: enduring regression test for the T-300 dirty-set incremental
# rescore wired into the tbr_search SPR accept path (src/ts_tbr.cpp,
# ~lines 1138-1180).
#
# On an accepted SPR move the search does NOT call full_rescore; instead it
# updates only the nodes on the union of paths nz->root and nx->root via
# fitch_dirty_downpass / fitch_dirty_uppass (EW) or the NA-block variants
# (fitch_na_dirty_*), then derives the score incrementally.  Four code paths
# exist: EW, IW, NA, and NA-IW (is_spr && !has_na | is_spr && has_na, each
# crossed with use_iw).
#
# The DEBUG_RESCORE / DEBUG_NA_RESCORE / DEBUG_NNI_RESCORE cross-checks that
# originally guarded this were removed (commits 5b210fdd, 44a4ebeb, 2be8228d),
# and an earlier incremental attempt regressed with a systematic delta = -3
# and had to be reverted (b7303ee5).  This test is the permanent guard: it
# drives MANY accepted SPR moves (small tips, weak signal, high maxHits) and
# asserts that the score the search reports equals an independent full
# recomputation.  If the dirty-set rescore ever drifts from the authoritative
# score, result$score != ts_score(result_tree, ds) and these fail.

ts_tbr <- function(tree, ds, maxHits = 20L, concavity = -1.0,
                   min_steps = integer(0)) {
  TreeSearch:::ts_tbr_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, min_steps = min_steps,
                             concavity = concavity)
}

result_tree <- function(result, ref_tree) {
  rt <- ref_tree
  rt$edge <- result$edge
  rt
}

test_that("TBR dirty-set rescore matches full rescore (EW, many accepts)", {
  # 12 tips, 6 random multistate characters -> weak signal, so the search
  # accepts a long chain of SPR moves, each exercising the EW dirty-set path.
  set.seed(4471)
  mat <- matrix(sample(0:3, 12 * 6, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (start in c(1, 17, 88, 256, 777)) {
    tree <- as.phylo(start, 12)
    set.seed(1000 + start)
    result <- ts_tbr(tree, ds, maxHits = 50L)

    rt <- result_tree(result, tree)
    independent_score <- ts_score(rt, ds)
    expect_equal(result$score, independent_score,
                 info = paste("EW start =", start))
    validate_result(result, 12L)
  }
})

test_that("TBR dirty-set rescore matches full rescore (IW, many accepts)", {
  set.seed(5529)
  mat <- matrix(sample(0:2, 12 * 6, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  for (start in c(1, 17, 88, 256, 777)) {
    tree <- as.phylo(start, 12)
    set.seed(2000 + start)
    result <- ts_tbr(tree, ds, maxHits = 50L, concavity = 10,
                     min_steps = minSteps)

    rt <- result_tree(result, tree)
    independent_score <- ts_score(rt, ds, concavity = 10, min_steps = minSteps)
    expect_equal(result$score, independent_score, tolerance = 1e-10,
                 info = paste("IW start =", start))
    validate_result(result, 12L)
  }
})

test_that("TBR dirty-set rescore matches full rescore (NA dataset, many accepts)", {
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  for (start in c(1, 42, 100, 314)) {
    tree <- as.phylo(start, n_tip)
    set.seed(3000 + start)
    result <- ts_tbr(tree, ds, maxHits = 20L)

    rt <- result_tree(result, tree)
    independent_score <- ts_score(rt, ds)
    expect_equal(result$score, independent_score,
                 info = paste("NA start =", start))
    validate_result(result, n_tip)
  }
})

test_that("TBR dirty-set rescore matches full rescore (NA-IW dataset, many accepts)", {
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  for (start in c(1, 42, 100, 314)) {
    tree <- as.phylo(start, n_tip)
    set.seed(4000 + start)
    result <- ts_tbr(tree, ds, maxHits = 20L, concavity = 10,
                     min_steps = minSteps)

    rt <- result_tree(result, tree)
    independent_score <- ts_score(rt, ds, concavity = 10, min_steps = minSteps)
    expect_equal(result$score, independent_score, tolerance = 1e-10,
                 info = paste("NA-IW start =", start))
    validate_result(result, n_tip)
  }
})

test_that("dirty-set rescore matches full rescore on TBR-REROOTING accepts", {
  # Issue #38: the dirty-set accept path originally covered SPR-classified
  # accepts only; accepts that rerooted the clipped fragment fell back to
  # full_rescore.  Extending it adds a third dirty seed at clip_node, because
  # apply_tbr_move reverses the parent/child links along
  # clip_node..reroot_parent and so gives every node on that path new children.
  #
  # The four tests above cannot guard this arm: they assert score identity but
  # have no way to tell whether a rerooting accept ever occurred, so they would
  # pass just as happily if the arm were never entered.  `n_reroot_accepts`
  # (src/ts_data.h) is what makes this one non-vacuous -- it is asserted
  # positive, so losing coverage fails the test rather than silently voiding it.
  data("inapplicable.phyData", package = "TreeSearch")
  minSteps <- function(dataset) {
    as.integer(MinimumLength(dataset, compress = TRUE))
  }

  set.seed(6273)
  mat <- matrix(sample(0:3, 20 * 8, replace = TRUE),
                nrow = 20, dimnames = list(paste0("t", 1:20), NULL))
  random20 <- MatrixToPhyDat(mat)
  vinther <- inapplicable.phyData[["Vinther2008"]]

  cases <- list(
    list(label = "EW", dataset = random20, concavity = -1, score_conc = Inf),
    list(label = "IW", dataset = random20, concavity = 10, score_conc = 10),
    list(label = "NA", dataset = vinther, concavity = -1, score_conc = Inf),
    list(label = "NA-IW", dataset = vinther, concavity = 10, score_conc = 10)
  )

  for (case in cases) {
    ds <- make_ts_data(case$dataset)
    n_tip <- length(case$dataset)
    ms <- if (is.finite(case$score_conc)) minSteps(case$dataset) else integer(0)
    n_reroot <- 0

    for (start in c(3, 29, 131, 512, 900)) {
      tree <- as.phylo(start, n_tip)
      set.seed(7000 + start)
      result <- ts_tbr(tree, ds, maxHits = 50L, concavity = case$concavity,
                       min_steps = ms)
      n_reroot <- n_reroot + result$n_reroot_accepts

      rt <- result_tree(result, tree)
      independent <- ts_score(rt, ds, concavity = case$score_conc,
                              min_steps = ms)
      expect_equal(result$score, independent, tolerance = 1e-10,
                   info = paste(case$label, "start =", start))
      validate_result(result, n_tip)
    }

    expect_gt(n_reroot, 0)  # coverage: the rerooting arm was actually entered
  }
})

test_that("XPIWE x4 + dirty-region opts are byte-identical to opts-off (port guard)", {
  # Regression guard for the IW->XPIWE opt port (src/ts_tbr.cpp `iw_family`
  # gate): the x4 reroot batch + extract_char_steps dirty-region must produce
  # byte-identical scores to the opts-off scalar path on the PRODUCTION XPIWE
  # path (MaximizeParsimony defaults to extended IW => ScoringMode::XPIWE).
  # Before the port these opts were gated to plain ScoringMode::IW and so never
  # ran under MaximizeParsimony; this asserts the widening did not perturb
  # XPIWE scores. Requires:
  #   - pure-XPIWE: recode "-"->"?" so has_na = FALSE (the opts are !has_na-gated)
  #   - ratchetCycles >= 3: a perturbation can then fully deactivate a block,
  #     the regime that surfaced the nx_cs/active_mask consistency bug (the
  #     dirty-region's per-clip internal invariant is itself guarded by the C++
  #     TS_IW_DIRTYCHK oracle; this test guards the opts' externally-visible
  #     byte-identity).
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")
  m <- PhyDatToMatrix(inapplicable.phyData[["Vinther2008"]], ambigNA = FALSE)
  m[m == "-"] <- "?"                      # pure-XPIWE: has_na = FALSE
  d <- MatrixToPhyDat(m)
  ctrl <- SearchControl(ratchetCycles = 4L, xssRounds = 0L, rssRounds = 0L,
                        cssRounds = 0L, driftCycles = 0L)
  run <- function(opts_on) {
    if (opts_on) { Sys.unsetenv("TS_IW_NOX4");   Sys.unsetenv("TS_IW_NODIRTY") }
    else         { Sys.setenv(TS_IW_NOX4 = "1"); Sys.setenv(TS_IW_NODIRTY = "1") }
    set.seed(909)
    r <- suppressWarnings(MaximizeParsimony(
      d, concavity = 10, maxReplicates = 1L, nThreads = 1L,
      verbosity = 0L, control = ctrl))
    min(attr(r, "score"))
  }
  withr::defer({ Sys.unsetenv("TS_IW_NOX4"); Sys.unsetenv("TS_IW_NODIRTY") })
  score_on  <- run(TRUE)
  score_off <- run(FALSE)
  expect_equal(score_on, score_off, tolerance = 0)
})

test_that("NA-IW x4 reroot batch is byte-identical to scalar (NA port guard)", {
  # Regression guard for indirect_na_iw_cached_flat_x4 (src/ts_fitch.cpp) wired
  # into the iw_family scan branch (src/ts_tbr.cpp): the 4-wide NA-IW reroot
  # batch must produce byte-identical scores to the one-at-a-time scalar
  # indirect_na_iw_length_cached on the IW+NA / XPIWE+NA path.  This is the
  # complement of the port guard above: that one RECODES "-"->"?" (has_na =
  # FALSE, exercising the no-NA IW x4); this one KEEPS the inapplicable
  # characters (has_na = TRUE) so the NA-aware batch kernel actually fires.
  # TS_IW_NOX4 toggles the batch off (falls through to the scalar `else`).
  skip_if_not_installed("TreeSearch")
  data("inapplicable.phyData", package = "TreeSearch")
  d <- inapplicable.phyData[["Vinther2008"]]   # keep "-"; has_na = TRUE
  ctrl <- SearchControl(ratchetCycles = 4L, xssRounds = 0L, rssRounds = 0L,
                        cssRounds = 0L, driftCycles = 0L)
  run <- function(x4_on) {
    if (x4_on) Sys.unsetenv("TS_IW_NOX4") else Sys.setenv(TS_IW_NOX4 = "1")
    set.seed(717)
    r <- suppressWarnings(MaximizeParsimony(
      d, concavity = 10, maxReplicates = 1L, nThreads = 1L,
      verbosity = 0L, control = ctrl))
    min(attr(r, "score"))
  }
  withr::defer(Sys.unsetenv("TS_IW_NOX4"))
  score_x4  <- run(TRUE)
  score_scalar <- run(FALSE)
  expect_equal(score_x4, score_scalar, tolerance = 0)
})
