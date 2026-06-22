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
  on.exit({ Sys.unsetenv("TS_IW_NOX4"); Sys.unsetenv("TS_IW_NODIRTY") }, add = TRUE)
  score_on  <- run(TRUE)
  score_off <- run(FALSE)
  expect_equal(score_on, score_off, tolerance = 0)
})
