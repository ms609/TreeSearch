# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# =========================================================================
# Degenerate dataset shapes, driven through every dataset-taking entry point
# =========================================================================
# These tests exist to be RUN UNDER A SANITIZER, not for their expectations.
# Two CI legs give the coverage:
#
#   * `glibcxx-assertions` (agent-check.yml) builds with
#     `-D_GLIBCXX_ASSERTIONS`, which aborts on out-of-range container
#     *address formation* -- `&v[0]` on an empty vector, `v.back()` on an
#     empty vector -- with no load or store required.  Plain builds and
#     ASan are both blind to that.
#   * `gcc-ASAN` runs UBSan, which reports `memcpy`'s `nonnull` parameters
#     receiving the null `.data()` of an empty vector.  Hardened libstdc++
#     is blind to *that*.
#
# Every historical instance of both classes was reached by a dataset that
# left the Fitch kernel nothing to do -- `DataSet::total_words == 0` and
# `n_blocks == 0`, so the per-word state vectors are empty:
# agent-issues/TreeSearch#51 (wagner_tree), #60 (no assertions leg at all),
# #124 (TbrSnapshot save/restore), #151 (ts_bench_tbr_phases).  Each was
# found by a sanitizer stumbling into the shape, never by a test aimed at
# it; this file aims at it.
#
# The expectations below are contract checks, deliberately shallow: an
# entry point that swallows a degenerate dataset silently is as much a bug
# as one that aborts, so each call asserts something the API *promises*.
# The deep check is that the process is still alive to run them.

Rp <- function(x, n) rep(x, n)

DegenerateDat <- function(nTip, cols) {
  mat <- do.call(cbind, cols)
  rownames(mat) <- paste0("t", seq_len(nTip))
  MatrixToPhyDat(mat)
}

# Bridge arguments in the form every ts_* entry point takes.
DatBits <- function(dataset) {
  at <- attributes(dataset)
  list(contrast = at$contrast,
       tipData = matrix(unlist(dataset, use.names = FALSE),
                        nrow = length(dataset), byrow = TRUE),
       weight = at$weight, levels = at$levels)
}

nTip <- 8L
zeroWordSets <- list(
  allConstant = DegenerateDat(
    nTip, list(Rp("0", nTip), Rp("0", nTip), Rp("0", nTip))),
  allConstantOnes = DegenerateDat(nTip, list(Rp("1", nTip), Rp("1", nTip))),
  allAutapomorphic = DegenerateDat(
    nTip, list(c("1", Rp("0", nTip - 1)), c(Rp("0", nTip - 1), "1"))),
  allInapplicable = DegenerateDat(nTip, list(Rp("-", nTip), Rp("-", nTip))),
  ambiguousPlusConstant = DegenerateDat(
    nTip, list(Rp("?", nTip), Rp("0", nTip))),
  minimumTips = DegenerateDat(3L, list(Rp("0", 3L), Rp("0", 3L)))
)


test_that("uninformative datasets really do leave the Fitch kernel empty", {
  # The premise the rest of this file rests on.  If simplification ever
  # stops collapsing these to zero blocks, the shapes below stop testing
  # what they were written to test, and this expectation says so loudly
  # rather than letting the file quietly go vacuous.
  for (nm in names(zeroWordSets)) {
    dataset <- zeroWordSets[[nm]]
    tree <- FixedTree(dataset, 1L)
    bits <- DatBits(dataset)
    phases <- TreeSearch:::ts_bench_tbr_phases(
      tree$edge, bits$contrast, bits$tipData, bits$weight, bits$levels)
    expect_equal(phases$total_words, 0L, info = nm)
    expect_equal(phases$n_blocks, 0L, info = nm)
  }
})


test_that("zero-Fitch-word datasets score correctly", {
  # A zero-word dataset is not necessarily a zero-score dataset:
  # autapomorphies are dropped from the kernel but each still contributes
  # its one inevitable step to the total.
  Score <- function(dataset) {
    bits <- DatBits(dataset)
    TreeSearch:::ts_tbr_search(FixedTree(dataset, 1L)$edge, bits$contrast,
                               bits$tipData, bits$weight, bits$levels)$score
  }
  expect_equal(Score(zeroWordSets$allConstant), 0)
  expect_equal(Score(zeroWordSets$allInapplicable), 0)
  expect_equal(Score(zeroWordSets$ambiguousPlusConstant), 0)
  # Two autapomorphic characters, one step each, on any topology.
  expect_equal(Score(zeroWordSets$allAutapomorphic), 2)

  autapomorphic <- zeroWordSets$allAutapomorphic
  tree <- FixedTree(autapomorphic, 1L)
  expect_equal(TreeLength(tree, autapomorphic), 2)
  expect_equal(unname(CharacterLength(tree, autapomorphic)), c(1, 1))
  expect_equal(unname(CharacterLength(FixedTree(zeroWordSets$allConstant, 1L),
                                      zeroWordSets$allConstant)),
               c(0, 0, 0))
})


test_that("every dataset-taking entry point survives zero Fitch words", {
  for (nm in names(zeroWordSets)) {
    dataset <- zeroWordSets[[nm]]
    tree <- FixedTree(dataset, 1L)
    edge <- tree$edge
    bits <- DatBits(dataset)
    co <- bits$contrast
    tp <- bits$tipData
    wt <- bits$weight
    lv <- bits$levels
    nTaxa <- length(dataset)

    for (concavity in c(-1, 10)) {
      expect_type(TreeSearch:::ts_fitch_score(edge, co, tp, wt, lv,
                                              concavity = concavity), "double")
      expect_type(TreeSearch:::ts_tbr_search(edge, co, tp, wt, lv,
                                             concavity = concavity), "list")
      expect_type(TreeSearch:::ts_spr_search(edge, co, tp, wt, lv,
                                             concavity = concavity), "list")
      expect_type(TreeSearch:::ts_nni_search(edge, co, tp, wt, lv,
                                             concavity = concavity), "list")
      expect_type(TreeSearch:::ts_ratchet_search(edge, co, tp, wt, lv,
                                                 nCycles = 2L,
                                                 concavity = concavity), "list")
      expect_type(TreeSearch:::ts_drift_search(edge, co, tp, wt, lv,
                                               nCycles = 2L,
                                               concavity = concavity), "list")
      expect_type(TreeSearch:::ts_rss_search(edge, co, tp, wt, lv,
                                             ratchetCycles = 2L,
                                             concavity = concavity), "list")
      expect_type(TreeSearch:::ts_xss_search(edge, co, tp, wt, lv,
                                             ratchetCycles = 2L,
                                             concavity = concavity), "list")
      expect_type(TreeSearch:::ts_tbr_diagnostics(edge, co, tp, wt, lv,
                                                  concavity = concavity),
                  "list")
      expect_type(TreeSearch:::ts_wagner_tree(co, tp, wt, lv,
                                              concavity = concavity), "list")
      expect_type(TreeSearch:::ts_random_wagner_tree(co, tp, wt, lv,
                                                     concavity = concavity),
                  "list")
      expect_type(TreeSearch:::ts_tree_fuse(edge, co, tp, wt, lv,
                                            pool_edges = list(edge, edge),
                                            pool_scores = c(0, 0),
                                            max_rounds = 2L,
                                            concavity = concavity), "list")
      expect_type(TreeSearch:::ts_resample_search(co, tp, wt, lv,
                                                  maxReplicates = 2L,
                                                  ratchetCycles = 1L,
                                                  concavity = concavity),
                  "list")
      expect_type(TreeSearch:::ts_successive_approx(co, tp, wt, lv,
                                                    maxSAIter = 2L,
                                                    maxReplicates = 2L,
                                                    ratchetCycles = 1L,
                                                    concavity = concavity),
                  "list")
      expect_type(TreeSearch:::ts_parallel_resample(co, tp, wt, lv,
                                                    nReplicates = 2L,
                                                    nThreads = 2L,
                                                    maxReplicates = 2L,
                                                    ratchetCycles = 1L,
                                                    concavity = concavity),
                  "list")
    }

    expect_type(TreeSearch:::ts_char_steps(edge, co, tp, wt, lv), "integer")
    expect_type(TreeSearch:::ts_na_char_steps(edge, co, tp, wt, lv), "list")
    expect_type(TreeSearch:::ts_simplify_diag(co, tp, wt, lv), "list")
    expect_type(TreeSearch:::ts_collapsed_flags_debug(edge, co, tp, wt, lv,
                                                      TRUE), "list")
    expect_type(TreeSearch:::ts_collapsed_flags_debug(edge, co, tp, wt, lv,
                                                      FALSE), "list")
    expect_type(TreeSearch:::ts_debug_clip(edge, co, tp, wt, lv, nTaxa + 2L),
                "list")
    expect_type(TreeSearch:::ts_sector_diag(edge, co, tp, wt, lv, nTaxa + 2L),
                "list")
    expect_type(TreeSearch:::ts_pool_test(list(edge, edge), c(0, 0), nTaxa),
                "list")

    # The benchmark harness reports no per-phase work rather than timing
    # zero-byte copies over empty buffers (agent-issues/TreeSearch#151).
    phases <- TreeSearch:::ts_bench_tbr_phases(edge, co, tp, wt, lv)
    expect_equal(phases$n_clips, 0L, info = nm)
    expect_equal(phases$n_snapshot_iters, 0L, info = nm)
  }
})


test_that("the user-facing API survives zero Fitch words", {
  # Three taxa have exactly one topology, so search is refused outright.
  expect_error(
    MaximizeParsimony(zeroWordSets$minimumTips,
                      tree = FixedTree(zeroWordSets$minimumTips, 1L),
                      maxReplicates = 2L, verbosity = 0L),
    "at least 4 taxa"
  )

  for (nm in setdiff(names(zeroWordSets), "minimumTips")) {
    dataset <- zeroWordSets[[nm]]
    tree <- FixedTree(dataset, 1L)
    for (concavity in list(Inf, 10, "profile")) {
      trees <- suppressWarnings(
        MaximizeParsimony(dataset, tree = tree, concavity = concavity,
                          maxReplicates = 2L, verbosity = 0L))
      expect_s3_class(trees[[1]], "phylo")
      expect_setequal(trees[[1]]$tip.label, names(dataset))
      expect_true(is.numeric(TreeLength(tree, dataset, concavity = concavity)))
    }
    # One length per character, not per distinct pattern.
    expect_length(CharacterLength(tree, dataset),
                  sum(attr(dataset, "weight")))
    expect_type(Consistency(dataset, tree), "double")
    if (sum(attr(dataset, "weight")) >= 3) {
      # Jackknifing a two-character matrix deletes nothing and is rejected
      # before any C++ runs, so it would test nothing here.
      expect_s3_class(suppressWarnings(Resample(dataset, tree,
                                                maxReplicates = 2L)),
                      "multiPhylo")
    }
  }
})


test_that("all-hierarchy data (zero Fitch words) survives HSJ and XFORM", {
  # Complements test-ts-hsj.R's search-level regression: HSJ and XFORM both
  # zero-weight every hierarchy character, so a dataset whose characters are
  # ALL hierarchical empties the equal-weights kernel while the search must
  # still work from the hierarchy term alone.
  mat <- matrix(c(
    "1", "0", "0", "-", "-",
    "1", "1", "1", "0", "1",
    "0", "-", "1", "1", "0",
    "1", "0", "1", "0", "0",
    "0", "-", "0", "-", "-",
    "1", "1", "0", "-", "-"
  ), nrow = 6, byrow = TRUE, dimnames = list(paste0("t", 1:6), NULL))
  dataset <- phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                              ambiguity = "?")
  hierarchy <- CharacterHierarchy("1" = 2L, "3" = 4:5)
  expect_length(setdiff(seq_len(5L), HierarchyChars(hierarchy)), 0L)
  tree <- FixedTree(dataset, 1L)

  for (mode in c("hsj", "xform")) {
    expect_true(is.numeric(TreeLength(tree, dataset, hierarchy = hierarchy,
                                      inapplicable = mode)))
    # The x-transformation's score is rooting-dependent, so on all-hierarchy
    # data its MPT set can hold trees of differing length at a common
    # rooting; the search is required to say so rather than report the set
    # as if it were homogeneous (see ?MaximizeParsimony).
    trees <- suppressWarnings(
      MaximizeParsimony(dataset, tree = tree, hierarchy = hierarchy,
                        inapplicable = mode, hsj_alpha = 1,
                        maxReplicates = 3L, targetHits = 2L, verbosity = 0L))
    expect_s3_class(trees[[1]], "phylo")
    expect_setequal(trees[[1]]$tip.label, names(dataset))
  }
})


test_that("L3b incremental edge sets survive degenerate data", {
  # `l3b_active` (ts_tbr.cpp) needs `tree.n_tip >= 150` unless
  # TS_L3B_INCREMENTAL forces it, so its six memcpy sites and its
  # `&edge_set_buf[db]` address formations are unreachable from any
  # ordinarily-sized test.  Force the path on a small tree instead of
  # paying for 150 tips.
  original <- Sys.getenv("TS_L3B_INCREMENTAL", unset = NA)
  on.exit({
    if (is.na(original)) {
      Sys.unsetenv("TS_L3B_INCREMENTAL")
    } else {
      Sys.setenv(TS_L3B_INCREMENTAL = original)
    }
  }, add = TRUE)
  Sys.setenv(TS_L3B_INCREMENTAL = "1")

  informative <- DegenerateDat(12L, list(
    c(Rp("0", 6), Rp("1", 6)), c(Rp("0", 5), Rp("1", 7)),
    c(Rp("1", 3), Rp("0", 9)), rep(c("0", "1"), 6)
  ))
  for (dataset in list(informative, zeroWordSets$allConstant,
                       zeroWordSets$allAutapomorphic)) {
    tree <- FixedTree(dataset, 1L)
    bits <- DatBits(dataset)
    result <- TreeSearch:::ts_tbr_search(tree$edge, bits$contrast,
                                         bits$tipData, bits$weight,
                                         bits$levels)
    expect_true(is.numeric(result$score))
    expect_type(TreeSearch:::ts_ratchet_search(tree$edge, bits$contrast,
                                               bits$tipData, bits$weight,
                                               bits$levels, nCycles = 2L),
                "list")
  }
})
