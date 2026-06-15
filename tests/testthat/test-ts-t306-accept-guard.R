# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# =========================================================================
# T-306 regression: HSJ/XFORM SPR/NNI accept-paths must include the
# hierarchy-DP / Sankoff contribution.
# =========================================================================
# Bug (T-306): the SPR accept path in tbr_search() (src/ts_tbr.cpp) and the
# accept path in nni_search() (src/ts_search.cpp) updated `best_score` with a
# Fitch-only (EW) or IW/profile incremental delta.  For HSJ and XFORM scoring,
# score_tree() additionally adds a topology-dependent hierarchy-DP (HSJ) or
# Sankoff (XFORM) term that the delta omits, so the search's internal
# accept/reject decisions tracked the wrong objective.  The fix gates the
# incremental fast path on the scoring mode (EW/IW/XPIWE/PROFILE) and falls
# back to a full score_tree() rescore for HSJ/XFORM.
#
# This bug is SILENT at the reported-score level: run_single_replicate always
# recomputes the final score via score_tree() before pooling, and the driven
# pipeline's TBR rerooting moves (which always full_rescore) reliably recover
# the true optimum even with the buggy accept path.  An empirical sweep
# confirmed a buggy build still reaches the optimum from every Fitch-optimal
# "trap" start tree.  A black-box MaximizeParsimony() test therefore cannot
# discriminate buggy from fixed (this matches the T-303 sibling finding).
#
# What these tests DO lock in is that the HSJ/XFORM accept path — now routed
# through full_rescore() — produces a correct, internally-consistent, and
# deterministic driven search: it reaches the true global optimum, every
# returned tree's independent score equals the reported score (no best_score
# desync), and identical seeds give identical results.  A regression that
# broke the full_rescore fallback (stale state, wrong score, crash, or
# nondeterminism) would fail here.

library("TreeTools")

make_t306_dat <- function(mat) {
  MatrixToPhyDat(mat)
}

# 7-tip dataset chosen so the full tree space (945 unrooted binaries) is
# brute-forceable and the hierarchy contribution sharpens the optimum:
# the HSJ/XFORM optimum is a STRICT subset of the Fitch optimum, so reaching
# it exercises the hierarchy-aware scoring, not just the Fitch component.
#   columns: primary  sec2  sec3  nh4  nh5  nh6  nh7  nh8
t306_mat <- matrix(unlist(strsplit(c(
  "0--00110",
  "0--01101",
  "0--10011",
  "10010110",
  "10101001",
  "11011100",
  "11110011"
), "")), nrow = 7, byrow = TRUE,
dimnames = list(paste0("t", 1:7), NULL))


test_that("MaximizeParsimony HSJ reaches the true HSJ optimum (T-306)", {
  ds <- make_t306_dat(t306_mat)
  h <- CharacterHierarchy("1" = 2:3)

  # Brute force the whole tree space with the authoritative scorer
  # (TreeLength uses the same path/units as the driven search's score_tree).
  all_trees <- as.phylo(0:944, nTip = 7, tipLabels = names(ds))
  hsj_sc <- TreeLength(all_trees, ds, hierarchy = h,
                       inapplicable = "hsj", hsj_alpha = 1.0)
  fit_sc <- TreeLength(all_trees, ds, concavity = Inf)
  opt <- min(hsj_sc)
  hsj_opt <- which(abs(hsj_sc - opt) < 1e-9)
  fit_opt <- which(fit_sc == min(fit_sc))

  # The dataset retains its discriminating structure: a sharp HSJ optimum
  # that is a strict subset of the (flatter) Fitch optimum.
  expect_lt(length(hsj_opt), 10L)
  expect_true(all(hsj_opt %in% fit_opt))
  expect_false(setequal(hsj_opt, fit_opt))

  for (s in c(1L, 7L, 42L, 256L)) {
    set.seed(s)
    res <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "hsj",
                             hsj_alpha = 1.0, maxReplicates = 4L,
                             targetHits = 1L, verbosity = 0L)
    expect_s3_class(res[[1]], "phylo")
    expect_equal(length(res[[1]]$tip.label), 7L)

    reported <- attr(res, "score")
    # The driven search finds the true global HSJ optimum.
    expect_equal(reported, opt)
    # No best_score desync: each returned tree's independent HSJ score
    # equals the reported optimum.
    recomputed <- TreeLength(res, ds, hierarchy = h,
                             inapplicable = "hsj", hsj_alpha = 1.0)
    expect_true(all(abs(recomputed - reported) < 1e-9))
  }

  # Determinism: identical seeds yield identical optima and tree counts.
  set.seed(99)
  a <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "hsj",
                         hsj_alpha = 1.0, maxReplicates = 4L,
                         targetHits = 1L, verbosity = 0L)
  set.seed(99)
  b <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "hsj",
                         hsj_alpha = 1.0, maxReplicates = 4L,
                         targetHits = 1L, verbosity = 0L)
  expect_equal(attr(a, "score"), attr(b, "score"))
  expect_equal(length(a), length(b))
})


test_that("MaximizeParsimony XFORM reaches the true Sankoff optimum (T-306)", {
  ds <- make_t306_dat(t306_mat)
  h <- CharacterHierarchy("1" = 2:3)

  all_trees <- as.phylo(0:944, nTip = 7, tipLabels = names(ds))
  xf_sc <- TreeLength(all_trees, ds, hierarchy = h, inapplicable = "xform")
  fit_sc <- TreeLength(all_trees, ds, concavity = Inf)
  opt <- min(xf_sc)
  xf_opt <- which(abs(xf_sc - opt) < 1e-9)

  # The Sankoff recoding is genuinely engaged (xform landscape differs from a
  # plain-Fitch landscape) and yields a non-trivial, reasonably sharp optimum.
  expect_false(isTRUE(all.equal(xf_sc, fit_sc)))
  expect_gt(length(unique(round(xf_sc, 6))), 1L)
  expect_lt(length(xf_opt), 30L)

  for (s in c(1L, 7L, 42L, 256L)) {
    set.seed(s)
    res <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "xform",
                             maxReplicates = 4L, targetHits = 1L,
                             verbosity = 0L)
    expect_s3_class(res[[1]], "phylo")
    expect_equal(length(res[[1]]$tip.label), 7L)

    reported <- attr(res, "score")
    expect_equal(reported, opt)
    recomputed <- TreeLength(res, ds, hierarchy = h, inapplicable = "xform")
    expect_true(all(abs(recomputed - reported) < 1e-9))
  }

  set.seed(99)
  a <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "xform",
                         maxReplicates = 4L, targetHits = 1L, verbosity = 0L)
  set.seed(99)
  b <- MaximizeParsimony(ds, hierarchy = h, inapplicable = "xform",
                         maxReplicates = 4L, targetHits = 1L, verbosity = 0L)
  expect_equal(attr(a, "score"), attr(b, "score"))
  expect_equal(length(a), length(b))
})
