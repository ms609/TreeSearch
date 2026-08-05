# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Regression tests for three guard clauses at the HSJ/XFORM bridge
# (src/ts_rcpp.cpp: unpack_hsj(), unpack_xform()) and the collapse gate
# (src/ts_collapsed.cpp).
#
# T-398 (#14): unpack_hsj() enabled HSJ scoring (scoring_mode = HSJ) even when
# hsjTipLabels was present-but-NULL, leaving ds.tip_labels empty and
# segfaulting score_hierarchy_block(). Fixed with an explicit presence/NULL
# check that Rcpp::stop()s instead of silently skipping population. Review
# also surfaced a second entrance to the same crash class -- a non-NULL but
# too-narrow hsjTipLabels not covering every block's primary/secondary index
# -- closed with the same guard.
#
# T-397 (#13): unpack_xform() never validated cost_matrix dimensions,
# combo_grid row count, tip_sec_known dimensions, or tip_states range
# against n_states, unlike its sibling ts_sankoff_test() (guarded since
# 0856748f). A same-length-but-wrong-shape cost matrix reads garbage with no
# warning at all (Rcpp's Matrix::operator() only bounds-checks the linear
# offset, not (row, col)).
#
# T-408 (#22): the collapse guards in ts_collapsed.cpp keyed on
# ds.scoring_mode alone, disabling collapse for an HSJ/XFORM config with NO
# hierarchy data -- a case collapse is provably safe for. Fixed to gate on
# hierarchy-data presence (matching DataSet::topology_independent()).

library("TreeTools")

.HierarchyToBlocks <- TreeSearch:::.HierarchyToBlocks
.BuildTipLabels <- TreeSearch:::.BuildTipLabels
.HSJAbsentState <- TreeSearch:::.HSJAbsentState
.NonHierarchyWeights <- TreeSearch:::.NonHierarchyWeights
ts_driven_search <- TreeSearch:::ts_driven_search
ts_collapse_pool <- TreeSearch:::ts_collapse_pool

make_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}


# =========================================================================
# T-398 / #14: hsjTipLabels omitted (left at compat-wrapper default NULL)
# must error cleanly, not segfault.
#
# A segfault kills the test process outright, so it cannot be asserted from
# inside testthat. The pre-fix segfault was reproduced separately with a
# standalone Rscript (exit code 139) -- see the PR body for that output.
# =========================================================================
test_that("HSJ scoring with omitted hsjTipLabels errors instead of segfaulting", {
  mat <- matrix(c(
    "0", "0",
    "0", "0",
    "1", "0",
    "1", "1",
    "1", "1",
    "1", "-"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("2" = integer(0))
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)
  blocks <- .HierarchyToBlocks(h)

  # hsjTipLabels intentionally omitted -- exercises the compat wrapper's own
  # default (R/ts-driven-compat.R), which previously reached unpack_hsj()
  # with hsjConfig$hsjTipLabels == NULL.
  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      hierarchyBlocks = blocks,
      hsjAlpha = 1.0,
      hsjAbsentState = .HSJAbsentState(ds),
      maxReplicates = 1L
    ),
    "hsjTipLabels"
  )
})

test_that("HSJ scoring with too-narrow hsjTipLabels errors", {
  mat <- matrix(c(
    "0", "0",
    "0", "0",
    "1", "0",
    "1", "1",
    "1", "1",
    "1", "-"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("2" = integer(0))
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)
  blocks <- .HierarchyToBlocks(h)
  tl <- .BuildTipLabels(ds)

  # A non-NULL but too-narrow hsjTipLabels (missing the primary's column)
  # must be caught the same way an omitted one is: score_hierarchy_block()
  # indexes tip_labels at [t * n_orig_chars + block$primary], so dropping
  # the column covering that index reads past the vector.
  narrow_tl <- tl[, -ncol(tl), drop = FALSE]

  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      hierarchyBlocks = blocks,
      hsjTipLabels = narrow_tl,
      hsjAlpha = 1.0,
      hsjAbsentState = .HSJAbsentState(ds),
      maxReplicates = 1L
    ),
    "hsjTipLabels"
  )
})


# =========================================================================
# T-397 / #13: mis-shaped cost matrix through unpack_xform() must error
# with the expected/actual dimensions, matching ts_sankoff_test()'s style.
# Covers the same-length-but-wrong-shape case (1x9 for a 3x3) that Rcpp's
# own indexing does not warn about.
# =========================================================================
test_that("Xform bridge errors on mis-shaped cost matrix", {
  # Two informative secondary states (sec = 0, 1) -> n_states = 1(absent) +
  # 2(present combos) = 3, giving a 3x3 cost matrix.
  mat <- matrix(c(
    "0", "-",
    "1", "0",
    "1", "1",
    "1", "0"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  recoded <- RecodeHierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  expect_equal(dim(blk$cost_matrix), c(3, 3))

  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  # Same-length-but-wrong-shape: 1x9 carries the same 9 values as the 3x3
  # matrix, so Rcpp's linear-offset indexing reads it with no warning.
  bad_blk <- blk
  bad_blk$cost_matrix <- matrix(as.vector(blk$cost_matrix), nrow = 1, ncol = 9)

  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      xformChars = list(bad_blk),
      maxReplicates = 1L
    ),
    "cost_matrix has dimensions 1 x 9.*3 states"
  )
})

test_that("Xform bridge errors on undersized combo_grid", {
  mat <- matrix(c(
    "1", "0", "0",
    "1", "0", "0",
    "1", "1", "1",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  recoded <- RecodeHierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  bad_blk <- blk
  bad_blk$combo_grid <- blk$combo_grid[-1, , drop = FALSE]

  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      xformChars = list(bad_blk),
      maxReplicates = 1L
    ),
    "combo_grid has"
  )
})

test_that("Xform bridge errors on mis-shaped tip_sec_known", {
  mat <- matrix(c(
    "1", "0", "0",
    "1", "0", "0",
    "1", "1", "1",
    "1", "1", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2:3)

  recoded <- RecodeHierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  bad_blk <- blk
  bad_blk$tip_sec_known <- blk$tip_sec_known[-1, , drop = FALSE]

  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      xformChars = list(bad_blk),
      maxReplicates = 1L
    ),
    "tip_sec_known has"
  )
})

test_that("Xform bridge errors on out-of-range tip_states", {
  mat <- matrix(c(
    "0", "-", "0",
    "1", "0", "1",
    "1", "0", "0",
    "1", "0", "1"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(paste0("t", 1:4), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("1" = 2L)

  recoded <- RecodeHierarchy(ds, h)
  blk <- recoded$sankoff_chars[[1]]
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)

  bad_blk <- blk
  bad_blk$tip_states[1] <- blk$n_states + 5L  # out of [0, n_states)

  expect_error(
    ts_driven_search(
      contrast = at$contrast,
      tip_data = tip_data,
      weight = adj_w,
      levels = at$levels,
      xformChars = list(bad_blk),
      maxReplicates = 1L
    ),
    "out of range"
  )
})


# =========================================================================
# T-408 / #22: an HSJ config with empty hierarchyBlocks must collapse
# zero-length branches exactly as the no-hsjConfig (EW) case does -- the
# guard must not disable collapse when no hierarchy data is actually
# present, only when it is.
# =========================================================================
test_that("Collapse fires for an HSJ config with no hierarchy blocks", {
  # Reuses the T-330 reproducing configuration (char1 zero-length-supports
  # the (t4,t5) clade once char2's weight is zeroed) but with an HSJ config
  # whose hierarchy_blocks is empty -- scoring_mode == HSJ, yet no hierarchy
  # data exists for the collapse kernel to be blind to.
  mat <- matrix(c(
    "0", "0",
    "0", "0",
    "1", "0",
    "1", "1",
    "1", "1",
    "1", "-"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  ds <- make_dat(mat)
  h <- CharacterHierarchy("2" = integer(0))
  tr <- Preorder(RenumberTips(
    ape::read.tree(text = "(((t1,t2),(t3,(t4,t5))),t6);"), names(ds)))
  at <- attributes(ds)
  adj_w <- as.integer(.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)
  scoringConfig <- list(min_steps = integer(0), concavity = Inf,
                        xpiwe = FALSE, xpiwe_r = 0.5, xpiwe_max_f = 5.0,
                        obs_count = integer(0), infoAmounts = NULL)

  hsjConfig <- list(
    hierarchyBlocks = list(),
    hsjAlpha = 1.0,
    hsjTipLabels = matrix(integer(0), nrow = length(ds), ncol = 0),
    hsjAbsentState = 0L)

  n_in <- nrow(tr$edge)
  cp_ew <- ts_collapse_pool(
    list(tr$edge), at$contrast, tip_data, adj_w, at$levels,
    scoringConfig, NULL, NULL, NULL)
  cp_hsj <- ts_collapse_pool(
    list(tr$edge), at$contrast, tip_data, adj_w, at$levels,
    scoringConfig, hsjConfig, NULL, NULL)

  # Sanity: the (t4,t5) clade genuinely is collapsible under plain EW.
  expect_lt(nrow(cp_ew$trees[[1]]), n_in)
  # The empty-hierarchy HSJ config must collapse identically -- not be
  # blocked by the scoring_mode-only guard.
  expect_equal(nrow(cp_hsj$trees[[1]]), nrow(cp_ew$trees[[1]]))
})
