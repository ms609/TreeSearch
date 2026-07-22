# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Regression for red-team T-330 (P1).
#
# The zero-length-branch collapse (compute_collapsed_flags[_aggressive] in
# src/ts_collapsed.cpp) decides min-length-0 branches by iterating ONLY the
# standard ds.blocks[].  HSJ stores a clade's support in ds.hierarchy_blocks and
# XFORM in ds.sankoff_* — both zero-weighted out of the standard blocks (by
# .NonHierarchyWeights) and dropped by simplify_patterns.  So a clade supported
# solely by a hierarchy/Sankoff character looked unsupported and was contracted:
#   1. final collapse via ts_collapse_pool over-collapsed the result trees;
#   2. MPT-enumeration dedup wrongly merged HSJ/XFORM-distinct topologies.
# The fix no-ops the collapse flags (all-zero == nothing collapses) whenever
# scoring_mode is HSJ or XFORM, in the kernel so every call site is covered.

library("TreeTools")

# A 6-tip matrix mixing ONE standard character with an HSJ hierarchy primary:
#   char 1 (standard):    0,0,1,1,1,1  -> constant (state 1) across {t3,t4,t5,t6}
#   char 2 (HSJ primary): 0,0,0,1,1,-  -> present in {t4,t5}, absent in t3
# On (((t1,t2),(t3,(t4,t5))),t6) the (t4,t5) clade has ZERO char-1 length (all
# state 1) yet IS supported by the HSJ primary.  A standard character must exist
# (else every char is zero-weighted, total_words == 0, and collapse no-ops for a
# different reason); this is the finder's reproducing configuration.
t330_dat <- function() {
  mat <- matrix(c(
    "0", "0",
    "0", "0",
    "1", "0",
    "1", "1",
    "1", "1",
    "1", "-"
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  phangorn::phyDat(mat, type = "USER", levels = c("-", "0", "1"),
                   ambiguity = "?")
}

t330_start <- function(ds) {
  Preorder(RenumberTips(
    ape::read.tree(text = "(((t1,t2),(t3,(t4,t5))),t6);"), names(ds)))
}


# =========================================================================
# Manifestation #1, deterministic: ts_collapse_pool must not contract an
# HSJ-supported zero-standard-length clade.  Driving the kernel directly (no
# search) pins the behaviour without RNG.  The WITHOUT-hsjConfig (EW) arm is the
# control: it proves the clade genuinely has zero *standard* support (so the
# guard is what saves it) AND that the guard did not disable EW collapse.
# =========================================================================
test_that("ts_collapse_pool keeps HSJ-supported clade; EW still collapses", {
  ds <- t330_dat()
  h <- CharacterHierarchy("2" = integer(0))
  tr <- t330_start(ds)
  at <- attributes(ds)

  adj_w <- as.integer(TreeSearch:::.NonHierarchyWeights(ds, h))
  tip_data <- matrix(unlist(ds, use.names = FALSE),
                     nrow = length(ds), byrow = TRUE)
  scoringConfig <- list(min_steps = integer(0), concavity = Inf,
                        xpiwe = FALSE, xpiwe_r = 0.5, xpiwe_max_f = 5.0,
                        obs_count = integer(0), infoAmounts = NULL)
  hsjConfig <- list(
    hierarchyBlocks = TreeSearch:::.HierarchyToBlocks(h),
    hsjAlpha = 1.0,
    hsjTipLabels = TreeSearch:::.BuildTipLabels(ds),
    hsjAbsentState = TreeSearch:::.HSJAbsentState(ds))

  n_in <- nrow(tr$edge)
  cp_ew <- TreeSearch:::ts_collapse_pool(
    list(tr$edge), at$contrast, tip_data, adj_w, at$levels,
    scoringConfig, NULL, NULL, NULL)
  cp_hsj <- TreeSearch:::ts_collapse_pool(
    list(tr$edge), at$contrast, tip_data, adj_w, at$levels,
    scoringConfig, hsjConfig, NULL, NULL)

  # EW (guard not engaged): the zero-standard-length clade is contracted.
  expect_lt(nrow(cp_ew$trees[[1]]), n_in)
  # HSJ (guard engaged): nothing is contracted — the clade survives intact.
  expect_equal(nrow(cp_hsj$trees[[1]]), n_in)
})


# =========================================================================
# Manifestations #1 + #2 end-to-end (HSJ): collapse=TRUE must preserve both the
# resolution and the MPT count that collapse=FALSE reports.
# =========================================================================
test_that("HSJ MaximizeParsimony collapse=TRUE preserves resolution and MPTs", {
  ds <- t330_dat()
  h <- CharacterHierarchy("2" = integer(0))
  tr <- t330_start(ds)

  set.seed(1)
  res_t <- MaximizeParsimony(ds, tree = tr, inapplicable = "hsj", hierarchy = h,
                             maxReplicates = 3L, verbosity = 0L, collapse = TRUE)
  set.seed(1)
  res_f <- MaximizeParsimony(ds, tree = tr, inapplicable = "hsj", hierarchy = h,
                             maxReplicates = 3L, verbosity = 0L, collapse = FALSE)

  # collapse is post-hoc and, for HSJ, must be a no-op: as resolved as
  # collapse = FALSE (pre-fix it dropped Nnode 5 -> 3 by contracting the clade).
  expect_equal(res_t[[1]]$Nnode, res_f[[1]]$Nnode)
  expect_equal(res_t[[1]]$Nnode, NTip(res_t[[1]]) - 1L)
  # collapse must not change the MPT count (post-hoc no-op under HSJ).
  expect_equal(attr(res_t, "n_topologies"), attr(res_f, "n_topologies"))
  # Manifestation #2: the collect_pool dedup path (ts_tbr.cpp) had no
  # scoring-mode guard, so it merged HSJ-distinct MPTs into one during search.
  # Post-fix the genuinely distinct HSJ optima (this dataset has several) are
  # retained; pre-fix the count collapsed to 1.
  expect_gte(attr(res_f, "n_topologies"), 2L)
})


# =========================================================================
# XFORM end-to-end: a Sankoff-recoded hierarchy character must likewise survive.
# =========================================================================
test_that("XFORM MaximizeParsimony collapse=TRUE preserves resolution", {
  ds <- t330_dat()
  h <- CharacterHierarchy("2" = integer(0))
  tr <- t330_start(ds)

  set.seed(1)
  res_t <- suppressWarnings(MaximizeParsimony(
    ds, tree = tr, inapplicable = "xform", hierarchy = h,
    maxReplicates = 3L, verbosity = 0L, collapse = TRUE))
  set.seed(1)
  res_f <- suppressWarnings(MaximizeParsimony(
    ds, tree = tr, inapplicable = "xform", hierarchy = h,
    maxReplicates = 3L, verbosity = 0L, collapse = FALSE))

  expect_equal(res_t[[1]]$Nnode, res_f[[1]]$Nnode)
  expect_equal(attr(res_t, "n_topologies"), attr(res_f, "n_topologies"))
})
