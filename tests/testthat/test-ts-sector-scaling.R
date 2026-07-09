# Tier 2: skipped on CRAN; see tests/testing-strategy.md
# Exercises the TNT-style sector-size scaling knobs added to rss_search():
#   (a) nTip-scaled cap  -- env TS_SECT_MAXFRAC / TS_SECT_THRESHOLD / TS_SECT_MAXSIZE
#   (b) adaptive growth  -- env TS_SECT_GROW(+_INC/_SELFACT/_MOVEON/_START)
# Both are default-OFF and additive: with no env override the search must be
# byte-identical to the pre-feature path.  These tests guard that additivity
# and the safety of the opt-in paths (valid, non-worse trees); the anytime /
# reach benefit is measured on Hamilton, not here.
skip_on_cran()
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_rss <- function(tree, ds, minSize = 6L, maxSize = 50L, rssPicks = 0L,
                   ratchetCycles = 2L, maxHits = 1L) {
  TreeSearch:::ts_rss_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             minSectorSize = minSize,
                             maxSectorSize = maxSize,
                             acceptEqual = FALSE,
                             rssPicks = rssPicks,
                             ratchetCycles = ratchetCycles,
                             maxHits = maxHits)
}

# Set env vars for the duration of `code`, restoring prior state (incl. unset).
with_env <- function(vars, code) {
  nms <- names(vars)
  old <- Sys.getenv(nms, unset = NA)
  do.call(Sys.setenv, as.list(vars))
  on.exit(for (i in seq_along(nms)) {
    if (is.na(old[[i]])) Sys.unsetenv(nms[[i]])
    else do.call(Sys.setenv, stats::setNames(list(old[[i]]), nms[[i]]))
  })
  force(code)
}

# 40-tip structured dataset: large enough that internal clades span a range of
# sizes >= sectorMinSize, so both the size window and the growth ramp engage.
set.seed(913)
scal_mat <- matrix(sample(0:1, 40 * 30, replace = TRUE), nrow = 40,
                   dimnames = list(paste0("t", 1:40), NULL))
scal_ds <- make_ts_data(MatrixToPhyDat(scal_mat))
scal_tree <- as.phylo(1, 40)
scal_start <- ts_score(scal_tree, scal_ds)


test_that("scaling knobs are OFF by default (deterministic baseline)", {
  set.seed(1); a <- ts_rss(scal_tree, scal_ds)
  set.seed(1); b <- ts_rss(scal_tree, scal_ds)
  expect_equal(a$score, b$score)
  expect_equal(a$edge, b$edge)
  validate_result(a, 40L)
  expect_lte(a$score, scal_start)
})

test_that("TS_SECT_GROW=0 kill-switch is identical to default-off", {
  set.seed(1); base <- ts_rss(scal_tree, scal_ds)
  # Growth knobs are supplied but the kill-switch must neutralise them entirely.
  set.seed(1)
  off <- with_env(c(TS_SECT_GROW = "0", TS_SECT_GROW_INC = "50",
                    TS_SECT_GROW_SELFACT = "40", TS_SECT_GROW_START = "6"),
                  ts_rss(scal_tree, scal_ds))
  expect_equal(off$score, base$score)
  expect_equal(off$edge, base$edge)
})

test_that("frac=0 / no-op env leaves the search byte-identical", {
  set.seed(1); base <- ts_rss(scal_tree, scal_ds)
  set.seed(1)
  same <- with_env(c(TS_SECT_MAXFRAC = "0"),   # 0 = disabled -> baseline
                   ts_rss(scal_tree, scal_ds))
  expect_equal(same$score, base$score)
  expect_equal(same$edge, base$edge)
})

test_that("adaptive growth engages (more sector picks), valid and non-worse", {
  set.seed(1); base <- ts_rss(scal_tree, scal_ds)
  set.seed(1)
  grown <- with_env(c(TS_SECT_GROW = "1", TS_SECT_GROW_INC = "50",
                      TS_SECT_GROW_SELFACT = "40", TS_SECT_GROW_START = "6",
                      TS_SECT_GROW_MOVEON = "0"),
                    ts_rss(scal_tree, scal_ds))
  validate_result(grown, 40L)
  expect_lte(grown$score, scal_start)
  # Growth ramps many bands of selections up to the cap, well past the flat
  # ~2*nTip/avgSize auto budget -- a no-op implementation would tie here.
  expect_gt(grown$n_sectors_searched, base$n_sectors_searched)
})

test_that("growth `moveon` bounds the run but does NOT suppress the size ramp", {
  set.seed(1); base <- ts_rss(scal_tree, scal_ds)
  set.seed(1)
  mo <- with_env(c(TS_SECT_GROW = "1", TS_SECT_GROW_INC = "50",
                   TS_SECT_GROW_START = "6", TS_SECT_GROW_MOVEON = "2"),
                 ts_rss(scal_tree, scal_ds))
  validate_result(mo, 40L)
  expect_lte(mo$score, scal_start)
  # Regression guard: a small moveon must not preempt the ramp.  Growth escalates
  # to the ceiling BEFORE moveon can fire, so far more sectors are searched than
  # the flat auto budget.  (Pre-fix, a global moveon tripped in the first band
  # and left growth inert -- this would tie with baseline.)
  expect_gt(mo$n_sectors_searched, base$n_sectors_searched)
})

test_that("nTip-scaled cap ENLARGES sector selection (upward binding)", {
  # The load-bearing path: the fraction must be able to raise the effective cap
  # ABOVE the preset ceiling so larger clades enter `eligible`.  Baseline caps
  # sectors at 12 (eligible clades in [6, 12]); raising the ceiling to 200 lets
  # frac set eff_max = min(200, ceil(0.9*40)=36) = 36, so clades of size 13..36
  # now also become selectable -> the search must behave differently.
  set.seed(1); base <- ts_rss(scal_tree, scal_ds, maxSize = 12L)
  set.seed(1)
  up <- with_env(c(TS_SECT_MAXSIZE = "200", TS_SECT_MAXFRAC = "0.9",
                   TS_SECT_THRESHOLD = "10"),
                 ts_rss(scal_tree, scal_ds, maxSize = 12L))
  validate_result(up, 40L)
  expect_lte(up$score, scal_start)
  expect_false(isTRUE(all.equal(
    c(up$score, up$n_sectors_searched, as.integer(up$edge)),
    c(base$score, base$n_sectors_searched, as.integer(base$edge)))))
})

test_that("nTip-scaled cap RESTRICTS when the fraction is small", {
  # threshold 10 < 40 tips, frac 0.5 -> eff_max = min(50, 20) = 20, narrowing
  # the window below the preset ceiling.  Must stay valid and non-worse.
  set.seed(1)
  capped <- with_env(c(TS_SECT_MAXFRAC = "0.5", TS_SECT_THRESHOLD = "10"),
                     ts_rss(scal_tree, scal_ds))
  validate_result(capped, 40L)
  expect_lte(capped$score, scal_start)
})

test_that("threshold gates the fraction: below threshold == baseline", {
  set.seed(1); base <- ts_rss(scal_tree, scal_ds)
  # threshold 999 > 40 tips -> fraction never applies -> baseline behaviour,
  # even with a raised MAXSIZE (cap_base only matters once frac binds; the flat
  # path still filters on the preset maxSize=50 via eff_max=cap_base... so here
  # we keep MAXSIZE at the preset to assert pure gating).
  set.seed(1)
  gated <- with_env(c(TS_SECT_MAXFRAC = "0.5", TS_SECT_THRESHOLD = "999"),
                    ts_rss(scal_tree, scal_ds))
  expect_equal(gated$score, base$score)
  expect_equal(gated$edge, base$edge)
})
