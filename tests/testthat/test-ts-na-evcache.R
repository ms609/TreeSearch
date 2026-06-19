# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Regression guard for the exact_verify_sweep optimum cache (src/ts_tbr.cpp).
#
# exact_verify_sweep is the NA convergence certifier: at every NA convergence it
# sweeps the full unrooted-TBR neighbourhood and either applies an improver or
# declares a genuine optimum.  Declaring an optimum is expensive, so FALSE
# ("genuine optimum") verdicts are memoized.  Such a verdict is valid ONLY under
# the weighting regime in force when it was recorded -- but the parsimony ratchet
# mutates the regime IN PLACE mid-search (active_mask / upweight_mask /
# pattern_freq; see save_perturb_state in ts_ratchet.cpp) and runs NA TBR under
# both base and perturbed weights within one cycle.  If the cache key omits the
# regime, a base-regime "optimal" verdict leaks into a perturbed pass and the
# search silently skips the improving moves the ratchet exists to find.
#
# This regression is INVISIBLE to the NA oracle (it never ratchets, so the regime
# is constant and the key collapses) and to final scores (always recomputed).
# So we test the key directly: ts_ev_cache_key_probe returns the EXACT key the
# cache uses (via the shared ts::exact_verify_cache_key helper) and the
# perturbation flags reproduce the three ways the ratchet changes the regime.
# None of the flags touch ds_fingerprint's inputs (n_tips/n_blocks/tip_states),
# so an observed key change is attributable to the weighting-regime term alone --
# this pins that weight_fingerprint is XORed INTO the composite key, not merely
# that the function exists.  Dropping `^ weight_fingerprint` from the key, or the
# known near-miss "key on upweight_mask alone" (a no-op for the DEFAULT ZERO_ONLY
# ratchet), fails the zero_active assertion below.

# Helpers from helper-ts.R: make_ts_data
make_na_ds <- function() {
  # 8 tips, 4 characters, several inapplicable ("-") tokens -> the NA path.
  mat <- matrix(
    c("-", "-", "-", "1", "1", "2", "2", "2",   # character 1
      "1", "1", "-", "-", "2", "2", "0", "0",   # character 2
      "0", "0", "1", "1", "-", "-", "2", "2",   # character 3
      "1", "2", "1", "2", "1", "2", "1", "2"),  # character 4
    nrow = 8, ncol = 4,
    dimnames = list(paste0("t", 1:8), NULL)
  )
  make_ts_data(MatrixToPhyDat(mat))
}

ev_key <- function(edge, ds, ...) {
  TreeSearch:::ts_ev_cache_key_probe(
    edge, ds$contrast, ds$tip_data, ds$weight, ds$levels, ...
  )
}


test_that("exact_verify cache key is sensitive to every weighting-regime field", {
  ds   <- make_na_ds()
  edge <- as.phylo(42, 8)$edge

  base <- ev_key(edge, ds)

  # Deterministic: identical inputs -> identical key (no spurious differences
  # that would make every cache lookup miss).
  expect_identical(base, ev_key(edge, ds))

  # PRIMARY guard.  Zeroing one active_mask bit IS the default (ZERO_ONLY)
  # ratchet perturbation.  Must change the key -- catches both "weight_fp dropped
  # from the key" and the near-miss "key on upweight_mask alone" (which would be
  # a silent no-op here, since ZERO_ONLY never touches upweight_mask).
  expect_false(identical(base, ev_key(edge, ds, zero_active = TRUE)))

  # upweight_mask (UPWEIGHT_ONLY / MIXED ratchet modes) must also be in the key.
  expect_false(identical(base, ev_key(edge, ds, set_upweight = TRUE)))

  # pattern_freq (IW ratchet upweighting) must also be in the key.
  expect_false(identical(base, ev_key(edge, ds, bump_pattern_freq = TRUE)))
})


test_that("exact_verify cache key separates topologies and datasets", {
  ds   <- make_na_ds()
  edge <- as.phylo(42, 8)$edge
  base <- ev_key(edge, ds)

  # Topology must be in the key: a different tree under the same data/regime must
  # not collide (else a true optimum for tree A suppresses the sweep on tree B).
  edge2 <- as.phylo(99, 8)$edge
  expect_false(identical(base, ev_key(edge2, ds)))

  # A genuine dataset switch (different tip states) must change the key -- this
  # is also the cache's clear-trigger (ds_fingerprint), so entries cannot carry
  # over between datasets.
  mat2 <- matrix(
    c("1", "1", "1", "2", "2", "2", "2", "2",
      "1", "1", "2", "2", "2", "2", "0", "0",
      "0", "0", "1", "1", "2", "2", "2", "2",
      "1", "2", "1", "2", "1", "2", "1", "2"),
    nrow = 8, ncol = 4, dimnames = list(paste0("t", 1:8), NULL)
  )
  ds2 <- make_ts_data(MatrixToPhyDat(mat2))
  expect_false(identical(base, ev_key(edge, ds2)))
})


test_that("TS_EV_AUDIT re-verifies cache hits without false alarms on a clean cache", {
  # Exercises the live tripwire end-to-end: with TS_EV_AUDIT set, a cache hit is
  # distrusted and the full sweep is re-run; on a correct (regime-keyed) cache it
  # must confirm the optimum and NOT abort.  Run TBR to convergence twice on the
  # same NA tree+data in one process: the 1st call populates the optimum cache,
  # the 2nd hits it and triggers the audit re-verification.
  skip_on_cran()
  ds   <- make_na_ds()
  edge <- as.phylo(42, 8)$edge

  old <- Sys.getenv("TS_EV_AUDIT", unset = NA)
  Sys.setenv(TS_EV_AUDIT = "1")
  on.exit({
    if (is.na(old)) Sys.unsetenv("TS_EV_AUDIT") else Sys.setenv(TS_EV_AUDIT = old)
  }, add = TRUE)

  r1 <- TreeSearch:::ts_tbr_diagnostics(
    edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    concavity = -1, unrooted = TRUE
  )
  # Second descent from the converged tree should re-reach (and cache-hit on) the
  # same optimum; the audit must verify it cleanly rather than error.
  r2 <- expect_error(
    TreeSearch:::ts_tbr_diagnostics(
      r1$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      concavity = -1, unrooted = TRUE
    ),
    NA
  )
  expect_equal(r2$score, r1$score)
})
