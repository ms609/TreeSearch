# Regression guard for the ratchet floor bug (ratchetCycles = 0 was not a no-op).
#
# Three stacked floors in ts_driven.cpp used to force >= 1 ratchet cycle even
# when the user asked for 0:
#   1. ratchet_per = max(1, ...)            -- ceiling division had no zero-case
#   2. the ratchet block was called unconditionally (and ratchet_search() runs
#      an initial TBR pass before its cycle loop, so n_cycles = 0 still perturbs)
#   3. the adaptive_level re-floor  max(1, base_ratchet_cycles * scale)  silently
#      resurrected a disabled ratchet on the *default* (adaptive) strategy.
# With all three guarded, ratchetCycles = 0 skips the ratchet phase entirely and
# ratchet_ms is exactly 0.  These run on the serial path (nThreads = 1L), where
# the guards live; the parallel path is tracked separately.

data("inapplicable.phyData", package = "TreeSearch")
ds <- inapplicable.phyData[["Vinther2008"]]

test_that("ratchetCycles = 0 disables ratchet under the default (adaptive) strategy", {
  set.seed(108)
  result <- MaximizeParsimony(ds, ratchetCycles = 0L, maxReplicates = 3L,
                              targetHits = 1L, nThreads = 1L, verbosity = 0L)
  timings <- attr(result, "timings")
  # The ratchet phase never ran, so no time is attributed to it.  (Pre-fix this
  # was strongly positive because the floors forced >= 1 cycle every replicate.)
  expect_identical(unname(timings[["ratchet_ms"]]), 0)
  expect_s3_class(result, "multiPhylo")
})

test_that("ratchet still runs when ratchetCycles > 0", {
  # Complement: the zero-guards must not disable an enabled ratchet.
  set.seed(108)
  result <- MaximizeParsimony(ds, ratchetCycles = 12L, maxReplicates = 3L,
                              targetHits = 1L, nThreads = 1L, verbosity = 0L)
  timings <- attr(result, "timings")
  expect_gt(timings[["ratchet_ms"]], 0)
})
