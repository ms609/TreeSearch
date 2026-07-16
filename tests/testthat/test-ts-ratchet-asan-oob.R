# DIAGNOSTIC repro for the ratchet multistate heap-corruption bug
# (memory node: ratchet-multistate-segfault). Deliberately NOT gated by
# skip_extended(), so the gcc-ASAN workflow (ASan.yml, which runs plain
# `tests`) executes it and halts on the first out-of-bounds access, giving
# an exact src/*.cpp file:line trace.
#
# Cumulative per-call heap corruption: a single isolated ratchet call is
# fine; the OOB write only lands somewhere fatal after many repeated calls
# in one process. Worse (fewer calls to crash) with more chars/states/tips.
# On native Windows/MinGW the process dies with exit 139 partway through
# the loop; under ASan the offending access is reported deterministically.
#
# EXPECTED to fail (ASan abort) until the underlying OOB is fixed. Once the
# bug is fixed this becomes a passing regression test.

# --- Config 1 (PRIMARY): ZERO_ONLY perturbation, 3-state EW -----------------
# 30 tips x 25 chars, states 0:2, equal weights, perturbMode = 0 (ZERO_ONLY).
# This is the config that crashes soonest in the documented repro.

test_that("ratchet ZERO_ONLY 3-state EW survives repeated calls (ASan)", {
  set.seed(6284)
  m <- matrix(sample(0:2, 30 * 25, replace = TRUE), nrow = 30,
              dimnames = list(paste0("t", 1:30), NULL))
  dataset <- MatrixToPhyDat(m)
  ds <- make_ts_data(dataset)
  tr <- as.phylo(1, 30)

  ok <- TRUE
  for (i in 1:30) {
    set.seed(i)
    r <- TreeSearch:::ts_ratchet_search(
      tr$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      nCycles = 5L, perturbMode = 0L, perturbProb = 0.15)
    ok <- ok && is.finite(r$score) && r$score >= 0
  }
  expect_true(ok)
})

# --- Config 2 (SECONDARY): UPWEIGHT_ONLY perturbation, IW (concavity 3) ------
# perturbMode = 1 with min_steps + concavity engages the IW pattern_freq
# upweight path. Only reached if Config 1 is clean (ASan halts on first hit).

test_that("ratchet UPWEIGHT_ONLY IW survives repeated calls (ASan)", {
  set.seed(6284)
  m <- matrix(sample(0:2, 30 * 25, replace = TRUE), nrow = 30,
              dimnames = list(paste0("t", 1:30), NULL))
  dataset <- MatrixToPhyDat(m)
  ds <- make_ts_data(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
  tr <- as.phylo(1, 30)

  ok <- TRUE
  for (i in 1:30) {
    set.seed(i)
    r <- TreeSearch:::ts_ratchet_search(
      tr$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      nCycles = 5L, perturbMode = 1L, perturbProb = 0.15,
      min_steps = minSteps, concavity = 3.0)
    ok <- ok && is.finite(r$score) && r$score >= 0
  }
  expect_true(ok)
})
