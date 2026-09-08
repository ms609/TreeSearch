# Memory-safety stress guard for the C++ ratchet on multistate data.
# Deliberately NOT gated by skip_extended(), so the gcc-ASAN workflow
# (ASan.yml runs plain `tests`) executes it: many repeated ts_ratchet_search
# calls in one process on a 30-tip x 25-char 3-state matrix, which is exactly
# the shape that would surface a heap out-of-bounds write in the TBR/fitch
# undo/rescore machinery if one were ever introduced.
#
# Origin: a 2026-07-16 report of an exit-139 heap-corruption crash here
# (memory node: ratchet-multistate-segfault). Investigation could NOT
# reproduce it on ANY clean build — Linux gcc-ASAN clean, and a fresh local
# Windows -O2 build ran 500 iterations without a crash — while the search
# kernels were byte-identical to the reported commit. The crash was most
# consistent with a stale-object / mixed-ABI incremental build (see memory
# node stale-object-abi-gotcha), not a defect in committed source.
#
# So: this test is EXPECTED TO PASS. If it ever fails, first REBUILD CLEAN
# (rm src/*.o; CCACHE_DISABLE=1 R CMD INSTALL --preclean) and re-run before
# treating it as a genuine source regression.

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
