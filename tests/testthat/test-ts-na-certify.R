# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Guards for TBRParams::certify_unrooted -- the gate on exact_verify_sweep, the
# NA convergence certifier that is 97.7% of tbr_search wall when it fires
# (dev/profiling/na-exact-verify-dominates.md).
#
# Three things need pinning, and none of them is visible in a final score alone:
#
# 1. The gate is OPT-IN (TS_NA_NOCERTIFY).  With the variable unset every caller
#    certifies exactly as before, whatever `certify_unrooted` says.  If that
#    inverts, a quality-sensitive change ships silently.
# 2. On the skip path, `best_score` must still describe the RETURNED tree.
#    exact_verify_sweep re-synced it on both exits (full_rescore at entry, and
#    again before returning false); skipping the sweep removes that sync.  A
#    score that drifts from its topology is the ts_tbr.cpp:660 bug class, and in
#    a floor-attainment panel it would read as a quality regression that is
#    really a reporting bug.
# 3. `naDiag$n_evs_skipped` must actually move.  do_reroot -- the gate on
#    exact_verify_sweep -- requires tabu_size == 0, and the shipped presets set
#    tabuSize = 100/200, so most call sites never reach the certifier at all.
#    Without this counter, "the flag changed nothing" is indistinguishable from
#    "the flag never fired", which is exactly the trap the panel must avoid.

naDataset <- function() {
  # 12 tips, 8 characters, inapplicable ("-") tokens in most -> the NA path.
  mat <- matrix(
    c("-", "-", "-", "-", "1", "1", "1", "2", "2", "2", "0", "0",
      "1", "1", "-", "-", "-", "2", "2", "0", "0", "1", "1", "2",
      "0", "0", "1", "1", "-", "-", "-", "2", "2", "0", "1", "1",
      "1", "2", "1", "2", "1", "2", "0", "1", "-", "-", "-", "0",
      "2", "2", "0", "0", "1", "1", "2", "-", "-", "1", "0", "1",
      "-", "1", "1", "2", "2", "0", "0", "1", "1", "2", "-", "-",
      "0", "1", "2", "0", "1", "2", "-", "-", "0", "1", "2", "0",
      "1", "0", "0", "1", "-", "1", "1", "0", "2", "2", "1", "0"),
    nrow = 12, ncol = 8,
    dimnames = list(paste0("t", 1:12), NULL)
  )
  MatrixToPhyDat(mat)
}

# Recode "-" -> "?" : identical pattern structure, but has_na is FALSE, so the
# certifier is structurally unreachable.  The control arm.
noNaDataset <- function() {
  mat <- matrix(
    c("-", "-", "-", "-", "1", "1", "1", "2", "2", "2", "0", "0",
      "1", "1", "-", "-", "-", "2", "2", "0", "0", "1", "1", "2",
      "0", "0", "1", "1", "-", "-", "-", "2", "2", "0", "1", "1",
      "1", "2", "1", "2", "1", "2", "0", "1", "-", "-", "-", "0",
      "2", "2", "0", "0", "1", "1", "2", "-", "-", "1", "0", "1",
      "-", "1", "1", "2", "2", "0", "0", "1", "1", "2", "-", "-",
      "0", "1", "2", "0", "1", "2", "-", "-", "0", "1", "2", "0",
      "1", "0", "0", "1", "-", "1", "1", "0", "2", "2", "1", "0"),
    nrow = 12, ncol = 8,
    dimnames = list(paste0("t", 1:12), NULL)
  )
  mat[mat == "-"] <- "?"
  MatrixToPhyDat(mat)
}

# tabuSize = 0 is REQUIRED for any of this to be observable: do_reroot gates on
# it, and the certifier hangs off do_reroot.  A run at the shipped tabuSize = 100
# would report n_evs == 0 in every arm and prove nothing.
searchNa <- function(dat, noCertify, seed = 1L, ...) {
  if (noCertify) {
    withr::local_envvar(c(TS_NA_NOCERTIFY = "1"))
  } else {
    withr::local_envvar(c(TS_NA_NOCERTIFY = NA))
  }
  set.seed(seed)
  MaximizeParsimony(dat, .rung = "default", verbosity = 0L,
                    nThreads = 1L, maxReplicates = 2L, tabuSize = 0L, ...)
}

test_that("the certification gate is opt-in via TS_NA_NOCERTIFY", {
  dat <- naDataset()
  on <- searchNa(dat, noCertify = FALSE)
  # Certifier reached, nothing skipped: the shipped default certifies.
  expect_gt(attr(on, "naDiag")$n_evs, 0)
  expect_equal(attr(on, "naDiag")$n_evs_skipped, 0)
})

test_that("clearing certify_unrooted skips certifications", {
  dat <- naDataset()
  off <- searchNa(dat, noCertify = TRUE)
  # The gate reached live call sites, so a wall/reach result from it is
  # attributable rather than vacuous.
  expect_gt(attr(off, "naDiag")$n_evs_skipped, 0)
})

test_that("the skipped-certification score describes the returned tree", {
  dat <- naDataset()
  off <- searchNa(dat, noCertify = TRUE)
  reported <- attr(off, "score")
  # Score by LABEL via TreeLength -- never edge + tip_data by index, which
  # RenumberTips permutes (na-validation-alignment-gotcha).
  rescored <- vapply(off, TreeLength, double(1), dataset = dat)
  expect_equal(min(rescored), reported)
  expect_true(all(abs(rescored - reported) < 1e-8))
})

test_that("certification is inert without inapplicables", {
  dat <- noNaDataset()
  for (noCertify in c(FALSE, TRUE)) {
    res <- searchNa(dat, noCertify = noCertify)
    # EW/IW indirect scans are exact, so has_na is FALSE and neither branch of
    # the gate is reachable.  Guards against the flag leaking into try_root_
    # edge_moves, which every non-NA search depends on.
    expect_equal(attr(res, "naDiag")$n_evs, 0)
    expect_equal(attr(res, "naDiag")$n_evs_skipped, 0)
  }
})

test_that("skipping certification does not change reachable scores here", {
  # A 12-tip matrix is small enough that both arms should find the same
  # optimum; this is a smoke check on the mechanism, NOT the floor-attainment
  # gate (which needs real matrices, several seeds, and per-matrix aggregation
  # -- see dev/profiling/na-certify-gate.md).
  dat <- naDataset()
  on <- searchNa(dat, noCertify = FALSE)
  off <- searchNa(dat, noCertify = TRUE)
  expect_equal(attr(off, "score"), attr(on, "score"))
})
