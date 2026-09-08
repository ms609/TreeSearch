# Tier 3: opt-in (TREESEARCH_EXTENDED_TESTS=true); see tests/testing-strategy.md
skip_on_cran()

# Regression guard for the exact_verify_sweep root-edge completeness fix.
#
# exact_verify_sweep certifies NA convergence by sweeping the TBR neighbourhood,
# but its clip loop structurally skips root-child clips, so the one unrooted edge
# the display root sits on (cL-cR) went unchecked.  On poor NA starts that left a
# root-edge improver undetected: the kernel declared convergence one or more steps
# above the true unrooted-TBR optimum.  The fix enumerates that edge exactly (via
# try_root_edge_moves_rescore) before declaring an optimum.
#
# Canonical documented failure (dev/benchmarks/tbr_oracle_na.R, Zanol2014):
# start #14 converged at 1323 with a 1320 neighbour (a 3-step miss).  After the
# fix it reaches 1320.  We assert the kernel reaches the true optimum, which a
# root-edge regression (back to 1323) would fail.  Tier 3: a 74-tip search to
# convergence is ~10s.

test_that("NA unrooted-TBR reaches the root-edge optimum on Zanol2014 start #14", {
  skip_extended()
  library(TreeSearch)
  library(TreeTools)

  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Zanol2014"]]
  at <- attributes(dataset)
  d <- list(contrast = at$contrast,
            tip_data = matrix(unlist(dataset, use.names = FALSE),
                              nrow = length(dataset), byrow = TRUE),
            weight = at$weight, levels = at$levels, labels = names(dataset))

  # Reproduce oracle start #14 exactly: same start tree and search seed.
  set.seed(7000L + 14L)
  start <- RandomTree(dataset, root = TRUE)
  edge <- Preorder(RenumberTips(start, d$labels))[["edge"]]
  set.seed(14L)
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L,
    concavity = -1, unrooted = TRUE
  )

  # Pre-fix this converged at 1323 (root-edge improver to 1320 missed).
  expect_lte(res$score, 1320)

  # The reported score must match an independent rescore of the returned tree
  # (guards against a score/tree mismatch in the root-edge apply path).
  conv <- structure(list(edge = res$edge, Nnode = length(d$labels) - 1L,
                         tip.label = d$labels), class = "phylo")
  rescored <- TreeSearch:::ts_fitch_score(
    Preorder(RenumberTips(conv, d$labels))[["edge"]],
    d$contrast, d$tip_data, d$weight, d$levels, concavity = -1
  )
  expect_equal(rescored, res$score)
})
