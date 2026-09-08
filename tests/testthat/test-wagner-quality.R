# Regression guard for the Wagner insertion-cost fix.
#
# A greedy Wagner (random-addition-sequence) tree should land only a few percent
# above the most-parsimonious length.  A historical bug scored candidate
# insertion edges with the UNION of the two endpoints' final-state sets
# (final[A] | final[D]); that set is too permissive, undercounts insertion cost,
# and degrades greedy placement to near-random -- producing Wagner trees ~+30%
# over the optimum.  The fix scores each edge with the exact DIRECTIONAL edge set
# combine(down[D], up[D]) (see compute_insertion_edge_sets / ts_fitch.cpp).
#
# This test pins the fixed behaviour: the mean of several addition trees must sit
# within a few percent of the known MPT length.  Fixed Wagner is ~+4%; the bug
# was ~+30%, so the threshold separates them with wide margin.

test_that("addition trees are within a few % of the MPT length (EW Fitch)", {
  skip_if_not_installed("TreeTools")
  data("inapplicable.phyData", package = "TreeSearch")

  # Equal-weights Fitch: treat inapplicable tokens as fully ambiguous ('?').
  Fitchify <- function(p) {
    m <- TreeTools::PhyDatToMatrix(p, ambigNA = FALSE)
    m[m == "-"] <- "?"
    TreeTools::MatrixToPhyDat(m)
  }

  # Known most-parsimonious lengths (EW Fitch), established by thorough search
  # (TNT 1.6 and TreeSearch agree on the ?-recoded matrices).
  cases <- list(
    Zanol2014 = 1261,
    Zhu2013   = 624
  )

  for (nm in names(cases)) {
    phy <- Fitchify(inapplicable.phyData[[nm]])
    mpt <- cases[[nm]]

    # Genuine random-addition-sequence Wagner trees (full random order).
    scores <- vapply(seq_len(8), function(s) {
      set.seed(s)
      TreeLength(AdditionTree(phy, sequence = sample(seq_along(phy))), phy)
    }, double(1))

    meanRatio <- mean(scores) / mpt
    # Fixed Wagner sits ~+4% (Zanol) to ~+6% (Zhu) over the MPT; the historical
    # union-of-finals bug produced ~+30%.  An 8% bound is the regression guard.
    expect_lt(meanRatio, 1.08)            # within a few % of the MPT
    expect_gte(min(scores), mpt)          # cannot beat the established optimum
  }
})
