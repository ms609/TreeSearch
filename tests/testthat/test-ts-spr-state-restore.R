# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Targeted test for SPR state restoration after rejected regrafts.
#
# After spr_regraft + full_rescore + rejection + spr_unregraft + spr_unclip,
# scoring arrays may be stale for nodes on the regraft-to-root path that
# aren't on the clip-to-root path. This test verifies that the final output
# of spr_search always matches an independent full rescore, even when many
# rejections occur (which maximizes stale-array exposure).

ts_spr <- function(tree, ds, maxHits = 20L, concavity = Inf) {
  TreeSearch:::ts_spr_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, concavity = concavity)
}

test_that("SPR final score matches independent rescore (EW, many rejections)", {
  # Use a 20-tip tree with weak signal to maximize rejected regrafts:
  # most indirect evaluations will look promising but full rescore rejects.
  set.seed(5872)
  mat <- matrix(sample(0:3, 20 * 4, replace = TRUE),
                nrow = 20, dimnames = list(paste0("t", 1:20), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (start in c(1, 42, 200, 500, 1000)) {
    tree <- as.phylo(start, 20)
    set.seed(3719 + start)
    result <- ts_spr(tree, ds, maxHits = 5L)

    # Build tree from result and independently rescore
    rt <- tree
    rt$edge <- result$edge
    independent_score <- ts_score(rt, ds)
    expect_equal(result$score, independent_score,
                 info = paste("start =", start))
    validate_result(result, 20L)
  }
})

test_that("SPR final score matches independent rescore (IW, many rejections)", {
  set.seed(6413)
  mat <- matrix(sample(0:2, 18 * 5, replace = TRUE),
                nrow = 18, dimnames = list(paste0("t", 1:18), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)

  for (start in c(1, 50, 300)) {
    tree <- as.phylo(start, 18)
    set.seed(2847 + start)
    result <- ts_spr(tree, ds, maxHits = 5L, concavity = 10.0)

    rt <- tree
    rt$edge <- result$edge
    independent_score <- ts_score(rt, ds, concavity = 10.0)
    expect_equal(result$score, independent_score, tolerance = 1e-10,
                 info = paste("IW start =", start))
    validate_result(result, 18L)
  }
})

test_that("SPR final score matches independent rescore (NA dataset)", {
  skip_if_not_installed("TreeSearch")
  data(inapplicable.phyData, package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  n_tip <- length(dataset)

  for (start in c(1, 42, 100)) {
    tree <- as.phylo(start, n_tip)
    set.seed(9054 + start)
    result <- ts_spr(tree, ds, maxHits = 3L)

    rt <- tree
    rt$edge <- result$edge
    independent_score <- ts_score(rt, ds)
    expect_equal(result$score, independent_score,
                 info = paste("NA start =", start))
    validate_result(result, n_tip)
  }
})
