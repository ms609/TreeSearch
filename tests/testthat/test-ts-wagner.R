# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Helper to prepare phyDat for C++ engine
prep_pd <- function(pd) {
  list(
    contrast = attr(pd, "contrast"),
    tip_data = t(vapply(pd, I, pd[[1]])),
    weight = attr(pd, "weight"),
    levels = attr(pd, "levels")
  )
}

test_that("ts_wagner_tree produces valid tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_true(is.list(result))
  expect_equal(ncol(result$edge), 2L)
  expect_equal(nrow(result$edge), 2L * (n_tip - 1L))
  expect_true(result$score > 0)
})

test_that("Wagner score matches ts_fitch_score", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(result$score, fitch_check)
})

test_that("Wagner score matches TreeLength", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  result_tree <- structure(list(
    edge = result$edge,
    tip.label = names(pd),
    Nnode = n_tip - 1L
  ), class = "phylo")

  expect_equal(result$score, TreeLength(result_tree, pd))
})

test_that("All tips present exactly once", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  child_nodes <- result$edge[, 2]
  tips_found <- sort(child_nodes[child_nodes <= n_tip])
  expect_equal(tips_found, seq_len(n_tip))
})

test_that("Same addition order gives same tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  order_seq <- seq_len(n_tip)
  r1 <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels, order_seq
  )
  r2 <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels, order_seq
  )

  expect_identical(r1$edge, r2$edge)
  expect_identical(r1$score, r2$score)
})

test_that("Random Wagner trees vary", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(7263)
  r1 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  set.seed(1984)
  r2 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_false(identical(r1$edge, r2$edge))
})

test_that("Random Wagner score verified", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(5511)
  result <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(result$score, fitch_check)
})

test_that("Small tree (5 tips) is correct", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd5 <- congreveLamsdellMatrices[[1]][1:5]
  d <- prep_pd(pd5)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(nrow(result$edge), 8L)
  expect_true(result$score > 0)

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)
})

test_that("Medium tree (20 tips) completes without error", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd20 <- congreveLamsdellMatrices[[1]][1:20]
  d <- prep_pd(pd20)

  expect_no_error({
    result <- TreeSearch:::ts_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )
  })

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)
})

test_that("Multiple datasets produce verified scores", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  for (i in 1:3) {
    pd <- congreveLamsdellMatrices[[i]]
    d <- prep_pd(pd)

    set.seed(3000 + i)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check,
                 info = paste("Dataset", i))
  }
})

# --- New tests for incremental scoring correctness ---

test_that("Wagner on inapplicable datasets matches fitch_score", {
  data("inapplicable.phyData", package = "TreeSearch")
  na_datasets <- c("Vinther2008", "Longrich2010", "Sansom2010")
  for (ds_name in na_datasets) {
    pd <- inapplicable.phyData[[ds_name]]
    d <- make_ts_data(pd)

    set.seed(4217)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check, info = ds_name)
  }
})

test_that("Wagner on NA datasets is deterministic", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)

  set.seed(6193)
  r1 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  set.seed(6193)
  r2 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})

test_that("Wagner on NA + IW matches fitch_score", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)

  for (k in c(3, 10)) {
    set.seed(8514)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels,
      concavity = k
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels,
      concavity = k
    )
    expect_equal(result$score, fitch_check, tolerance = 1e-6,
                 info = paste("IW k =", k))
  }
})

test_that("Wagner NA tree has valid topology", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)
  n_tip <- length(pd)

  set.seed(2917)
  result <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(nrow(result$edge), 2L * (n_tip - 1L))
  tips <- sort(result$edge[result$edge[, 2] <= n_tip, 2])
  expect_equal(tips, seq_len(n_tip))
})

test_that("Wagner with many addition orders all verify", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  for (s in c(1234, 5678, 9012)) {
    set.seed(s)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check, info = paste("seed", s))
  }
})

test_that("Driven search still finds good scores after Wagner optimization", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(8822)
  result <- TreeSearch:::ts_driven_search(
    d$contrast, d$tip_data, d$weight, d$levels,
    maxReplicates = 2L, targetHits = 1L,
    tbrMaxHits = 1L, ratchetCycles = 2L,
    ratchetPerturbProb = 0.04, driftCycles = 0L,
    driftAfdLimit = 3L, driftRfdLimit = 0.1,
    xssRounds = 0L, xssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 3L, fuseAcceptEqual = FALSE,
    poolMaxSize = 10L, poolSuboptimal = 0.0,
    maxSeconds = 0, verbosity = 0L
  )

  expect_true(result$best_score <= 200)
  expect_true(result$best_score > 0)
})
