# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Regression test for the IW NNI incremental-rescore bug fixed in 3df90882
# ("fix(nni): correct IW score computation in incremental rescore").
#
# Before the fix, nni_search() in src/ts_search.cpp computed
#   new_score = best_score + delta
# unconditionally, where `delta` is an integer EW step count from
# fitch_incremental_downpass.  Under IW / profile parsimony, `best_score`
# is a float weighted score, so the addition mixed units and produced
# garbage; accept/reject comparisons were essentially random and the
# `score` returned by ts_nni_search did not match the IW score of the
# returned tree.
#
# These tests pin the contract that ts_nni_search's reported score equals
# the IW score of the returned topology, recomputed independently.

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_nni <- function(tree, ds, maxHits = 20L, concavity = Inf,
                   min_steps = integer(0)) {
  TreeSearch:::ts_nni_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits,
                             min_steps = min_steps,
                             concavity = concavity)
}

test_that("NNI: IW score matches independent recompute (binary)", {
  set.seed(4815)
  mat <- matrix(sample(0:1, 12 * 20, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  min_steps <- as.integer(MinimumLength(dataset, compress = TRUE))

  tree <- as.phylo(7, 12)

  result <- ts_nni(tree, ds, concavity = 10.0, min_steps = min_steps,
                   maxHits = 5L)
  validate_result(result, 12L)

  rt <- tree
  rt$edge <- result$edge

  # C++ recompute under same IW configuration
  c_score <- ts_score(rt, ds, concavity = 10.0, min_steps = min_steps)
  expect_equal(result$score, c_score, tolerance = 1e-8)

  # R-level TreeLength as second independent oracle
  r_score <- TreeLength(rt, dataset, concavity = 10)
  expect_equal(result$score, r_score, tolerance = 1e-8)
})

test_that("NNI: IW score matches independent recompute (multistate)", {
  set.seed(2306)
  mat <- matrix(sample(0:3, 12 * 18, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  min_steps <- as.integer(MinimumLength(dataset, compress = TRUE))

  tree <- as.phylo(123, 12)

  result <- ts_nni(tree, ds, concavity = 10.0, min_steps = min_steps,
                   maxHits = 5L)
  validate_result(result, 12L)

  rt <- tree
  rt$edge <- result$edge

  c_score <- ts_score(rt, ds, concavity = 10.0, min_steps = min_steps)
  expect_equal(result$score, c_score, tolerance = 1e-8)

  r_score <- TreeLength(rt, dataset, concavity = 10)
  expect_equal(result$score, r_score, tolerance = 1e-8)
})

test_that("NNI: IW score matches independent recompute across concavities", {
  # The bug returned garbage regardless of concavity value; sweep a few to
  # make sure the fix holds for both tight and loose weighting.
  set.seed(9182)
  mat <- matrix(sample(0:1, 12 * 20, replace = TRUE),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  min_steps <- as.integer(MinimumLength(dataset, compress = TRUE))
  tree <- as.phylo(31, 12)

  for (k in c(3.0, 10.0, 100.0)) {
    result <- ts_nni(tree, ds, concavity = k, min_steps = min_steps,
                     maxHits = 5L)
    validate_result(result, 12L)

    rt <- tree
    rt$edge <- result$edge
    c_score <- ts_score(rt, ds, concavity = k, min_steps = min_steps)
    expect_equal(result$score, c_score, tolerance = 1e-8,
                 label = paste0("concavity=", k))
  }
})

test_that("NNI: IW score matches independent recompute with NA tokens", {
  # Combine IW with inapplicable tokens — exercises the IW accept-path on
  # the NA-aware scoring branch.
  set.seed(7401)
  mat <- matrix(sample(c("0", "1", "-"), 12 * 20, replace = TRUE,
                       prob = c(0.45, 0.45, 0.10)),
                nrow = 12, dimnames = list(paste0("t", 1:12), NULL))
  dataset <- MatrixToPhyDat(mat)
  ds <- make_ts_data(dataset)
  min_steps <- as.integer(MinimumLength(dataset, compress = TRUE))

  tree <- as.phylo(55, 12)

  result <- ts_nni(tree, ds, concavity = 10.0, min_steps = min_steps,
                   maxHits = 5L)
  validate_result(result, 12L)

  rt <- tree
  rt$edge <- result$edge

  c_score <- ts_score(rt, ds, concavity = 10.0, min_steps = min_steps)
  expect_equal(result$score, c_score, tolerance = 1e-8)
})
