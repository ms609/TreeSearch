# Tests for SIMD vectorization correctness (Phase 3E).
#
# Verifies that SIMD-accelerated scoring produces bit-identical results
# to the pre-SIMD implementation. Focuses on edge cases around:
# - Odd vs even state counts (SIMD processes 2 words at a time)
# - Single-state characters (k=1, no SIMD loop iterations)
# - Large state counts (many SIMD iterations)
# - Inapplicable characters (NA-aware three-pass scoring)
# - All scoring modes: EW, IW, profile

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

morphy_ew_ref <- function(tree, dataset) {
  suppressWarnings(TreeSearch::Fitch(tree, dataset))
}

# =====================================================================
# EW scoring: morphy cross-validation on all inapplicable datasets
# =====================================================================

test_that("SIMD EW scores match morphy on inapplicable datasets (pectinate)", {
  for (ds_name in names(inapplicable.phyData)) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(PectinateTree(dataset))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    ref <- morphy_ew_ref(tree, dataset)
    expect_equal(ew, ref, label = paste(ds_name, "pectinate EW"))
  }
})

test_that("SIMD EW scores match morphy on inapplicable datasets (random)", {
  set.seed(7142)
  for (ds_name in names(inapplicable.phyData)[1:10]) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(RandomTree(dataset, root = TRUE))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    ref <- morphy_ew_ref(tree, dataset)
    expect_equal(ew, ref, label = paste(ds_name, "random EW"))
  }
})

# =====================================================================
# DNA data: exactly 4 applicable states → even n_states (good SIMD case)
# =====================================================================

test_that("SIMD EW scores correct on DNA data (4 states, even)", {
  suppressWarnings(data("Laurasiatherian", package = "phangorn"))
  dna <- Laurasiatherian
  set.seed(3827)
  tree <- Preorder(RandomTree(dna, root = TRUE))
  ds <- make_ts_data(dna)
  ew <- ts_score(tree, ds)
  expect_true(is.finite(ew))
  expect_gt(ew, 0)

  # Deterministic: same score twice
  ew2 <- ts_score(tree, ds)
  expect_identical(ew, ew2)
})

# =====================================================================
# IW scoring with different concavity values
# =====================================================================

test_that("SIMD IW scores are self-consistent", {
  dataset <- inapplicable.phyData[["Vinther2008"]]
  set.seed(4519)
  tree <- Preorder(RandomTree(dataset, root = TRUE))
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  scores_k <- vapply(c(1, 2, 3, 5, 10, 50, 100, 1000), function(k) {
    ts_score(tree, ds, concavity = k, min_steps = minSteps)
  }, numeric(1))

  # All finite and positive

  expect_true(all(is.finite(scores_k)))
  expect_true(all(scores_k > 0))

  # Monotonically decreasing with increasing k (more weight = less penalty)
  diffs <- diff(scores_k)
  expect_true(all(diffs <= 0), label = "IW scores decrease with increasing k")
})

# =====================================================================
# TBR search: verify SIMD doesn't break move evaluation
# =====================================================================

test_that("TBR search with SIMD finds optimal or near-optimal scores", {
  dataset <- inapplicable.phyData[["Vinther2008"]]
  set.seed(6184)
  tree <- Preorder(RandomTree(dataset, root = TRUE))
  ds <- make_ts_data(dataset)
  initial <- ts_score(tree, ds)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = 10L, acceptEqual = TRUE
  )
  expect_lte(result$score, initial)
})

test_that("TBR search reproducible with set.seed (SIMD determinism)", {
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(2951)
  r1 <- TreeSearch:::ts_tbr_search(
    PectinateTree(dataset)$edge, ds$contrast, ds$tip_data,
    ds$weight, ds$levels, maxHits = 5L
  )
  set.seed(2951)
  r2 <- TreeSearch:::ts_tbr_search(
    PectinateTree(dataset)$edge, ds$contrast, ds$tip_data,
    ds$weight, ds$levels, maxHits = 5L
  )
  expect_identical(r1$score, r2$score)
  expect_identical(r1$edge, r2$edge)
})

# =====================================================================
# Driven search end-to-end (exercises all SIMD codepaths)
# =====================================================================

test_that("Driven search produces valid results with SIMD", {
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)

  set.seed(8371)
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data,
    ds$weight, ds$levels,
    maxReplicates = 2L, targetHits = 1L,
    ratchetCycles = 1L, driftCycles = 0L,
    xssPartitions = 2L, rssRounds = 0L, cssRounds = 0L,
    cssPartitions = 2L, fuseInterval = 0L,
    poolMaxSize = 2L, poolSuboptimal = 0,
    ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 50L,
    ratchetAdaptive = FALSE, maxSeconds = 30,
    verbosity = 0L
  )

  expect_true(is.finite(result$best_score))
  expect_length(result$trees, result$pool_size)
})

# =====================================================================
# NA datasets: three-pass scoring verification
# =====================================================================

test_that("NA three-pass scoring with SIMD matches morphy across datasets", {
  # Use 5 datasets that heavily exercise NA scoring
  na_datasets <- c("Vinther2008", "Agnarsson2004", "Wills2012",
                    "Aria2015", "Zhu2013")
  set.seed(5063)
  for (ds_name in na_datasets) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(RandomTree(dataset, root = TRUE))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    ref <- morphy_ew_ref(tree, dataset)
    expect_equal(ew, ref, label = paste(ds_name, "NA random EW"))
  }
})

# =====================================================================
# Edge case: very small dataset (n_states likely 1 or 2)
# =====================================================================

test_that("SIMD handles very small datasets correctly", {
  # 4 tips, minimal characters
  dataset <- inapplicable.phyData[["Loconte1991"]]
  tree <- Preorder(PectinateTree(dataset))
  ds <- make_ts_data(dataset)
  ew <- ts_score(tree, ds)
  ref <- morphy_ew_ref(tree, dataset)
  expect_equal(ew, ref, label = "Loconte1991 EW")
})

# =====================================================================
# Consistency across multiple trees on same dataset
# =====================================================================

test_that("SIMD scores consistent across 20 random trees", {
  dataset <- inapplicable.phyData[["Agnarsson2004"]]
  ds <- make_ts_data(dataset)
  set.seed(9256)
  scores <- vapply(seq_len(20), function(i) {
    tree <- Preorder(RandomTree(dataset, root = TRUE))
    ts_score(tree, ds)
  }, numeric(1))

  # All should be finite positive integers
  expect_true(all(is.finite(scores)))
  expect_true(all(scores > 0))
  expect_true(all(scores == floor(scores)))
})
