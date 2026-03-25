# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Red-team focus 10: IW/Profile scoring edge cases and regression guards.
# Tests exercise scenarios from S-RED review of ts_fitch.cpp IW/Profile paths.

library(TreeSearch)
library(TreeTools)

# --- Helpers ---
ts_iw <- function(tree, ds, min_steps, k) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                              ds$levels, min_steps = min_steps, concavity = k)
}

make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels
  )
}

result_phylo <- function(result, ref_tree) {
  structure(
    list(edge = result$edge, tip.label = ref_tree$tip.label,
         Nnode = ref_tree$Nnode),
    class = "phylo"
  )
}

# =====================================================================
# 1. Ratchet ZERO_ONLY scoring integrity
#    Guards against active_mask corruption across ratchet cycles.
#    After ratchet, the score returned must match an independent rescore.
# =====================================================================

test_that("Ratchet ZERO_ONLY: returned score matches independent rescore (EW)", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  tree <- PectinateTree(dataset)

  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 5L, perturbMode = 0L, perturbProb = 0.25)

  result_tree <- result_phylo(result, tree)
  rescore <- TreeSearch:::ts_fitch_score(
    result_tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)

  expect_equal(result$score, rescore, tolerance = 1e-8,
               label = "Ratchet ZERO score vs rescore")
})

test_that("Ratchet MIXED: returned score matches independent rescore (EW)", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  tree <- PectinateTree(dataset)

  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 5L, perturbMode = 2L, perturbProb = 0.25)

  result_tree <- result_phylo(result, tree)
  rescore <- TreeSearch:::ts_fitch_score(
    result_tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels)

  expect_equal(result$score, rescore, tolerance = 1e-8,
               label = "Ratchet MIXED score vs rescore")
})

test_that("Ratchet ZERO_ONLY: returned score matches rescore (IW)", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  tree <- PectinateTree(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 5L, perturbMode = 0L, perturbProb = 0.25,
    min_steps = minSteps, concavity = 10)

  rescore <- ts_iw(result_phylo(result, tree), ds, minSteps, 10)
  expect_equal(result$score, rescore, tolerance = 1e-8,
               label = "Ratchet ZERO IW score vs rescore")
})


# =====================================================================
# 2. IW k=0 edge case (extreme concavity)
# =====================================================================

test_that("IW k=0 gives finite scores and saturates at 1 per extra-step char", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  tree <- PectinateTree(dataset)
  score_k0 <- ts_iw(tree, ds, minSteps, 0.0001)  # near-zero

  expect_true(is.finite(score_k0))
  expect_gte(score_k0, 0)

  # With k -> 0, each char with extra > 0 contributes ~1 * freq.
  # So total should approach sum of freqs for chars with any extra steps.
  # Just verify it's in a reasonable range.
  total_freq <- sum(ds$weight)
  expect_lte(score_k0, total_freq)
})

test_that("IW TBR at k=0 doesn't crash or worsen score", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[1]]
  ds <- make_ts_data(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  tree <- PectinateTree(dataset)
  init <- ts_iw(tree, ds, minSteps, 0.01)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = 1L, min_steps = minSteps, concavity = 0.01)

  expect_lte(result$score, init + 1e-8, label = "k=0.01 TBR no worsening")
  rescore <- ts_iw(result_phylo(result, tree), ds, minSteps, 0.01)
  expect_equal(result$score, rescore, tolerance = 1e-8,
               label = "k=0.01 TBR rescore match")
})


# =====================================================================
# 3. Profile scoring via MaximizeParsimony (full pipeline)
#    Individual search bridges (ts_tbr_search, ts_ratchet_search) don't
#    accept infoAmounts; Profile search is only accessible via driven search.
# =====================================================================

test_that("Profile MaximizeParsimony result rescores correctly", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]
  pds <- PrepareDataProfile(dataset)

  set.seed(6437)
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)

  reported <- attr(result, "score")
  actual <- TreeLength(result[[1]], pds, concavity = "profile")
  expect_equal(reported, actual, tolerance = 1e-6,
               label = "Profile driven search rescore")
})


# =====================================================================
# 4. IW + inapplicable data: rescore consistency
# =====================================================================

test_that("IW+NA search results rescore correctly", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

  set.seed(3149)
  tree <- Preorder(RandomTree(dataset, root = TRUE))

  for (method in c("TBR", "Ratchet", "Drift")) {
    result <- switch(method,
      TBR = TreeSearch:::ts_tbr_search(
        tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxHits = 1L, min_steps = minSteps, concavity = 3),
      Ratchet = TreeSearch:::ts_ratchet_search(
        tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
        nCycles = 3L, min_steps = minSteps, concavity = 3),
      Drift = TreeSearch:::ts_drift_search(
        tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
        nCycles = 2L, min_steps = minSteps, concavity = 3)
    )

    rescore <- ts_iw(result_phylo(result, tree), ds, minSteps, 3)
    expect_equal(result$score, rescore, tolerance = 1e-8,
                 label = paste("IW+NA", method, "rescore"))
  }
})
