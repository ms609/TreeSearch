# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for implied weights (IW) scoring in the C++ engine.
# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

ts_iw <- function(tree, ds, min_steps, k) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                              ds$levels, min_steps = min_steps, concavity = k)
}

result_phylo <- function(result, ref_tree) {
  structure(
    list(edge = result$edge, tip.label = ref_tree$tip.label,
         Nnode = ref_tree$Nnode),
    class = "phylo"
  )
}

# =====================================================================
# Hard-coded reference IW scores (C++ engine, verified against morphy
# during development — see AGENTS.md Phase 4 and IW tests).
#
# 6 representative datasets × pectinate tree × k = 3, 10, 100
# Plus random tree (seed 5729) for 3 datasets.
# =====================================================================

# Reference values recomputed 2026-03-19 (T-131) after T-113 NA ambiguity
# bit-stripping fix — `build_dataset()` now strips applicable bits from
# partial {-,X} tokens and preserves full missing data.
#
# `rand` values use BalancedTree (deterministic, platform-independent).
# Previously used RandomTree(seed=5729) which diverged on macOS.
iw_ref <- list(
  Vinther2008 = list(
    ew_pect = 139,
    pect = c(`3` = 15.8714285714, `10` = 6.4955044955, `100` = 0.7641911074),
    ew_rand = 142,
    rand = c(`3` = 16.0142857143, `10` = 6.6543456543, `100` = 0.7921859906)
  ),
  Agnarsson2004 = list(
    ew_pect = 1081,
    pect = c(`3` = 102.7370533878, `10` = 51.8654533610, `100` = 7.4988685706),
    ew_rand = 1117,
    rand = c(`3` = 101.7186105561, `10` = 52.7578507727, `100` = 7.8058662874)
  ),
  Wills2012 = list(
    ew_pect = 499,
    pect = c(`3` = 40.2243589744, `10` = 21.4272698288, `100` = 3.3652284459),
    ew_rand = 516,
    rand = c(`3` = 41.4493714619, `10` = 22.2303189147, `100` = 3.5166229085)
  ),
  Aria2015 = list(
    ew_pect = 184,
    pect = c(`3` = 18.8750000000, `10` = 8.7590840532, `100` = 1.1607895426),
    ew_rand = 196,
    rand = c(`3` = 20.2827380952, `10` = 9.5382957835, `100` = 1.2740063632)
  ),
  Zhu2013 = list(
    ew_pect = 2150,
    pect = c(`3` = 164.0737728728, `10` = 97.9724526960, `100` = 17.1087479711),
    ew_rand = 2186,
    rand = c(`3` = 161.6812051805, `10` = 96.8232223149, `100` = 17.2979327526)
  ),
  Loconte1991 = list(
    ew_pect = 1081,
    pect = c(`3` = 67.0555935288, `10` = 42.2927501720, `100` = 8.1801157451),
    ew_rand = 1055,
    rand = c(`3` = 65.9395171481, `10` = 41.2621955048, `100` = 7.9468418755)
  )
)

# Hard-coded per-pattern step counts (pectinate tree)
steps_ref <- list(
  # Recomputed 2026-03-19 (T-131) after T-113 NA ambiguity bit-stripping fix
  Vinther2008 = as.integer(c(0, 2, 1, 2, 1, 1, 1, 2, 1, 2, 3, 2, 3, 2, 2,
                  4, 4, 3, 3, 5, 2, 2, 2, 0, 3, 3, 3, 5, 3, 2, 2, 4, 2,
                  4, 3, 2, 2, 4, 3, 1, 0, 3, 0, 6, 2, 2, 2, 4, 3, 2)),
  Aria2015 = as.integer(c(2, 7, 2, 2, 9, 2, 3, 3, 6, 2, 4, 3, 2, 5, 2, 2,
               3, 2, 1, 3, 4, 5, 6, 4, 2, 3, 17, 8, 5, 2, 1, 2, 2, 2, 3,
               2, 6, 2, 4, 3, 2, 3, 5, 2, 1, 5, 5, 8, 3, 2))
)


# =====================================================================
# Scoring tests — hard-coded reference agreement
# =====================================================================

test_that("IW pectinate scores match reference for 6 datasets", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in names(iw_ref)) {
    dataset <- inapplicable.phyData[[nm]]
    tree <- TreeTools::PectinateTree(dataset)
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    for (k_str in c("3", "10", "100")) {
      k <- as.numeric(k_str)
      score <- ts_iw(tree, ds, minSteps, k)
      expect_equal(score, iw_ref[[nm]]$pect[[k_str]], tolerance = 1e-8,
                   label = paste(nm, "pect k =", k))
    }
  }
})

test_that("IW balanced-tree scores match reference for 6 datasets", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in names(iw_ref)) {
    dataset <- inapplicable.phyData[[nm]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    # BalancedTree is deterministic and platform-independent (unlike RandomTree
    # whose output depends on the RNG implementation and diverges on macOS).
    tree <- TreeTools::Preorder(TreeTools::BalancedTree(dataset))

    for (k_str in c("3", "10", "100")) {
      k <- as.numeric(k_str)
      score <- ts_iw(tree, ds, minSteps, k)
      expect_equal(score, iw_ref[[nm]]$rand[[k_str]], tolerance = 1e-8,
                   label = paste(nm, "balanced k =", k))
    }
  }
})

test_that("IW k=Inf equals EW score for 6 datasets", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in names(iw_ref)) {
    dataset <- inapplicable.phyData[[nm]]
    tree <- TreeTools::PectinateTree(dataset)
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    ew_score <- ts_score(tree, ds)
    iw_inf <- ts_iw(tree, ds, minSteps, Inf)
    expect_equal(iw_inf, ew_score, label = paste(nm, "k=Inf vs EW"))
    expect_equal(ew_score, iw_ref[[nm]]$ew_pect, label = paste(nm, "EW"))
  }
})


# =====================================================================
# Per-pattern step counts
# =====================================================================

test_that("Per-pattern step counts match reference", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in names(steps_ref)) {
    dataset <- inapplicable.phyData[[nm]]
    tree <- TreeTools::Preorder(TreeTools::PectinateTree(dataset))
    at <- attributes(dataset)
    info <- TreeSearch:::ts_na_char_steps(
      tree$edge, at$contrast,
      matrix(unlist(dataset, use.names = FALSE),
             nrow = length(dataset), byrow = TRUE),
      at$weight, at$levels
    )
    expect_identical(info$steps, steps_ref[[nm]],
                     label = paste(nm, "per-pattern steps"))
  }
})


# =====================================================================
# Edge cases — extreme k values and monotonicity
# =====================================================================

test_that("Extreme k values return finite non-negative scores", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in c("Vinther2008", "Agnarsson2004", "Zhu2013")) {
    dataset <- inapplicable.phyData[[nm]]
    tree <- TreeTools::PectinateTree(dataset)
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    for (k in c(0.01, 0.1, 0.5, 1e4, 1e6)) {
      score <- ts_iw(tree, ds, minSteps, k)
      expect_true(is.finite(score) && score >= 0,
                  label = paste(nm, "k =", k))
    }
  }
})

test_that("IW score decreases monotonically as k increases", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  for (nm in c("Vinther2008", "Agnarsson2004", "Wills2012")) {
    dataset <- inapplicable.phyData[[nm]]
    tree <- TreeTools::PectinateTree(dataset)
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    k_series <- c(0.1, 1, 3, 10, 100, 1000, 1e6)
    scores <- vapply(k_series, function(k) ts_iw(tree, ds, minSteps, k),
                     double(1))
    expect_true(all(diff(scores) <= 1e-10),
                label = paste(nm, "monotonicity"))
  }
})


# =====================================================================
# Search tests — TBR, Ratchet, Drift under IW
# =====================================================================

test_that("IW TBR search improves score and rescores correctly", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  tree <- TreeTools::PectinateTree(dataset)
  initial_iw <- ts_iw(tree, ds, minSteps, 10)

  result <- TreeSearch:::ts_tbr_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxHits = 1L, min_steps = minSteps, concavity = 10)

  expect_lte(result$score, initial_iw)
  rescore <- ts_iw(result_phylo(result, tree), ds, minSteps, 10)
  expect_equal(result$score, rescore, tolerance = 1e-10)
})

test_that("IW ratchet search works", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  tree <- TreeTools::PectinateTree(dataset)
  initial_iw <- ts_iw(tree, ds, minSteps, 10)

  result <- TreeSearch:::ts_ratchet_search(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    nCycles = 3L, min_steps = minSteps, concavity = 10)

  expect_lte(result$score, initial_iw)
})

test_that("IW search results rescore correctly across 3 datasets", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  set.seed(2914)
  for (nm in c("Vinther2008", "Agnarsson2004", "Wills2012")) {
    dataset <- inapplicable.phyData[[nm]]
    ds <- make_ts_data(dataset)
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
    start_tree <- TreeTools::PectinateTree(dataset)

    for (method in c("TBR", "Ratchet", "Drift")) {
      result <- switch(method,
        TBR = TreeSearch:::ts_tbr_search(
          start_tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
          maxHits = 1L, min_steps = minSteps, concavity = 10),
        Ratchet = TreeSearch:::ts_ratchet_search(
          start_tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
          nCycles = 2L, min_steps = minSteps, concavity = 10),
        Drift = TreeSearch:::ts_drift_search(
          start_tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
          nCycles = 2L, min_steps = minSteps, concavity = 10)
      )

      ts_rescore <- ts_iw(result_phylo(result, start_tree), ds, minSteps, 10)
      expect_equal(result$score, ts_rescore, tolerance = 1e-8,
                   label = paste(nm, method, "rescore"))
    }
  }
})

test_that("IW TBR never worsens starting score", {
  skip_on_cran()
  data("inapplicable.phyData", package = "TreeSearch")

  set.seed(7851)
  for (nm in c("Vinther2008", "Agnarsson2004", "Wills2012")) {
    dataset <- inapplicable.phyData[[nm]]
    ds <- make_ts_data(dataset)
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

    for (i in 1:2) {
      tree <- TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))
      init_score <- ts_iw(tree, ds, minSteps, 10)

      result <- TreeSearch:::ts_tbr_search(
        tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxHits = 1L, min_steps = minSteps, concavity = 10)

      expect_lte(result$score, init_score + 1e-8,
                 label = paste(nm, "rep", i))
    }
  }
})
