# Tests for implied weights (IW) scoring in the C++ engine.
#
# Verifies that ts_fitch_score with finite concavity matches
# the reference morphy_iw implementation across all datasets,
# tree topologies, and concavity values.

library("TreeTools")

make_ts_data <- function(pd) {
  list(
    contrast = attr(pd, "contrast"),
    tip_data = t(vapply(pd, I, pd[[1]])),
    weight   = attr(pd, "weight"),
    levels   = attr(pd, "levels")
  )
}

ts_iw <- function(tree, ds, min_steps, k) {
  ts_fitch_score(tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
                 min_steps = min_steps, concavity = k)
}

# Compute morphy-based IW score for reference
morphy_iw_ref <- function(tree, dataset, k) {
  tree <- TreeTools::Preorder(tree)
  characters <- TreeTools::PhyToString(dataset, ps = "", useIndex = FALSE,
                                      byTaxon = FALSE, concatenate = FALSE)
  morphyObjs <- lapply(characters, SingleCharMorphy)
  on.exit(vapply(morphyObjs, UnloadMorphy, integer(1)), add = TRUE)
  weight <- attr(dataset, "weight")
  minLength <- MinimumLength(dataset, compress = TRUE)
  charSeq <- seq_along(characters) - 1L
  morphy_iw(tree$edge, morphyObjs, weight, minLength, charSeq,
            k, target = Inf)
}

# Build result phylo from search output
result_phylo <- function(result, ref_tree) {
  structure(
    list(edge = result$edge, tip.label = ref_tree$tip.label,
         Nnode = ref_tree$Nnode),
    class = "phylo"
  )
}

# =====================================================================
# Scoring tests — morphy reference agreement
# =====================================================================

test_that("IW score matches morphy_iw on Vinther2008", {
  dataset <- inapplicable.phyData[["Vinther2008"]]
  tree <- TreeTools::PectinateTree(dataset)

  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  for (k in c(2, 3, 10, 100)) {
    ts_val <- ts_iw(tree, ds, minSteps, k)
    ref_val <- morphy_iw_ref(tree, dataset, k)
    expect_equal(ts_val, ref_val, tolerance = 1e-8,
                 label = paste("Vinther2008 k =", k))
  }
})

test_that("IW score matches morphy_iw on Agnarsson2004", {
  dataset <- inapplicable.phyData[["Agnarsson2004"]]
  tree <- TreeTools::PectinateTree(dataset)

  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  for (k in c(3, 10)) {
    ts_val <- ts_iw(tree, ds, minSteps, k)
    ref_val <- morphy_iw_ref(tree, dataset, k)
    expect_equal(ts_val, ref_val, tolerance = 1e-8,
                 label = paste("Agnarsson2004 k =", k))
  }
})

test_that("IW score with k=Inf equals EW score for all datasets", {
  for (ds_name in names(inapplicable.phyData)) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- TreeTools::PectinateTree(dataset)
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    ew_score <- ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                               ds$weight, ds$levels)
    iw_inf <- ts_iw(tree, ds, minSteps, Inf)
    expect_equal(iw_inf, ew_score, label = paste(ds_name, "k=Inf vs EW"))
  }
})

test_that("IW matches morphy_iw on random trees across datasets", {
  skip_on_cran()

  datasets_to_test <- c("Vinther2008", "Agnarsson2004", "Wills2012")

  for (ds_name in datasets_to_test) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    set.seed(5729)
    for (i in 1:3) {
      tree <- TreeTools::RandomTree(dataset, root = TRUE)
      tree <- TreeTools::Preorder(tree)

      for (k in c(3, 10)) {
        ts_val <- ts_iw(tree, ds, minSteps, k)
        ref_val <- morphy_iw_ref(tree, dataset, k)
        expect_equal(ts_val, ref_val, tolerance = 1e-8,
                     label = paste(ds_name, "tree", i, "k =", k))
      }
    }
  }
})

test_that("IW matches morphy across all 30 datasets", {
  skip_on_cran()

  set.seed(3847)
  for (ds_name in names(inapplicable.phyData)) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)

    trees <- list(
      TreeTools::PectinateTree(dataset),
      TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE)),
      TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))
    )

    for (ti in seq_along(trees)) {
      for (k in c(2, 3, 10, 100)) {
        ts_val <- ts_iw(trees[[ti]], ds, minSteps, k)
        ref_val <- morphy_iw_ref(trees[[ti]], dataset, k)
        expect_equal(ts_val, ref_val, tolerance = 1e-8,
                     label = paste(ds_name, "tree", ti, "k =", k))
      }
    }
  }
})

# =====================================================================
# Per-pattern step counts match morphy
# =====================================================================

test_that("Per-pattern step counts match morphy across all datasets", {
  skip_on_cran()

  set.seed(6192)
  for (ds_name in names(inapplicable.phyData)) {
    dataset <- inapplicable.phyData[[ds_name]]
    at <- attributes(dataset)
    contrast <- at$contrast
    tip_data <- matrix(unlist(dataset, use.names = FALSE),
                       nrow = length(dataset), byrow = TRUE)
    wt <- at$weight
    lvls <- at$levels

    trees <- list(
      TreeTools::Preorder(TreeTools::PectinateTree(dataset)),
      TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))
    )

    characters <- TreeTools::PhyToString(dataset, ps = "", useIndex = FALSE,
                                         byTaxon = FALSE, concatenate = FALSE)

    for (ti in seq_along(trees)) {
      tree <- trees[[ti]]

      info <- ts_na_char_steps(tree$edge, contrast, tip_data, wt, lvls)
      ts_steps <- info$steps

      morphyObjs <- lapply(characters, SingleCharMorphy)
      morphy_steps <- preorder_morphy_by_char(tree$edge, morphyObjs)
      vapply(morphyObjs, UnloadMorphy, integer(1))

      expect_identical(
        ts_steps, morphy_steps,
        label = paste(ds_name, "tree", ti, "per-pattern steps")
      )
    }
  }
})

# =====================================================================
# Edge cases — extreme k values and monotonicity
# =====================================================================

test_that("Extreme k values return finite non-negative scores", {
  skip_on_cran()

  test_ds <- c("Vinther2008", "Agnarsson2004", "Conrad2008", "Zhu2013")
  set.seed(4523)

  for (ds_name in test_ds) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)
    tree <- TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))

    for (k in c(0.01, 0.1, 0.5, 1e4, 1e6)) {
      score <- ts_iw(tree, ds, minSteps, k)
      expect_true(is.finite(score) && score >= 0,
                  label = paste(ds_name, "k =", k))
    }
  }
})

test_that("IW score decreases monotonically as k increases", {
  skip_on_cran()

  test_ds <- c("Vinther2008", "Agnarsson2004", "Wills2012", "Conrad2008")

  for (ds_name in test_ds) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- MinimumLength(dataset, compress = TRUE)
    tree <- TreeTools::PectinateTree(dataset)

    k_series <- c(0.1, 1, 3, 10, 100, 1000, 1e6)
    scores <- vapply(k_series, function(k) ts_iw(tree, ds, minSteps, k),
                     double(1))

    # e/(k+e) decreases as k increases, so IW score decreases
    expect_true(all(diff(scores) <= 1e-10),
                label = paste(ds_name, "monotonicity"))
  }
})

# =====================================================================
# Search tests — TBR, Ratchet, Drift under IW
# =====================================================================

test_that("IW TBR search improves score", {
  skip_on_cran()

  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  tree <- TreeTools::PectinateTree(dataset)
  initial_iw <- ts_iw(tree, ds, minSteps, 10)

  result <- ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                          ds$levels, maxHits = 1L,
                          min_steps = minSteps, concavity = 10)

  expect_lte(result$score, initial_iw)

  # Verify the returned score matches a rescore of the result tree
  rescore <- ts_iw(result_phylo(result, tree), ds, minSteps, 10)
  expect_equal(result$score, rescore, tolerance = 1e-10)
})

test_that("IW TBR result score matches morphy_iw", {
  skip_on_cran()

  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  tree <- TreeTools::PectinateTree(dataset)
  result <- ts_tbr_search(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                          ds$levels, maxHits = 1L,
                          min_steps = minSteps, concavity = 10)

  ref_iw <- morphy_iw_ref(result_phylo(result, tree), dataset, 10)
  expect_equal(result$score, ref_iw, tolerance = 1e-8)
})

test_that("IW ratchet search works", {
  skip_on_cran()

  dataset <- inapplicable.phyData[["Vinther2008"]]
  ds <- make_ts_data(dataset)
  minSteps <- MinimumLength(dataset, compress = TRUE)

  tree <- TreeTools::PectinateTree(dataset)
  initial_iw <- ts_iw(tree, ds, minSteps, 10)

  result <- ts_ratchet_search(tree$edge, ds$contrast, ds$tip_data, ds$weight,
                              ds$levels, nCycles = 3L,
                              min_steps = minSteps, concavity = 10)

  expect_lte(result$score, initial_iw)
})

test_that("IW search results rescore correctly across datasets and methods", {
  skip_on_cran()

  test_datasets <- c("Vinther2008", "Agnarsson2004", "Wills2012",
                     "Aria2015", "Conrad2008", "Loconte1991",
                     "Griswold1999", "Zhu2013", "Schulze2007", "Sano2011")
  set.seed(2914)

  for (ds_name in test_datasets) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
    start_tree <- TreeTools::PectinateTree(dataset)

    for (method in c("TBR", "Ratchet", "Drift")) {
      result <- switch(method,
        TBR = ts_tbr_search(start_tree$edge, ds$contrast, ds$tip_data,
                            ds$weight, ds$levels, maxHits = 1L,
                            min_steps = minSteps, concavity = 10),
        Ratchet = ts_ratchet_search(start_tree$edge, ds$contrast, ds$tip_data,
                                    ds$weight, ds$levels, nCycles = 2L,
                                    min_steps = minSteps, concavity = 10),
        Drift = ts_drift_search(start_tree$edge, ds$contrast, ds$tip_data,
                                ds$weight, ds$levels, nCycles = 2L,
                                min_steps = minSteps, concavity = 10)
      )

      # Independent ts rescore must match
      ts_rescore <- ts_iw(result_phylo(result, start_tree), ds, minSteps, 10)
      expect_equal(result$score, ts_rescore, tolerance = 1e-8,
                   label = paste(ds_name, method, "ts rescore"))

      # morphy rescore must match
      ref <- morphy_iw_ref(result_phylo(result, start_tree), dataset, 10)
      expect_equal(result$score, ref, tolerance = 1e-8,
                   label = paste(ds_name, method, "morphy rescore"))
    }
  }
})

test_that("IW TBR never worsens starting score", {
  skip_on_cran()

  test_datasets <- c("Vinther2008", "Agnarsson2004", "Wills2012",
                     "Conrad2008", "Zhu2013")
  set.seed(7851)

  for (ds_name in test_datasets) {
    dataset <- inapplicable.phyData[[ds_name]]
    ds <- make_ts_data(dataset)
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))

    for (i in 1:3) {
      tree <- TreeTools::Preorder(TreeTools::RandomTree(dataset, root = TRUE))
      init_score <- ts_iw(tree, ds, minSteps, 10)

      result <- ts_tbr_search(tree$edge, ds$contrast, ds$tip_data,
                              ds$weight, ds$levels, maxHits = 1L,
                              min_steps = minSteps, concavity = 10)

      expect_lte(result$score, init_score + 1e-8,
                 label = paste(ds_name, "rep", i))
    }
  }
})
