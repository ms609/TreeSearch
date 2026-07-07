# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Internal C++ scoring function (not exported)
ts_fitch_score <- TreeSearch:::ts_fitch_score

test_that("C++ profile score matches R-level TreeLength", {
  data("congreveLamsdellMatrices", package = "TreeSearch")

  # Test across multiple datasets
  for (ds_idx in c(1, 5, 10, 20, 42)) {
    dataset <- congreveLamsdellMatrices[[ds_idx]]
    pds <- PrepareDataProfile(dataset)
    at <- attributes(pds)
    tip_data <- matrix(unlist(pds, use.names = FALSE),
                       nrow = length(pds), byrow = TRUE)

    for (seed in c(3017, 5539, 9281)) {
      set.seed(seed + ds_idx)
      tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)

      rScore <- TreeLength(tree, pds, concavity = "profile")

      tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(pds)))
      if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
        tree2 <- TreeTools::RootTree(tree2, 1L)
      }

      cScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                                at$weight, at$levels,
                                infoAmounts = at$info.amounts)

      expect_equal(cScore, rScore, tolerance = 1e-8,
                   label = paste0("ds=", ds_idx, " seed=", seed))
    }
  }
})

test_that("EW scoring unchanged when infoAmounts not provided", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]
  at <- attributes(dataset)
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)

  set.seed(6203)
  tree <- TreeTools::RootTree(TreeTools::RandomTree(dataset), 1L)
  tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(dataset)))
  if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
    tree2 <- TreeTools::RootTree(tree2, 1L)
  }

  # EW: no infoAmounts
  ewScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                             at$weight, at$levels)
  rScore <- TreeLength(tree, dataset)
  expect_equal(ewScore, rScore)
})

test_that("IW scoring unchanged when infoAmounts not provided", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]
  at <- attributes(dataset)
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  min_steps <- MinimumLength(dataset, compress = TRUE)

  set.seed(7814)
  tree <- TreeTools::RootTree(TreeTools::RandomTree(dataset), 1L)
  tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(dataset)))
  if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
    tree2 <- TreeTools::RootTree(tree2, 1L)
  }

  iwScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                             at$weight, at$levels,
                             min_steps = as.integer(min_steps),
                             concavity = 10.0)
  rScore <- TreeLength(tree, dataset, concavity = 10)
  expect_equal(iwScore, rScore, tolerance = 1e-8)
})

test_that("MaximizeParsimony with profile returns valid trees", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]

  set.seed(4488)
  # collapse = FALSE: the score check below re-scores result[[1]] with the
  # binary-only TreeLength(), so it needs a fully-resolved tree (the default
  # collapse = TRUE returns polytomies).
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 3L, targetHits = 2L,
                               verbosity = 0L, collapse = FALSE)

  expect_s3_class(result, "multiPhylo")
  expect_true(length(result) >= 1L)
  expect_true(is.finite(attr(result, "score")))

  # Verify reported score matches TreeLength
  pds <- PrepareDataProfile(dataset)
  reported <- attr(result, "score")
  actual <- TreeLength(result[[1]], pds, concavity = "profile")
  expect_equal(reported, actual, tolerance = 1e-6)
})

test_that("Profile search improves or equals starting score", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]
  pds <- PrepareDataProfile(dataset)

  set.seed(2977)
  startTree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
  startScore <- TreeLength(startTree, pds, concavity = "profile")

  set.seed(2977)
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 5L, targetHits = 2L,
                               verbosity = 0L)

  searchScore <- attr(result, "score")
  expect_true(searchScore <= startScore + 1e-8)
})

test_that("Profile scoring handles inapplicable datasets", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (dsName in c("Vinther2008", "Sansom2010")) {
    dataset <- inapplicable.phyData[[dsName]]
    pds <- suppressMessages(PrepareDataProfile(dataset))
    at <- attributes(pds)

    if (length(at$weight) == 0L || attr(pds, "nr") == 0L) next

    tip_data <- matrix(unlist(pds, use.names = FALSE),
                       nrow = length(pds), byrow = TRUE)

    set.seed(1042)
    tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
    tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(pds)))
    if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
      tree2 <- TreeTools::RootTree(tree2, 1L)
    }

    rScore <- TreeLength(tree, pds, concavity = "profile")
    cScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                              at$weight, at$levels,
                              infoAmounts = at$info.amounts)

    expect_equal(cScore, rScore, tolerance = 1e-8,
                 label = paste0("dataset=", dsName))
  }
})

test_that("Profile driven search is reproducible with set.seed", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]

  set.seed(6701)
  result1 <- MaximizeParsimony(dataset, concavity = "profile",
                                maxReplicates = 2L, targetHits = 1L,
                                verbosity = 0L)
  set.seed(6701)
  result2 <- MaximizeParsimony(dataset, concavity = "profile",
                                maxReplicates = 2L, targetHits = 1L,
                                verbosity = 0L)

  expect_equal(attr(result1, "score"), attr(result2, "score"))
})

# --- Multi-state profile parsimony integration tests (T-104) ----------------

test_that("TreeLength with profile scoring on multi-state data", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Longrich2010"]]
  pds <- PrepareDataProfile(dataset)

  set.seed(8347)
  tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
  score <- TreeLength(tree, pds, concavity = "profile")

  expect_true(is.finite(score))
  expect_gt(score, 0)
})

test_that("C++ and R-level profile scores agree on multi-state data", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (dsName in c("Longrich2010", "Vinther2008")) {
    dataset <- inapplicable.phyData[[dsName]]
    pds <- PrepareDataProfile(dataset)
    at <- attributes(pds)

    if (length(at$weight) == 0L || attr(pds, "nr") == 0L) next

    tip_data <- matrix(unlist(pds, use.names = FALSE),
                       nrow = length(pds), byrow = TRUE)

    set.seed(5581)
    tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
    tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(pds)))
    if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
      tree2 <- TreeTools::RootTree(tree2, 1L)
    }

    rScore <- TreeLength(tree, pds, concavity = "profile")
    cScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                              at$weight, at$levels,
                              infoAmounts = at$info.amounts)

    expect_equal(cScore, rScore, tolerance = 1e-8,
                 label = paste0("dataset=", dsName))
  }
})

test_that("MaximizeParsimony profile search works with multi-state data", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Longrich2010"]]

  set.seed(3692)
  # collapse = FALSE: the score check below re-scores result[[1]] with the
  # binary-only TreeLength(), so it needs a fully-resolved tree (Longrich2010 has
  # zero-length branches, which the default collapse = TRUE would contract).
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 3L, targetHits = 2L,
                               verbosity = 0L, collapse = FALSE)

  expect_s3_class(result, "multiPhylo")
  expect_true(length(result) >= 1L)
  expect_true(is.finite(attr(result, "score")))

  pds <- PrepareDataProfile(dataset)
  reported <- attr(result, "score")
  actual <- TreeLength(result[[1]], pds, concavity = "profile")
  expect_equal(reported, actual, tolerance = 1e-6)
})

test_that("Profile search improves score on multi-state data", {
  data("inapplicable.phyData", package = "TreeSearch")
  dataset <- inapplicable.phyData[["Vinther2008"]]
  pds <- PrepareDataProfile(dataset)

  set.seed(7204)
  startTree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
  startScore <- TreeLength(startTree, pds, concavity = "profile")

  set.seed(7204)
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 3L, targetHits = 2L,
                               verbosity = 0L)

  searchScore <- attr(result, "score")
  expect_true(searchScore <= startScore + 1e-8)
})

test_that("Infeasible multi-state chars reduced to binary in PrepareDataProfile", {
  # Sun2018-like dataset with 3+ state characters and many tips.
  # Without the feasibility guard, PrepareDataProfile would hang.
  sun_file <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  sun <- suppressWarnings(ReadAsPhyDat(sun_file))

  # Should complete in reasonable time (< 10 s) with warning suppressed
  pds <- suppressWarnings(PrepareDataProfile(sun))

  expect_true(!is.null(attr(pds, "info.amounts")))
  expect_true(ncol(attr(pds, "info.amounts")) > 0)

  set.seed(1934)
  tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
  score <- TreeLength(tree, pds, concavity = "profile")
  expect_true(is.finite(score))
  expect_gt(score, 0)

  # Profile search also works end-to-end
  set.seed(1934)
  result <- suppressWarnings(MaximizeParsimony(
    sun, concavity = "profile",
    maxReplicates = 1L, targetHits = 1L,
    maxSeconds = 30, verbosity = 0L
  ))
  expect_true(is.finite(attr(result, "score")))
})

test_that("Binary-only dataset: profile scores unchanged by multi-state code", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  dataset <- congreveLamsdellMatrices[[10]]
  pds <- PrepareDataProfile(dataset)
  at <- attributes(pds)

  tip_data <- matrix(unlist(pds, use.names = FALSE),
                     nrow = length(pds), byrow = TRUE)

  set.seed(9412)
  tree <- TreeTools::RootTree(TreeTools::RandomTree(pds), 1L)
  tree2 <- TreeTools::Preorder(TreeTools::RenumberTips(tree, names(pds)))
  if (tree2[["edge"]][1, 2] > TreeTools::NTip(tree2)) {
    tree2 <- TreeTools::RootTree(tree2, 1L)
  }

  rScore <- TreeLength(tree, pds, concavity = "profile")
  cScore <- ts_fitch_score(tree2[["edge"]], at$contrast, tip_data,
                            at$weight, at$levels,
                            infoAmounts = at$info.amounts)

  expect_equal(cScore, rScore, tolerance = 1e-8)

  # Verify search also works identically
  set.seed(9412)
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 2L, targetHits = 1L,
                               verbosity = 0L)
  expect_true(is.finite(attr(result, "score")))
})

test_that("profile search with Monte Carlo info is reproducible under set.seed()", {
  # Regression: the Monte Carlo profile-information estimate (used when the
  # exact solver is infeasible) drew random trees from an unseeded generator,
  # so two MaximizeParsimony(concavity = "profile") runs with the same seed
  # returned different scores.  Forcing profile_approx = "mc" exercises that
  # path on a small dataset; the same seed must now give an identical result,
  # and a different seed a (near-certainly) different one.  Replicate-bounded
  # (maxSeconds = 0) so the only remaining randomness is the seeded RNG.
  data("inapplicable.phyData", package = "TreeSearch")
  phy0 <- inapplicable.phyData[["Vinther2008"]]
  m <- TreeTools::PhyDatToMatrix(phy0, ambigNA = FALSE)
  m[m == "-"] <- "?"                       # pure Fitch (has_na = FALSE)
  phy <- TreeTools::MatrixToPhyDat(m)

  searchScore <- function(seed) {
    set.seed(seed)
    min(attr(suppressWarnings(suppressMessages(
      MaximizeParsimony(phy, concavity = "profile", profile_approx = "mc",
                        maxReplicates = 2L, maxSeconds = 0, nThreads = 1L,
                        verbosity = 0L))), "score"))
  }

  s_a <- searchScore(99L)
  s_b <- searchScore(99L)
  s_c <- searchScore(42L)
  expect_equal(s_a, s_b)                          # same seed -> identical
  expect_false(isTRUE(all.equal(s_a, s_c)))       # different seed -> different
})
