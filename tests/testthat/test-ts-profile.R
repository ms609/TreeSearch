context("Profile parsimony (C++ engine)")

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
  result <- MaximizeParsimony(dataset, concavity = "profile",
                               maxReplicates = 3L, targetHits = 2L,
                               verbosity = 0L)

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
    pds <- PrepareDataProfile(dataset)
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
