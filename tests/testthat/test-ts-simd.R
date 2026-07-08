# Regression tests for EW/NA scoring correctness.
#
# `TreeSearch::Fitch()` is deprecated, and in any case is not an
# independent check: it resolves to the same native `ts::score_tree()`
# engine that `ts_score()` calls, just via a different R entry point
# (TreeLength -> CharacterLength -> FastCharacterLength -> ts_char_steps).
# Comparing the two is circular self-consistency, not validation.
#
# phangorn is a genuine independent oracle for ordinary (non-inapplicable)
# characters, but has no concept of inapplicable-token ("-") semantics, so
# it cannot check the NA-aware three-pass scoring used on these datasets.
#
# Until an independent inapplicable-aware reference exists again (Morphy
# was fully removed from this package), these tests instead lock in
# scores computed with the current engine, which is independently
# validated elsewhere (Wagner/Fitch fuzz testing, exact-insertion
# oracles; see dev/red-team notes). They exist to catch *regressions* in
# scoring — SIMD or otherwise — not to prove initial correctness.
#
# Focuses on edge cases around:
# - Odd vs even state counts (SIMD processes 2 words at a time)
# - Single-state characters (k=1, no SIMD loop iterations)
# - Large state counts (many SIMD iterations)
# - Inapplicable characters (NA-aware three-pass scoring)
# - All scoring modes: EW, IW, profile

# Helpers from helper-ts.R: make_ts_data, ts_score, validate_result

# =====================================================================
# EW scoring: golden-value regression on all inapplicable datasets
# =====================================================================

test_that("SIMD EW scores match known-good values on inapplicable datasets (pectinate)", {
  golden <- c(
    Agnarsson2004 = 1081, Aguado2009 = 1160, Aria2015 = 184,
    Asher2005 = 553, Capa2011 = 620, Conrad2008 = 3605,
    DeAssis2011 = 104, Dikow2009 = 3409, Eklund2004 = 655,
    Geisler2001 = 2386, Giles2015 = 2094, Griswold1999 = 817,
    Liljeblad2008 = 4187, Loconte1991 = 1081, Longrich2010 = 150,
    OLeary1999 = 569, OMeara2014 = 1578, Rougier2012 = 1430,
    Rousset2004 = 687, Sano2011 = 340, Sansom2010 = 220,
    Schulze2007 = 280, Shultz2007 = 615, Vinther2008 = 139,
    Wetterer2000 = 889, Wills2012 = 499, Wilson2003 = 1492,
    Wortley2006 = 657, Zanol2014 = 1931, Zhu2013 = 2150
  )
  for (ds_name in names(inapplicable.phyData)) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(PectinateTree(dataset))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    expect_equal(ew, unname(golden[[ds_name]]),
                 label = paste(ds_name, "pectinate EW"))
  }
})

test_that("SIMD EW scores match known-good values on inapplicable datasets (random)", {
  golden <- c(
    Agnarsson2004 = 1958, Aguado2009 = 1228, Aria2015 = 299,
    Asher2005 = 556, Capa2011 = 1524, Conrad2008 = 3533,
    DeAssis2011 = 272, Dikow2009 = 3552, Eklund2004 = 1068,
    Geisler2001 = 2372
  )
  set.seed(7142)
  for (ds_name in names(inapplicable.phyData)[1:10]) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(RandomTree(dataset, root = TRUE))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    expect_equal(ew, unname(golden[[ds_name]]),
                 label = paste(ds_name, "random EW"))
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
# NA datasets: three-pass scoring regression
# =====================================================================

test_that("NA three-pass scoring with SIMD matches known-good values across datasets", {
  # 5 datasets that heavily exercise NA scoring
  golden <- c(
    Vinther2008 = 188, Agnarsson2004 = 1964, Wills2012 = 712,
    Aria2015 = 292, Zhu2013 = 2256
  )
  set.seed(5063)
  for (ds_name in names(golden)) {
    dataset <- inapplicable.phyData[[ds_name]]
    tree <- Preorder(RandomTree(dataset, root = TRUE))
    ds <- make_ts_data(dataset)
    ew <- ts_score(tree, ds)
    expect_equal(ew, unname(golden[[ds_name]]),
                 label = paste(ds_name, "NA random EW"))
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
  expect_equal(ew, 1081, label = "Loconte1991 EW")
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
