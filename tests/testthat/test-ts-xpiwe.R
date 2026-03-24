# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Extended implied weighting (Goloboff 2014) — missing-entries correction.
# Tests verify per-character adjusted concavity via the extrapolation factor.

library(TreeSearch)
library(TreeTools)

# Shorthand for internal C++ bridge
ts_score <- function(...) TreeSearch:::ts_fitch_score(...)

# --- Helper: build a phyDat with controllable missing entries ---
# 8 taxa, 6 characters.
# Chars 1-3: fully observed; Chars 4-6: taxa G,H are "?".
# Tree chosen so chars 4-6 have homoplasy > 1 (needed for XPIWE to differ,
# since Φ-rescaling normalises the first step identically).
make_missing_data <- function() {
  mat <- matrix(c(
    # C1 C2 C3 C4 C5 C6
    "0", "0", "1", "0", "1", "0",  # A
    "0", "1", "0", "1", "0", "1",  # B
    "1", "0", "1", "0", "1", "0",  # C
    "1", "1", "0", "1", "0", "1",  # D
    "0", "0", "1", "1", "0", "1",  # E
    "1", "1", "0", "0", "1", "0",  # F
    "0", "1", "1", "?", "?", "?",  # G — missing C4-6
    "1", "0", "0", "?", "?", "?"   # H — missing C4-6
  ), nrow = 8, byrow = TRUE,
  dimnames = list(LETTERS[1:8], paste0("c", 1:6)))
  phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))
}

make_tree <- function() {
  ape::read.tree(text = "((((A,B),(C,D)),E),(F,(G,H)));")
}

test_that("No-missing regression: XPIWE = IW on complete data", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  ds <- congreveLamsdellMatrices[[1]]
  tree <- TreeTools::BalancedTree(ds)

  score_iw  <- TreeLength(tree, ds, concavity = 10, extended_iw = FALSE)
  score_xp  <- TreeLength(tree, ds, concavity = 10, extended_iw = TRUE)

  # No missing data: f = 1 for all chars => eff_k = k => XPIWE = IW
  expect_equal(score_xp, score_iw)
})

test_that("XPIWE differs from standard IW on data with missing entries", {
  ds <- make_missing_data()
  tree <- make_tree()

  score_iw  <- TreeLength(tree, ds, concavity = 10, extended_iw = FALSE)
  # With r=1: f = 1 + 1*2/6 = 1.333, eff_k = 10/1.333 = 7.5.
  # XPIWE only differs from IW for chars with homoplasy > 1 AND missing data.
  score_xp_r1 <- TreeLength(tree, ds, concavity = 10, extended_iw = TRUE,
                             xpiwe_r = 1.0)
  expect_false(isTRUE(all.equal(score_xp_r1, score_iw)))
})

test_that("Hand-computed XPIWE score matches R-level calculation", {
  ds <- make_missing_data()
  tree <- make_tree()
  k <- 10
  r <- 1.0
  max_f <- 5.0

  # Compute expected score from first principles
  obsCount <- TreeSearch:::.ObsCount(ds)
  nTaxa <- length(ds)

  # Goloboff (2014) Extension 3: f = 1 + r * missing / obs
  f <- pmin(pmax(1 + r * (nTaxa - obsCount) / obsCount, 1), max_f)
  eff_k <- k / f
  phi <- (1 + eff_k) / (1 + k)

  ds_iw <- TreeSearch:::PrepareDataIW(ds)
  at <- attributes(ds_iw)

  steps <- TreeSearch::CharacterLength(tree, ds_iw, compress = TRUE)
  minLength <- at[["min.length"]]
  weight <- at[["weight"]]
  homoplasies <- steps - minLength
  expected <- sum(homoplasies / (homoplasies + eff_k) * weight * phi)

  actual <- TreeLength(tree, ds, concavity = k, extended_iw = TRUE,
                       xpiwe_r = r, xpiwe_max_f = max_f)

  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("xpiwe_r and xpiwe_max_f affect scores", {
  ds <- make_missing_data()
  tree <- make_tree()

  s1 <- TreeLength(tree, ds, concavity = 10, extended_iw = TRUE,
                   xpiwe_r = 0.5)
  s2 <- TreeLength(tree, ds, concavity = 10, extended_iw = TRUE,
                   xpiwe_r = 1.0)

  # Both r=0.5 and r=1.0 produce f > 1 for chars with missing data, but
  # different magnitudes. Scores should differ when missing-data chars
  # have homoplasy.
  expect_false(isTRUE(all.equal(s1, s2)))

  # max_f clamping: set max_f=1 to force f clamped to 1, so eff_k = k always
  s_clamped <- TreeLength(tree, ds, concavity = 10, extended_iw = TRUE,
                          xpiwe_r = 1.0, xpiwe_max_f = 1)
  s_standard <- TreeLength(tree, ds, concavity = 10, extended_iw = FALSE)
  expect_equal(s_clamped, s_standard)
})

test_that("MaximizeParsimony accepts extended_iw", {
  data("inapplicable.datasets", package = "TreeSearch")
  ds <- inapplicable.phyData[["Vinther2008"]]
  result <- MaximizeParsimony(ds, concavity = 10, extended_iw = TRUE,
                              maxReplicates = 2L, targetHits = 1L,
                              verbosity = 0L)
  expect_true(inherits(result, "multiPhylo") || inherits(result, "phylo"))
})

test_that(".ObsCount returns correct counts", {
  ds <- make_missing_data()
  obsCount <- TreeSearch:::.ObsCount(ds)
  # 8 taxa. Chars 1-3: 8 observed. Chars 4-6: 6 observed (G, H missing).
  expect_true(all(obsCount >= 1))
  expect_true(all(obsCount <= 8))
  expect_true(any(obsCount < 8))
})

# --- TNT validation gallery ---
# Per-character concavity values from TNT 1.6 (piwe = 3, xpiwe(*0.5 <5).
# TNT's piwe& table shows 8 values per decade (chars 0-7, 10-17, ...),
# skipping chars ending in 8 or 9.
tnt_ref <- list()
tnt_ref[["Vinther2008"]] <- list(
  xpiwe_score = 4.04283,
  char_k = c(
    1.94118, 2.05714, 1.94118, 1.94118, 1.81818, 2.05714, 2.05714,
    1.81818, 2.16667, 1.94118, 2.05714, 2.05714, 2.79070, 2.05714,
    3.00000, 2.86364, 2.71429, 2.79070, 2.79070, 2.27027, 2.79070,
    2.79070, 2.79070, 2.79070, 2.36842, 2.36842, 2.36842, 2.46154,
    2.46154, 2.86364, 3.00000, 3.00000, 2.86364, 1.54839, 3.00000,
    3.00000, 1.81818, 2.86364, 2.93333, 3.00000, 2.05714, 3.00000,
    2.86364, 3.00000, 2.86364, 1.94118, 2.79070)
)
tnt_ref[["Sano2011"]] <- list(
  xpiwe_score = 20.51881,
  char_k = c(
    2.82353, 2.82353, 2.62500, 2.45902, 2.33898, 2.91429, 2.82353,
    2.91429, 2.91429, 2.95775, 2.57143, 2.72727, 2.51613, 2.67692,
    2.51613, 2.77612, 2.72727, 2.67692, 2.72727, 2.62500, 2.51613,
    2.77612, 2.72727, 3.00000, 2.45902, 2.45902, 2.40000, 2.82353,
    1.68000, 2.62500, 1.92453, 2.72727, 2.91429, 1.76471, 3.00000,
    2.95775, 2.95775, 2.95775, 2.95775, 2.91429, 1.59184, 2.86957)
)
tnt_ref[["Sansom2010"]] <- list(
  xpiwe_score = 16.24712,
  char_k = c(
    1.07143, 1.54839, 1.94118, 1.24138, 1.07143, 1.40000, 1.54839,
    1.07143, 1.07143, 1.07143, 1.81818, 1.07143, 1.81818, 1.40000,
    1.07143, 1.40000, 1.24138, 1.68750, 1.07143, 2.05714, 1.81818,
    2.05714, 2.55000, 2.71429, 2.79070, 2.27027, 1.07143, 1.81818,
    2.71429, 2.79070, 1.40000, 1.24138, 1.24138, 1.40000, 1.07143,
    1.68750, 2.55000, 2.71429, 2.71429, 2.16667, 2.93333, 2.36842,
    2.36842, 2.46154, 2.46154, 1.40000, 1.54839, 2.46154, 1.68750,
    1.40000, 1.24138, 1.94118, 3.00000, 1.07143, 1.40000, 1.40000,
    2.79070, 2.71429, 2.93333, 2.86364, 2.71429, 1.54839, 3.00000,
    2.79070, 1.24138, 1.54839, 2.79070, 2.86364, 2.93333, 1.81818,
    3.00000, 2.79070, 3.00000, 2.86364, 2.71429, 2.86364, 1.07143,
    1.07143, 1.07143, 1.07143, 1.54839, 1.24138, 1.07143, 0.88889,
    0.69231, 0.88889, 1.24138, 1.81818)
)

# Helper: map TNT 8-per-decade indices to 0-based character indices
tnt_char_indices <- function(n_char) {
  idx <- integer(0)
  decade <- 0L
  while (decade < n_char) {
    n_in_row <- min(8L, n_char - decade)
    idx <- c(idx, decade + seq_len(n_in_row))
    decade <- decade + 10L
  }
  idx
}

test_that("Per-character concavities match TNT 1.6", {
  data("inapplicable.phyData", package = "TreeSearch")

  for (ds_name in names(tnt_ref)) {
    ds <- inapplicable.phyData[[ds_name]]
    obsCount <- TreeSearch:::.ObsCount(ds)
    nTaxa <- length(ds)
    k <- 3; r <- 0.5; max_f <- 5

    f <- pmin(pmax(1 + r * (nTaxa - obsCount) / obsCount, 1), max_f)
    eff_k <- k / f

    # Expand from patterns to original characters
    eff_k_orig <- eff_k[attr(ds, "index")]
    n_char <- length(eff_k_orig)

    # Select the characters TNT reports (8 per decade)
    tnt_idx <- tnt_char_indices(n_char)
    our_k <- eff_k_orig[tnt_idx]
    ref_k <- tnt_ref[[ds_name]]$char_k

    expect_equal(length(our_k), length(ref_k),
                 info = paste(ds_name, "k-value count"))
    expect_equal(our_k, ref_k, tolerance = 1e-4,
                 info = paste(ds_name, "k-values vs TNT"))
  }
})

test_that(".ObsCount counts inapplicable (-) as missing", {
  # Dataset with both ? and -
  mat <- matrix(c(
    "0", "0", "0",
    "1", "1", "1",
    "0", "?", "-",
    "1", "?", "-"
  ), nrow = 4, byrow = TRUE,
  dimnames = list(LETTERS[1:4], paste0("c", 1:3)))
  ds <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1", "-"))

  obsCount <- TreeSearch:::.ObsCount(ds)
  # c1: 4 observed, c2: 2 observed (C,D are ?), c3: 2 observed (C,D are -)
  expect_true(all(obsCount[attr(ds, "index") == attr(ds, "index")[1]] == 4))
  # Both ? and - should reduce observed count
  expect_true(all(obsCount[attr(ds, "index") == attr(ds, "index")[2]] == 2))
  expect_true(all(obsCount[attr(ds, "index") == attr(ds, "index")[3]] == 2))
})
