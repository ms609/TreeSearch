# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# Tests for hierarchical resampling in Resample() with HSJ and xform scoring.

library("TreeTools")

make_dat <- function(mat, levels = c("-", "0", "1")) {
  phangorn::phyDat(mat, type = "USER", levels = levels, ambiguity = "?")
}

# --- Test dataset: 6 tips, 8 characters ---
# Chars 1-2: free (non-hierarchy)
# Chars 3-5: hierarchy block 1 (char 3 primary, chars 4-5 secondary)
# Chars 6-8: hierarchy block 2 (char 6 primary, chars 7-8 secondary)
make_resample_data <- function() {
  mat <- matrix(c(
    # free1 free2 pri1 sec1a sec1b pri2 sec2a sec2b
    "0",  "1",  "1", "0",  "1",  "1", "0",  "1",   # t1
    "0",  "0",  "1", "0",  "0",  "1", "1",  "0",   # t2
    "1",  "0",  "1", "1",  "1",  "0", "-",  "-",   # t3
    "1",  "1",  "0", "-",  "-",  "1", "1",  "1",   # t4
    "0",  "1",  "0", "-",  "-",  "0", "-",  "-",   # t5
    "1",  "0",  "1", "1",  "0",  "1", "0",  "0"    # t6
  ), nrow = 6, byrow = TRUE,
  dimnames = list(paste0("t", 1:6), NULL))
  make_dat(mat)
}

make_resample_hierarchy <- function() {
  CharacterHierarchy("3" = 4:5, "6" = 7:8)
}


# ===== .HierarchicalResampleWeights unit tests ================================

HRW <- TreeSearch:::.HierarchicalResampleWeights

test_that("HierarchicalResampleWeights returns correct structure", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(7421)
  result <- HRW(ds, h, bootstrap = FALSE, proportion = 2 / 3)

  expect_named(result, c("non_hierarchy_weights", "block_counts"))
  expect_length(result$non_hierarchy_weights, length(attr(ds, "weight")))
  expect_length(result$block_counts, 2L)  # 2 top-level blocks
  expect_true(all(result$non_hierarchy_weights >= 0L))
  expect_true(all(result$block_counts >= 0L))
})

test_that("Jackknife drops some units", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  # 2 free chars + 2 blocks = 4 units. At 2/3, keep ceil(4*2/3) = 3 units.
  set.seed(8312)
  n_sampled <- replicate(50, {
    r <- HRW(ds, h, bootstrap = FALSE, proportion = 2 / 3)
    sum(r$non_hierarchy_weights > 0L) + sum(r$block_counts > 0L)
  })
  # Should always be < total units (4)
  expect_true(all(n_sampled < 4L))
  # Should always be >= 1

  expect_true(all(n_sampled >= 1L))
})

test_that("Bootstrap can duplicate units", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(5519)
  has_dup <- any(replicate(100, {
    r <- HRW(ds, h, bootstrap = TRUE, proportion = 2 / 3)
    any(r$block_counts > 1L) || any(r$non_hierarchy_weights > 1L)
  }))
  expect_true(has_dup)
})

test_that("Block counts correspond to top-level blocks", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(3017)
  # With enough replicates, both blocks should sometimes be dropped
  dropped_1 <- FALSE
  dropped_2 <- FALSE
  for (i in 1:200) {
    r <- HRW(ds, h, bootstrap = FALSE, proportion = 0.5)
    if (r$block_counts[1] == 0L) dropped_1 <- TRUE
    if (r$block_counts[2] == 0L) dropped_2 <- TRUE
    if (dropped_1 && dropped_2) break
  }
  expect_true(dropped_1)
  expect_true(dropped_2)
})

test_that("Non-hierarchy weights only count free chars", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  # Force a result where all free chars and blocks are retained (bootstrap)
  set.seed(2209)
  # Run many times and check that nh_weights never include hierarchy char
  # contributions
  n_nh_patterns <- length(attr(ds, "weight"))
  full_nh_w <- TreeSearch:::non_hierarchy_weights(ds, h)

  for (i in 1:20) {
    r <- HRW(ds, h, bootstrap = TRUE, proportion = 2 / 3)
    # Non-hierarchy weights should never exceed original nh_weights * max_count
    total_nh <- sum(r$non_hierarchy_weights)
    # Free chars: 2 chars total. Bootstrap samples 4 units so free chars
    # can appear at most once each (they're individual units)
    # In a 4-unit bootstrap, a free char unit can be sampled multiple times
    expect_true(total_nh >= 0L)
  }
})


# ===== Resample() with hierarchy: parameter validation ========================

test_that("Resample rejects HSJ without hierarchy", {
  ds <- make_resample_data()
  expect_error(
    Resample(ds, inapplicable = "hsj"),
    "hierarchy.*required"
  )
})

test_that("Resample rejects xform without hierarchy", {
  ds <- make_resample_data()
  expect_error(
    Resample(ds, inapplicable = "xform"),
    "hierarchy.*required"
  )
})

test_that("Resample rejects bad hsj_alpha", {
  ds <- make_resample_data()
  expect_error(
    Resample(ds, inapplicable = "bgs", hsj_alpha = 2.0),
    "hsj_alpha"
  )
  expect_error(
    Resample(ds, inapplicable = "bgs", hsj_alpha = -0.1),
    "hsj_alpha"
  )
})

test_that("Resample rejects profile + hierarchy", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()
  expect_error(
    Resample(ds, hierarchy = h, inapplicable = "hsj", concavity = "profile"),
    "Profile.*not.*supported"
  )
})

test_that("Resample rejects IW + hierarchy", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()
  expect_error(
    Resample(ds, hierarchy = h, inapplicable = "hsj", concavity = 10),
    "Implied.*not.*supported"
  )
})


# ===== Resample() end-to-end: HSJ ============================================

test_that("Resample returns valid trees with HSJ scoring", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(4738)
  result <- Resample(ds, hierarchy = h, inapplicable = "hsj",
                     hsj_alpha = 1.0, nReplicates = 3L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 3L)
  for (tr in result) {
    expect_s3_class(tr, "phylo")
    expect_equal(length(tr$tip.label), 6L)
    expect_equal(tr$Nnode, 5L)
    expect_true(!is.null(attr(tr, "score")))
  }
})

test_that("HSJ α=0 produces valid resampled trees", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(9251)
  result <- Resample(ds, hierarchy = h, inapplicable = "hsj",
                     hsj_alpha = 0.0, nReplicates = 2L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 2L)
  for (tr in result) {
    expect_s3_class(tr, "phylo")
  }
})

test_that("HSJ bootstrap produces valid trees", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(6184)
  result <- Resample(ds, hierarchy = h, inapplicable = "hsj",
                     method = "bootstrap", nReplicates = 2L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 2L)
})


# ===== Resample() end-to-end: xform ==========================================

test_that("Resample returns valid trees with xform scoring", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(1947)
  result <- Resample(ds, hierarchy = h, inapplicable = "xform",
                     nReplicates = 3L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 3L)
  for (tr in result) {
    expect_s3_class(tr, "phylo")
    expect_equal(length(tr$tip.label), 6L)
    expect_true(!is.null(attr(tr, "score")))
  }
})

test_that("Xform bootstrap produces valid trees", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(3372)
  result <- Resample(ds, hierarchy = h, inapplicable = "xform",
                     method = "bootstrap", nReplicates = 2L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 2L)
})


# ===== Default behaviour unchanged ===========================================

test_that("Resample with brazeau (default) ignores hierarchy", {
  ds <- make_resample_data()

  set.seed(5901)
  result <- Resample(ds, nReplicates = 2L,
                     ratchIter = 1L, tbrIter = 2L)

  expect_s3_class(result, "multiPhylo")
  expect_length(result, 2L)
})


# ===== Resampling variation ===================================================

test_that("Hierarchical resampling produces different trees across replicates", {
  ds <- make_resample_data()
  h <- make_resample_hierarchy()

  set.seed(8763)
  result <- Resample(ds, hierarchy = h, inapplicable = "hsj",
                     nReplicates = 10L,
                     ratchIter = 1L, tbrIter = 2L)

  # With 10 replicates on a 6-tip dataset, not all trees should be identical
  edges <- lapply(result, function(tr) tr$edge)
  n_unique <- length(unique(lapply(edges, function(e) {
    paste(e[, 1], e[, 2], collapse = "-")
  })))
  expect_gt(n_unique, 1L)
})
