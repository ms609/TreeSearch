library("TreeTools")

test_that("Split count is T-3 for resolved tree", {
  for (n_tip in c(5, 8, 10, 20)) {
    tree <- as.phylo(1, n_tip)
    splits <- ts_compute_splits(tree$edge, n_tip)
    expect_equal(length(splits), n_tip - 3,
                 label = paste0("n_tip=", n_tip))
  }
})

test_that("Same topology rooted differently produces same splits", {
  tree1 <- as.phylo(42, 8)
  tree2 <- TreeTools::RootTree(tree1, 3)
  tree2 <- TreeTools::Preorder(tree2)

  expect_true(
    ts_trees_equal(tree1$edge, tree2$edge, 8)
  )
})

test_that("Known splits for a small tree", {
  # 5-tip tree: should have exactly 2 non-trivial splits
  tree <- as.phylo(1, 5)
  splits <- ts_compute_splits(tree$edge, 5)
  expect_equal(length(splits), 2)

  # Each split should contain at least 2 and at most 3 tips
  for (s in splits) {
    expect_true(length(s) >= 2 && length(s) <= 3)
  }
})

test_that("Identical trees have same hash / splits_equal = TRUE", {
  tree <- as.phylo(42, 10)
  expect_true(ts_trees_equal(tree$edge, tree$edge, 10))
})

test_that("Different topologies give splits_equal = FALSE", {
  tree1 <- as.phylo(1, 10)
  tree2 <- as.phylo(2, 10)
  expect_false(ts_trees_equal(tree1$edge, tree2$edge, 10))
})

test_that("NNI can produce a different topology", {
  tree <- TreeTools::Preorder(as.phylo(42, 10))
  n_tip <- 10L

  # Try NNI on each internal edge until we find one that changes the topology
  internal_edges <- which(tree$edge[, 2] > n_tip)
  found_different <- FALSE
  for (ie in internal_edges) {
    tree2 <- NNI(tree, ie)
    tree2 <- TreeTools::Preorder(tree2)
    if (!ts_trees_equal(tree$edge, tree2$edge, n_tip)) {
      found_different <- TRUE
      # Both should have the same number of splits
      s1 <- ts_compute_splits(tree$edge, n_tip)
      s2 <- ts_compute_splits(tree2$edge, n_tip)
      expect_equal(length(s1), length(s2))
      break
    }
  }
  expect_true(found_different, info = "At least one NNI should change topology")
})

test_that("Various tree sizes work", {
  for (n_tip in c(5, 10, 20, 50)) {
    tree <- as.phylo(1, n_tip)
    splits <- ts_compute_splits(tree$edge, n_tip)
    expect_equal(length(splits), n_tip - 3,
                 label = paste0("n_tip=", n_tip))
  }
})

test_that("Multi-word splits work (65+ tips)", {
  skip_if_not_installed("TreeTools", minimum_version = "1.9.0")

  for (n_tip in c(65, 100)) {
    tree <- TreeTools::RandomTree(n_tip, root = TRUE)
    tree <- TreeTools::Preorder(tree)
    splits <- ts_compute_splits(tree$edge, n_tip)
    expect_equal(length(splits), n_tip - 3,
                 label = paste0("n_tip=", n_tip, " multi-word"))

    # Same tree re-rooted should be equal
    tree2 <- TreeTools::RootTree(tree, 3)
    tree2 <- TreeTools::Preorder(tree2)
    expect_true(ts_trees_equal(tree$edge, tree2$edge, n_tip),
                label = paste0("n_tip=", n_tip, " reroot equality"))
  }
})

test_that("Splits are canonical (tip 0 always in 0 partition)", {

  tree <- as.phylo(42, 10)
  splits <- ts_compute_splits(tree$edge, 10)
  # Tip 1 (R 1-based = C++ tip 0) should NOT be in any split

  # (because canonical form ensures tip 0 is in the "0" = unset partition)
  for (s in splits) {
    expect_false(1 %in% s,
                 info = "Tip 1 (= C++ tip 0) should not be in canonical split")
  }
})

test_that("Cross-validate with TreeTools::as.Splits", {
  skip_if_not_installed("TreeTools")

  tree <- as.phylo(42, 10)
  our_splits <- ts_compute_splits(tree$edge, 10)

  # TreeTools gives splits as logical matrix rows
  tt_splits <- as.Splits(tree)
  tt_mat <- as.logical(tt_splits)
  if (!is.matrix(tt_mat)) tt_mat <- matrix(tt_mat, nrow = 1)

  # Both should have n_tip - 3 splits

  expect_equal(length(our_splits), nrow(tt_mat))

  # Convert our splits to logical vectors for comparison
  our_logical <- lapply(our_splits, function(tips) {
    v <- rep(FALSE, 10)
    v[tips] <- TRUE
    v
  })

  # For each of our splits, it should match one TreeTools split or its complement
  for (our_s in our_logical) {
    complement <- !our_s
    found <- FALSE
    for (i in seq_len(nrow(tt_mat))) {
      if (identical(as.logical(our_s), as.logical(tt_mat[i, ])) ||
          identical(as.logical(complement), as.logical(tt_mat[i, ]))) {
        found <- TRUE
        break
      }
    }
    expect_true(found, info = "Each C++ split should match a TreeTools split")
  }
})
