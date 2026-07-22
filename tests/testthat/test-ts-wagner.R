# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()
# Helper to prepare phyDat for C++ engine
prep_pd <- function(pd) {
  list(
    contrast = attr(pd, "contrast"),
    tip_data = t(vapply(pd, I, pd[[1]])),
    weight = attr(pd, "weight"),
    levels = attr(pd, "levels")
  )
}

test_that("ts_wagner_tree produces valid tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_true(is.list(result))
  expect_equal(ncol(result$edge), 2L)
  expect_equal(nrow(result$edge), 2L * (n_tip - 1L))
  expect_true(result$score > 0)
})

test_that("Wagner score matches ts_fitch_score", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(result$score, fitch_check)
})

test_that("Wagner score matches TreeLength", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  result_tree <- structure(list(
    edge = result$edge,
    tip.label = names(pd),
    Nnode = n_tip - 1L
  ), class = "phylo")

  expect_equal(result$score, TreeLength(result_tree, pd))
})

test_that("All tips present exactly once", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  child_nodes <- result$edge[, 2]
  tips_found <- sort(child_nodes[child_nodes <= n_tip])
  expect_equal(tips_found, seq_len(n_tip))
})

test_that("Same addition order gives same tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  order_seq <- seq_len(n_tip)
  r1 <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels, order_seq
  )
  r2 <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels, order_seq
  )

  expect_identical(r1$edge, r2$edge)
  expect_identical(r1$score, r2$score)
})

test_that("Random Wagner trees vary", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(7263)
  r1 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  set.seed(1984)
  r2 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_false(identical(r1$edge, r2$edge))
})

test_that("Random Wagner score verified", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(5511)
  result <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(result$score, fitch_check)
})

test_that("Small tree (5 tips) is correct", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd5 <- congreveLamsdellMatrices[[1]][1:5]
  d <- prep_pd(pd5)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(nrow(result$edge), 8L)
  expect_true(result$score > 0)

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)
})

test_that("Medium tree (20 tips) completes without error", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd20 <- congreveLamsdellMatrices[[1]][1:20]
  d <- prep_pd(pd20)

  expect_no_error({
    result <- TreeSearch:::ts_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )
  })

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)
})

test_that("Multiple datasets produce verified scores", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  for (i in 1:3) {
    pd <- congreveLamsdellMatrices[[i]]
    d <- prep_pd(pd)

    set.seed(3000 + i)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check,
                 info = paste("Dataset", i))
  }
})

# --- New tests for incremental scoring correctness ---

test_that("Wagner on inapplicable datasets matches fitch_score", {
  data("inapplicable.phyData", package = "TreeSearch")
  na_datasets <- c("Vinther2008", "Longrich2010", "Sansom2010")
  for (ds_name in na_datasets) {
    pd <- inapplicable.phyData[[ds_name]]
    d <- make_ts_data(pd)

    set.seed(4217)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check, info = ds_name)
  }
})

test_that("Wagner on NA datasets is deterministic", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)

  set.seed(6193)
  r1 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  set.seed(6193)
  r2 <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(r1$score, r2$score)
  expect_equal(r1$edge, r2$edge)
})

test_that("Wagner on NA + IW matches fitch_score", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)

  # T-322: pass min_steps exactly as MaximizeParsimony()'s IW path does, so the
  # cross-check exercises the production formula h = steps - min_steps (not the
  # degenerate h = steps - 0). Vinther2008 has inapplicable chars => non-zero.
  min_steps <- as.integer(MinimumLength(pd, compress = TRUE))

  for (k in c(3, 10)) {
    set.seed(8514)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels,
      min_steps = min_steps, concavity = k
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels,
      min_steps = min_steps, concavity = k
    )
    expect_equal(result$score, fitch_check, tolerance = 1e-6,
                 info = paste("IW k =", k))
  }
})

test_that("Wagner NA tree has valid topology", {
  data("inapplicable.phyData", package = "TreeSearch")
  pd <- inapplicable.phyData[["Vinther2008"]]
  d <- make_ts_data(pd)
  n_tip <- length(pd)

  set.seed(2917)
  result <- TreeSearch:::ts_random_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  expect_equal(nrow(result$edge), 2L * (n_tip - 1L))
  tips <- sort(result$edge[result$edge[, 2] <= n_tip, 2])
  expect_equal(tips, seq_len(n_tip))
})

test_that("Wagner with many addition orders all verify", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  for (s in c(1234, 5678, 9012)) {
    set.seed(s)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check, info = paste("seed", s))
  }
})

test_that("Wagner minimum case: 3 tips produces valid tree", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd3 <- congreveLamsdellMatrices[[1]][1:3]
  d <- prep_pd(pd3)

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels
  )

  # 3 tips: 4 edges, 2 internal nodes
  expect_equal(nrow(result$edge), 4L)
  expect_true(result$score > 0)

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)

  # All 3 tips present exactly once
  child_nodes <- result$edge[, 2]
  tips_found <- sort(child_nodes[child_nodes <= 3])
  expect_equal(tips_found, 1:3)
})

test_that("Driven search still finds good scores after Wagner optimization", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  set.seed(8822)
  result <- TreeSearch:::ts_driven_search(
    d$contrast, d$tip_data, d$weight, d$levels,
    maxReplicates = 2L, targetHits = 1L,
    tbrMaxHits = 1L, ratchetCycles = 2L,
    ratchetPerturbProb = 0.04, driftCycles = 0L,
    driftAfdLimit = 3L, driftRfdLimit = 0.1,
    xssRounds = 0L, xssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 3L, fuseAcceptEqual = FALSE,
    poolMaxSize = 10L, poolSuboptimal = 0.0,
    maxSeconds = 0, verbosity = 0L
  )

  expect_true(result$best_score <= 200)
  expect_true(result$best_score > 0)
})

# ---- Constrained Wagner tests (S-RED focus 9) ----
# Regression for boundary-edge false positive (S-RED round 9):
# wagner_edge_violates_constraint and regraft_violates_constraint both
# rejected the edge (parent_of_cn, cn) for outside tips/clades.  Inserting
# an outside element just above the constraint clade makes it a sibling of
# that clade and does NOT break monophyly.  Fixed with `&& below != cn`.

test_that("constrained random Wagner score is verified", {
  # Score check only: ts_random_wagner_tree without posthoc data uses the
  # primary wagner_edge_violates_constraint check but the retry loop only fires
  # when has_posthoc=TRUE.  For addition orders where inside tips land on both
  # sides of the rooted root in the 3-taxon start, cn==root and the check is
  # skipped (by design — can't enforce directionality from an unrooted root).
  # Full constraint satisfaction requires posthoc data; that path is exercised
  # via MaximizeParsimony in test-ts-constraint-small.R.
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]][1:6]
  d <- prep_pd(pd)

  cons_mat <- matrix(c(0L, 0L, 1L, 1L, 1L, 1L), nrow = 1L)

  for (seed in c(1771L, 2882L, 3993L)) {
    set.seed(seed)
    result <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels,
      consSplitMatrix = cons_mat
    )
    fitch_check <- TreeSearch:::ts_fitch_score(
      result$edge, d$contrast, d$tip_data, d$weight, d$levels
    )
    expect_equal(result$score, fitch_check,
                 info = paste("constrained Wagner seed", seed))
  }
})

test_that("constrained sequential Wagner boundary edge: outside tip adjacent to clade", {
  # Sequential addition order (t1 first, the two inside tips go in first).
  # After placing t1, t2, t3 as 3-taxon tree, cn = LCA(t1,t2).
  # Adding t4 (outside): with the boundary fix, edge (root, cn) is now
  # accepted (previously false-positive rejected).  Either way, the
  # constraint must be satisfied and the final score correct.
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]][1:6]
  d <- prep_pd(pd)
  n_tip <- length(pd)

  cons_mat <- matrix(c(0L, 0L, 1L, 1L, 1L, 1L), nrow = 1L)
  order_seq <- seq_len(n_tip)  # R 1-indexed; t1 added first

  result <- TreeSearch:::ts_wagner_tree(
    d$contrast, d$tip_data, d$weight, d$levels,
    addition_order = order_seq,
    consSplitMatrix = cons_mat
  )

  fitch_check <- TreeSearch:::ts_fitch_score(
    result$edge, d$contrast, d$tip_data, d$weight, d$levels
  )
  expect_equal(result$score, fitch_check)

  # R tips 1 and 2 must be sisters
  ec <- result$edge
  p1 <- ec[ec[, 2] == 1L, 1L]
  p2 <- ec[ec[, 2] == 2L, 1L]
  expect_equal(p1, p2)

  # All tips present exactly once
  child_tips <- sort(result$edge[result$edge[, 2] <= n_tip, 2L])
  expect_equal(child_tips, seq_len(n_tip))
})

test_that("tip_data values outside [1, n_tokens] error cleanly (T-328)", {
  # Rcpp-boundary guard: an out-of-range tip_data entry must not reach the
  # unguarded 1-based index into token_states (src/ts_simplify.cpp) -- it
  # should error, not silently misread or crash.
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)
  n_tokens <- nrow(d$contrast)

  bad_low <- d$tip_data
  bad_low[1, 1] <- 0L
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, bad_low, d$weight, d$levels),
    "tip_data"
  )

  bad_high <- d$tip_data
  bad_high[1, 1] <- n_tokens + 1L
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, bad_high, d$weight, d$levels),
    "tip_data"
  )

  bad_negative <- d$tip_data
  bad_negative[1, 1] <- -1L
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, bad_negative, d$weight, d$levels),
    "tip_data"
  )
})

test_that("out-of-range tip_data errors via other Wagner entry points (T-328)", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)
  n_tokens <- nrow(d$contrast)
  bad <- d$tip_data
  bad[1, 1] <- n_tokens + 1L

  expect_error(
    TreeSearch:::ts_random_wagner_tree(d$contrast, bad, d$weight, d$levels),
    "tip_data"
  )
  expect_error(
    TreeSearch:::ts_resample_search(d$contrast, bad, d$weight, d$levels),
    "tip_data"
  )
  expect_error(
    TreeSearch:::ts_parallel_resample(d$contrast, bad, d$weight, d$levels),
    "tip_data"
  )
  expect_error(
    TreeSearch:::ts_successive_approx(d$contrast, bad, d$weight, d$levels),
    "tip_data"
  )
  expect_error(
    TreeSearch:::ts_simplify_diag(d$contrast, bad, d$weight, d$levels),
    "tip_data"
  )
})

test_that("ts_simplify_diag validates weight/levels length (T-328 hardening)", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  expect_error(
    TreeSearch:::ts_simplify_diag(d$contrast, d$tip_data, d$weight[-1],
                                   d$levels),
    "weight"
  )
  expect_error(
    TreeSearch:::ts_simplify_diag(d$contrast, d$tip_data, d$weight,
                                   d$levels[-1]),
    "levels"
  )
})

test_that("addition_order length/range/duplicate errors cleanly (T-323)", {
  # Rcpp-boundary guard: wagner_tree() reads addition_order as a length-n_tip
  # permutation with no validation -- a short vector reads past its end
  # (segfault), an out-of-range or duplicated value corrupts the tree.
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)
  n_tip <- length(pd)

  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels,
                                 addition_order = 1L),
    "addition_order"
  )

  bad_range <- seq_len(n_tip)
  bad_range[1] <- n_tip + 1L
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels,
                                 addition_order = bad_range),
    "addition_order"
  )

  bad_zero <- seq_len(n_tip)
  bad_zero[1] <- 0L
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels,
                                 addition_order = bad_zero),
    "addition_order"
  )

  dup <- seq_len(n_tip)
  dup[2] <- dup[1]
  expect_error(
    TreeSearch:::ts_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels,
                                 addition_order = dup),
    "addition_order"
  )
})

# ===== All-uninformative data: zero Fitch words (regression) ================
# When every character is constant or an autapomorphy (no state shared by
# more than one taxon), simplify_patterns removes every character, leaving
# DataSet::total_words == 0 and empty per-word state vectors.
# wagner_goloboff_scores()/wagner_entropy_scores() took the address of
# element 0 of the (empty) ds.tip_states vector -- undefined behaviour that
# aborts under hardened libstdc++ assertions / ASan.

test_that("wagner_goloboff_scores/wagner_entropy_scores handle zero Fitch words", {
  n <- 10L
  set.seed(1)
  # 3 characters, each a distinct permutation of the single-character states
  # 0..(n-1): no state is shared by more than one taxon in any character, and
  # every token is unambiguous, so every character is parsimony-uninformative
  # and simplify_patterns removes all of them. (States 0..9 are used, NOT
  # 1..10 -- the token "10" would be parsed as the AMBIGUITY {0,1}, which makes
  # the character genuinely ambiguous rather than all-distinct-singletons.)
  # Multi-character avoids the single-character vapply/t() degenerate case in
  # prep_pd().
  mat <- vapply(seq_len(3L),
                function(i) as.character(sample(seq_len(n)) - 1L),
                character(n))
  rownames(mat) <- paste0("t", seq_len(n))
  pd <- MatrixToPhyDat(mat)
  d <- prep_pd(pd)

  result <- TreeSearch:::ts_wagner_bias_bench(
    d$contrast, d$tip_data, d$weight, d$levels,
    min_steps = integer(0), concavity = -1,
    bias = 1L, temperature = 1, n_reps = 1L, run_tbr = FALSE
  )

  expect_true(all(result$goloboff_scores == 0))
  expect_true(all(result$entropy_scores == 0))
})
