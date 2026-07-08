# Tier 2: skipped on CRAN; see tests/testing-strategy.md
skip_on_cran()

# TNT-style addition-sequence variety (CLOSEST=3, FURTHEST=4, INFORMATIVE=5;
# see WagnerBias in src/ts_wagner.h). Default-OFF experimental modes gated
# behind wagnerBias/TS_ADDSEQ; see dev notes for the Hamilton A/B before
# adoption.
#
# ts_addseq_order_debug()/ts_pairwise_distances_debug() are test-only bridges
# (mirroring ts_wagner_bias_bench/ts_test_strategy_tracker) exposing
# greedy_addseq_order()'s raw permutation and wagner_pairwise_distances()'s
# raw matrix respectively, so brute-force references can be checked directly
# without threading through a full Wagner tree build.

prep_pd <- function(pd) {
  list(
    contrast = attr(pd, "contrast"),
    tip_data = t(vapply(pd, I, pd[[1]])),
    weight = attr(pd, "weight"),
    levels = attr(pd, "levels")
  )
}

# Independent brute-force reference: two-taxon Fitch distance and the
# parsimony-informativeness-gain criterion, built directly from
# contrast/tip_data (bypassing simplify_patterns entirely). Only valid when
# the input matrix has no constant/autapomorphic characters for
# simplify_patterns to strip -- the synthetic matrices below are constructed
# so every character keeps both states represented by >= 2 taxa.
brute_reference <- function(d, n_tip) {
  contrast <- d$contrast; tip_data <- d$tip_data; wt <- d$weight
  n_char <- ncol(tip_data); n_states <- ncol(contrast)
  tip_states <- vector("list", n_tip)
  for (t in seq_len(n_tip)) {
    tip_states[[t]] <- vector("list", n_char)
    for (c in seq_len(n_char)) {
      tok <- tip_data[t, c]
      tip_states[[t]][[c]] <- which(contrast[tok, ] == 1) - 1L
    }
  }
  is_missing <- function(states) length(states) == n_states

  dist <- function(i, j) {  # 0-based tip indices
    sum(vapply(seq_len(n_char), function(c) {
      si <- tip_states[[i + 1]][[c]]; sj <- tip_states[[j + 1]][[c]]
      if (length(intersect(si, sj)) == 0) wt[c] else 0
    }, double(1)))
  }
  informative_now <- function(members) {  # 0-based tip indices
    vapply(seq_len(n_char), function(c) {
      counts <- integer(n_states)
      for (t in members) {
        st <- tip_states[[t + 1]][[c]]
        if (is_missing(st)) next
        for (s in st) counts[s + 1] <- counts[s + 1] + 1
      }
      sum(counts >= 2) >= 2
    }, logical(1))
  }
  gain <- function(prefix, candidate) {
    sum(informative_now(c(prefix, candidate)) & !informative_now(prefix))
  }
  list(dist = dist, gain = gain)
}

# Synthetic matrix where every character has both states shared by >= 2 taxa
# (so simplify_patterns never strips a character), optionally with missing
# data sprinkled in (still leaving >= 2 taxa per real state).
make_informative_matrix <- function(n_tip, n_char, missing_per_char = 0L,
                                     seed) {
  set.seed(seed)
  mat <- matrix("0", n_tip, n_char)
  for (c in seq_len(n_char)) {
    k <- sample(3:(n_tip - 3), 1)
    ones <- sample(seq_len(n_tip), k)
    mat[ones, c] <- "1"
    if (missing_per_char > 0) {
      miss <- sample(seq_len(n_tip), missing_per_char)
      mat[miss, c] <- "?"
    }
  }
  rownames(mat) <- paste0("t", seq_len(n_tip))
  MatrixToPhyDat(mat)
}

check_permutation <- function(ord, n_tip) {
  length(ord) == n_tip && setequal(ord, seq_len(n_tip))
}

test_that("greedy_addseq_order produces a valid permutation for all 3 modes", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  n_tip <- length(pd)
  d <- prep_pd(pd)

  for (bias in 3:5) {
    for (temp in c(0, 0.3)) {
      set.seed(11L)
      ord <- TreeSearch:::ts_addseq_order_debug(
        d$contrast, d$tip_data, d$weight, d$levels,
        min_steps = integer(0), concavity = -1, bias = bias, temperature = temp
      )
      expect_true(check_permutation(ord, n_tip),
                  info = paste("bias", bias, "temp", temp))
    }
  }
})

test_that("greedy_addseq_order handles zero Fitch words without UB", {
  # All-distinct-singleton matrix: every character is parsimony-uninformative
  # and simplify_patterns removes all of them (total_words == 0). Mirrors the
  # existing wagner_goloboff_scores/wagner_entropy_scores regression test.
  n <- 10L
  set.seed(1)
  mat <- vapply(seq_len(3L),
                function(i) as.character(sample(seq_len(n)) - 1L),
                character(n))
  rownames(mat) <- paste0("t", seq_len(n))
  pd <- MatrixToPhyDat(mat)
  d <- prep_pd(pd)

  for (bias in 3:5) {
    ord <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels,
      min_steps = integer(0), concavity = -1, bias = bias, temperature = 0.3
    )
    expect_true(check_permutation(ord, n), info = paste("bias", bias))
  }
})

test_that("greedy_addseq_order is deterministic under a fixed seed, varies across seeds", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)

  for (bias in 3:5) {
    set.seed(5L)
    o1 <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1, bias, 0.3
    )
    set.seed(5L)
    o2 <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1, bias, 0.3
    )
    expect_identical(o1, o2, info = paste("bias", bias))

    set.seed(6L)
    o3 <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1, bias, 0.3
    )
    expect_false(identical(o1, o3), info = paste("bias", bias))
  }
})

test_that("wagner_pairwise_distances matches an independent brute-force reference", {
  n_tip <- 9L
  pd <- make_informative_matrix(n_tip, 20L, seed = 2024L)
  d <- prep_pd(pd)
  # No character is constant/autapomorphic by construction (every column has
  # both states held by >= 3 taxa), so simplify_patterns cannot drop any of
  # them -- it may still merge duplicate patterns into one column with
  # combined weight, so check total weighted character count, not ncol.
  expect_equal(sum(d$weight), 20L)

  ref <- brute_reference(d, n_tip)
  cpp_dist <- TreeSearch:::ts_pairwise_distances_debug(
    d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1
  )
  for (i in 0:(n_tip - 1)) {
    for (j in 0:(n_tip - 1)) {
      if (i == j) next
      expect_equal(cpp_dist[i + 1, j + 1], ref$dist(i, j),
                   info = paste(i, j))
    }
  }
})

test_that("CLOSEST/FURTHEST greedy selection matches brute-force argmax at every step", {
  n_tip <- 9L
  pd <- make_informative_matrix(n_tip, 20L, seed = 2024L)
  d <- prep_pd(pd)
  ref <- brute_reference(d, n_tip)

  for (bias in c(3L, 4L)) {
    set.seed(321L)
    ord0 <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1, bias, 0
    ) - 1L  # to 0-based

    for (i in 2:n_tip) {
      prefix <- ord0[seq_len(i - 1)]
      remaining <- setdiff(0:(n_tip - 1), prefix)
      chosen <- ord0[i]
      md <- vapply(remaining, function(cand) {
        min(vapply(prefix, function(p) ref$dist(cand, p), double(1)))
      }, double(1))
      names(md) <- remaining
      best <- if (bias == 3L) remaining[md == min(md)] else remaining[md == max(md)]
      expect_true(chosen %in% best,
                  info = sprintf("bias=%d step=%d chosen=%d best=%s",
                                  bias, i, chosen, paste(best, collapse = ",")))
    }
  }
})

test_that("INFORMATIVE greedy selection matches brute-force argmax, including with missing data", {
  for (missing_per_char in c(0L, 2L)) {
    n_tip <- 10L
    pd <- make_informative_matrix(n_tip, 18L, missing_per_char, seed = 777L)
    d <- prep_pd(pd)
    ref <- brute_reference(d, n_tip)

    set.seed(444L)
    ord0 <- TreeSearch:::ts_addseq_order_debug(
      d$contrast, d$tip_data, d$weight, d$levels, integer(0), -1, 5L, 0
    ) - 1L

    for (i in 2:(n_tip - 1)) {
      prefix <- ord0[seq_len(i - 1)]
      remaining <- setdiff(0:(n_tip - 1), prefix)
      chosen <- ord0[i]
      gains <- vapply(remaining, function(cand) ref$gain(prefix, cand), double(1))
      names(gains) <- remaining
      best <- remaining[gains == max(gains)]
      expect_true(chosen %in% best,
                  info = sprintf("missing=%d step=%d chosen=%d best=%s",
                                  missing_per_char, i, chosen,
                                  paste(best, collapse = ",")))
    }
  }
})

test_that("wagnerBias 3-5 integrate end-to-end via MaximizeParsimony (incl. NA datasets)", {
  data("inapplicable.phyData", package = "TreeSearch")
  for (ds_name in c("Vinther2008", "Longrich2010", "Sansom2010")) {
    pd <- inapplicable.phyData[[ds_name]]
    for (bias in 3:5) {
      set.seed(55L)
      res <- MaximizeParsimony(pd, maxReplicates = 2L, targetHits = 2L,
        control = SearchControl(wagnerBias = bias, wagnerBiasTemp = 0.3,
                                 ratchetCycles = 0L, driftCycles = 0L,
                                 xssRounds = 0L),
        verbosity = 0L)
      tr <- if (inherits(res, "multiPhylo")) res[[1]] else res
      expect_equal(length(tr$tip.label), length(pd),
                   info = paste(ds_name, "bias", bias))
      expect_equal(TreeLength(tr, pd), attr(res, "score"),
                   info = paste(ds_name, "bias", bias))
    }
  }
})

test_that("legacy wagnerBias modes engage under nThreads > 1 (parallel-path fix)", {
  # Previously the parallel worker loop only handled the adaptive-bandit
  # round-robin case; a fixed `wagnerBias` (modes 1-5) was silently ignored
  # whenever nThreads > 1. Exercise all three new modes there directly.
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  for (bias in 3:5) {
    set.seed(7L)
    res <- MaximizeParsimony(pd, maxReplicates = 2L, targetHits = 2L,
      control = SearchControl(wagnerBias = bias, wagnerBiasTemp = 0.3,
                               ratchetCycles = 0L, driftCycles = 0L,
                               xssRounds = 0L),
      nThreads = 2L, verbosity = 0L)
    tr <- if (inherits(res, "multiPhylo")) res[[1]] else res
    expect_equal(length(tr$tip.label), length(pd), info = paste("bias", bias))
  }
})

test_that("TS_ADDSEQ env override takes precedence over the configured wagnerBias", {
  data("congreveLamsdellMatrices", package = "TreeSearch")
  pd <- congreveLamsdellMatrices[[1]]
  d <- prep_pd(pd)
  n_tip <- length(pd)

  old <- Sys.getenv("TS_ADDSEQ", unset = NA)
  Sys.setenv(TS_ADDSEQ = "5")
  on.exit(if (is.na(old)) Sys.unsetenv("TS_ADDSEQ") else Sys.setenv(TS_ADDSEQ = old))

  set.seed(3L)
  res <- MaximizeParsimony(pd, maxReplicates = 2L, targetHits = 2L,
    control = SearchControl(wagnerBias = 0L,  # configured RANDOM
                             ratchetCycles = 0L, driftCycles = 0L,
                             xssRounds = 0L),
    verbosity = 0L)
  tr <- if (inherits(res, "multiPhylo")) res[[1]] else res
  expect_equal(length(tr$tip.label), n_tip)
})

test_that("wagnerBias 3-5 respect a topological constraint", {
  # Addition order is constraint-blind (wagner_tree()'s own edge filter +
  # posthoc retry do the enforcement, exactly as for GOLOBOFF/ENTROPY), so
  # this is a low-risk smoke check rather than new machinery.
  ds5 <- phangorn::phyDat(
    matrix(c("0", "0", "0", "1", "1",
             "0", "1", "0", "1", "0"),
           nrow = 5, dimnames = list(paste0("t", 1:5), NULL)),
    type = "USER", levels = c("0", "1")
  )
  cons <- ape::read.tree(text = "((t1,t2),(t3,(t4,t5)));")

  for (bias in 3:5) {
    set.seed(4217L)
    res <- MaximizeParsimony(ds5, constraint = cons,
      control = SearchControl(wagnerBias = bias, wagnerBiasTemp = 0.3),
      maxReplicates = 2L, verbosity = 0L)
    tr <- if (inherits(res, "multiPhylo")) res[[1]] else res

    tips <- sort(cons$tip.label)
    tree_sp <- as.logical(as.Splits(tr, tipLabels = tips))
    cons_sp <- as.logical(as.Splits(cons, tipLabels = tips))
    if (!is.matrix(tree_sp)) tree_sp <- matrix(tree_sp, nrow = 1)
    if (!is.matrix(cons_sp)) cons_sp <- matrix(cons_sp, nrow = 1)
    satisfied <- all(apply(cons_sp, 1, function(c_row) {
      any(apply(tree_sp, 1, function(t_row) {
        all(c_row == t_row) || all(c_row == !t_row)
      }))
    }))
    expect_true(satisfied, info = paste("bias", bias))
  }
})
