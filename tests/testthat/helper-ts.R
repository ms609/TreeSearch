# Shared helpers for ts-* test files.
# testthat auto-sources helper-*.R files before running tests.

#' Set environment variables for the duration of `code`, restoring (or
#' unsetting) the prior values on exit. Minimal single-purpose stand-in for
#' withr::with_envvar so the package doesn't need withr just for this.
with_envvar <- function(vars, code) {
  old <- Sys.getenv(names(vars), unset = NA, names = TRUE)
  on.exit({
    keep <- !is.na(old)
    if (any(keep)) do.call(Sys.setenv, as.list(old[keep]))
    if (any(!keep)) Sys.unsetenv(names(old)[!keep])
  })
  do.call(Sys.setenv, as.list(vars))
  code
}

#' Skip a test unless TREESEARCH_EXTENDED_TESTS=true is set.
#' Use inside test_that() or at file level for Tier 3 (stress/bench) tests.
#' See tests/testing-strategy.md for full tiering documentation.
skip_extended <- function() {
  testthat::skip_if(
    !identical(Sys.getenv("TREESEARCH_EXTENDED_TESTS"), "true"),
    "Set TREESEARCH_EXTENDED_TESTS=true to run extended tests"
  )
}

library("TreeTools")

suppressMessages(suppressPackageStartupMessages(
  requireNamespace("phangorn", quietly = TRUE)
))

#' Deterministic LCG (Numerical Recipes constants). Double arithmetic keeps
#' every intermediate exact (well within 2^53) and sidesteps the 32-bit
#' integer overflow a literal integer LCG would hit.
.lcg_next <- function(state) (1664525 * state + 1013904223) %% 2147483647

#' A fixed, non-pectinate tree for a dataset, built deterministically by
#' successive random-edge tip insertion (the same construction RandomTree()
#' uses), but driven by a self-contained LCG rather than R's set.seed() +
#' sample(). This keeps the topology stable across R versions and
#' platforms, since it never touches R's RNG stream — unlike
#' set.seed()+RandomTree(), whose output would silently change if R's
#' sampling algorithm ever changes.
#'
#' Attempting to reuse TreeTools::as.phylo.numeric() (TreeNumber) for this
#' was tried and rejected: small index values map to a near-pectinate
#' topology for these tip counts (Colless index ~1.0, i.e. no more mixed
#' than PectinateTree()), because the mixed-radix place values are
#' dominated by the first few leaves for any index much smaller than
#' NUnrooted(nTip). This LCG builder instead reliably lands at Colless
#' index ~0.2-0.4 (vs. 1.0 for pectinate, ~0 for a fully balanced tree)
#' across the datasets these tests use.
FixedTree <- function(dataset, index) {
  tipLabels <- names(dataset)
  n <- length(tipLabels)
  edge <- matrix(c(n + 1L, 1L, n + 1L, 2L), ncol = 2, byrow = TRUE)
  next_node <- n + 2L
  state <- index %% 2147483647
  for (k in seq_len(n - 2L) + 2L) {
    state <- .lcg_next(state)
    pick <- (state %% nrow(edge)) + 1L
    parent <- edge[pick, 1]
    child <- edge[pick, 2]
    edge[pick, ] <- c(parent, next_node)
    edge <- rbind(edge, c(next_node, child), c(next_node, k))
    next_node <- next_node + 1L
  }
  tree <- list(edge = edge, tip.label = tipLabels, Nnode = n - 1L)
  class(tree) <- "phylo"
  TreeTools::Preorder(tree)
}

#' Convert phyDat object to the list format expected by ts_* C++ bridges
make_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels
  )
}

#' Score a tree using the C++ Fitch engine
ts_score <- function(tree, ds, concavity = Inf, min_steps = integer(0),
                     infoAmounts = NULL) {
  TreeSearch:::ts_fitch_score(
    tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
    min_steps = min_steps, concavity = concavity,
    infoAmounts = infoAmounts
  )
}

#' Compare two phylogenetic trees topologically (handles order/attr differences).
#' In testthat edition 3, expect_equal() uses waldo which checks all attributes
#' including `order` (cladewise vs preorder). all.equal.phylo handles these.
expect_equal_tree <- function(actual, expected) {
  cmp <- all.equal(actual, expected)
  if (!isTRUE(cmp)) testthat::fail(paste(cmp, collapse = "\n")) else testthat::succeed()
}

#' Validate that a search result has correct tree topology
validate_result <- function(result, n_tip) {
  if ("trees" %in% names(result)) {
    edges <- result$trees[[1]]
  } else {
    edges <- result$edge
  }
  testthat::expect_equal(nrow(edges), 2L * (n_tip - 1L))
  children <- edges[, 2]
  tips <- sort(children[children <= n_tip])
  testthat::expect_equal(tips, seq_len(n_tip))
}
