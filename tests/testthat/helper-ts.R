# Shared helpers for ts-* test files.
# testthat auto-sources helper-*.R files before running tests.

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
