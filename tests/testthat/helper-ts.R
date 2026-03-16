# Shared helpers for ts-* test files.
# testthat auto-sources helper-*.R files before running tests.

library("TreeTools")

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
