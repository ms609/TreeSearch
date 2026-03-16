# Extracted from test-ts-simplify.R:275

# prequel ----------------------------------------------------------------------
library("TreeTools")
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
ts_score <- function(tree, ds, concavity = Inf, min_steps = integer(0)) {
  TreeSearch:::ts_fitch_score(tree$edge, ds$contrast, ds$tip_data,
                               ds$weight, ds$levels,
                               min_steps = min_steps,
                               concavity = concavity)
}
ts_diag <- function(ds) {
  TreeSearch:::ts_simplify_diag(ds$contrast, ds$tip_data,
                                 ds$weight, ds$levels)
}
ts_driven <- function(ds, ...) {
  TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 3L, targetHits = 1L,
    ratchetCycles = 1L, xssRounds = 0L,
    xssPartitions = 2L, fuseInterval = 10L,
    maxSeconds = 0, verbosity = 0L,
    ...
  )
}
info_mat <- matrix(c(
  0, 0, 0, 0, 0, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 1, 0, 0, 1, 1, 1,
  0, 1, 0, 1, 1, 0, 1, 0, 1, 1,
  1, 0, 0, 0, 1, 1, 1, 1, 0, 0
), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
info_dataset <- MatrixToPhyDat(info_mat)
info_ds <- make_ts_data(info_dataset)
autap_mat <- matrix(c(
  # Char 1: informative (0/1 split)
  0, 0, 0, 0, 1, 1, 1, 1,
  # Char 2: single autapomorphy (tip 1 = 2, rest = 0) -> uninformative
  2, 0, 0, 0, 0, 0, 0, 0,
  # Char 3: two autapomorphies (tips 1,2 unique) -> uninformative
  1, 2, 0, 0, 0, 0, 0, 0,
  # Char 4: informative (0/1 split)
  0, 0, 1, 1, 0, 0, 1, 1
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
autap_dataset <- MatrixToPhyDat(autap_mat)
autap_ds <- make_ts_data(autap_dataset)
singleton_mat <- matrix(c(
  0, 0, 0, 1, 1, 1, 2,
  0, 0, 1, 1, 0, 0, 1
), nrow = 7, dimnames = list(paste0("t", 1:7), NULL))
singleton_dataset <- MatrixToPhyDat(singleton_mat)
singleton_ds <- make_ts_data(singleton_dataset)
invar_mat <- matrix(c(
  0, 0, 0, 0, 0, 0, 0, 0,  # invariant
  0, 0, 0, 0, 1, 1, 1, 1   # informative
), nrow = 8, dimnames = list(paste0("t", 1:8), NULL))
invar_dataset <- MatrixToPhyDat(invar_mat)
invar_ds <- make_ts_data(invar_dataset)

# test -------------------------------------------------------------------------
set.seed(9201)
mat4 <- matrix(c(
    0, 0, 0, 0, 1, 1, 1, 1, 2, 3,  # states 2,3 are singletons
    0, 0, 1, 1, 0, 0, 1, 1, 0, 0,  # informative binary
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0,  # informative binary
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1   # invariant
  ), nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
ds4 <- MatrixToPhyDat(mat4)
ds4_data <- make_ts_data(ds4)
diag <- ts_diag(ds4_data)
expect_gte(diag$n_patterns_removed, 1L)
result <- ts_driven(ds4_data, maxReplicates = 3L)
