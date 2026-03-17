# Extracted from test-ts-spr-nni-opt.R:303

# prequel ----------------------------------------------------------------------
ts_spr <- function(tree, ds, maxHits = 20L, concavity = Inf) {
  TreeSearch:::ts_spr_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, concavity = concavity)
}
ts_nni <- function(tree, ds, maxHits = 20L, concavity = Inf) {
  TreeSearch:::ts_nni_search(tree$edge, ds$contrast, ds$tip_data,
                             ds$weight, ds$levels,
                             maxHits = maxHits, concavity = concavity)
}

# test -------------------------------------------------------------------------
set.seed(8115)
mat <- matrix(sample(0:1, 10 * 8, replace = TRUE),
                nrow = 10, dimnames = list(paste0("t", 1:10), NULL))
dataset <- MatrixToPhyDat(mat)
ds <- make_ts_data(dataset)
scores_nni <- numeric(3)
scores_spr <- numeric(3)
for (i in 1:3) {
    set.seed(1000 + i)
    tree <- as.phylo(i * 100, 10)
    scores_nni[i] <- ts_nni(tree, ds, maxHits = 10L)$score
    set.seed(1000 + i)
    scores_spr[i] <- ts_spr(tree, ds, maxHits = 10L)$score
  }
for (i in 1:3) {
    expect_true(scores_spr[i] <= scores_nni[i])
  }
