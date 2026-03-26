# VTune driver: TBR inner loop profiling at 88 tips (Dikow2009)
# Target: ~15-30s of pure TBR evaluation time
library(TreeSearch, lib.loc = ".vtune-lib")
library(TreeTools)

data(inapplicable.phyData)
dataset <- inapplicable.phyData$Dikow2009

at <- attributes(dataset)
contrast <- at$contrast
tip_data <- matrix(unlist(dataset, use.names = FALSE),
                   nrow = length(dataset), byrow = TRUE)
weight <- at$weight
levels <- at$levels

set.seed(5813)
n_tip <- length(dataset)
t0 <- proc.time()

# Many random starts: Wagner → NNI → TBR chain
# Each start exercises the full TBR evaluation pipeline
for (rep in seq_len(50)) {
  wag <- TreeSearch:::ts_random_wagner_tree(contrast, tip_data, weight, levels)
  nni <- TreeSearch:::ts_nni_search(wag$edge, contrast, tip_data, weight, levels)
  edge <- nni$edge
  for (pass in seq_len(20)) {
    res <- TreeSearch:::ts_tbr_search(
      edge, contrast, tip_data, weight, levels,
      maxHits = 1L, acceptEqual = FALSE
    )
    edge <- res$edge
  }
}

elapsed <- (proc.time() - t0)["elapsed"]
cat("Elapsed:", round(elapsed, 1), "s (", 50 * 20, "TBR passes)\n")
