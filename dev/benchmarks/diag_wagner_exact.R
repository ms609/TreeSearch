# Decisive test: does an EXACT-insertion RAS Wagner (try every edge, full
# TreeLength, pick the true argmin) reach TNT-like Wagner quality (~1300)?
# If yes, our fast insertion-cost formula (fitch_indirect_length:
# Y = final(A) | final(D)) is the bug behind the +30% Wagner deficit.
suppressMessages({
  library(TreeSearch, lib.loc = "C:/Users/pjjg18/GitHub/TS-selectem/.agent-selectem")
  library(TreeTools)
})
nm   <- Sys.getenv("DS", "Zanol2014")
phy  <- readRDS(sprintf("dev/benchmarks/t0/%s.phy.rds", nm))
taxa <- names(phy); n <- length(taxa)
NSEED <- as.integer(Sys.getenv("NSEED", "3"))

ExactWagner <- function(seed) {
  set.seed(seed)
  ord <- sample(taxa)
  tr  <- PectinateTree(ord[1:3])                     # unique 3-tip topology
  for (k in 4:n) {
    tip   <- ord[k]
    phySub <- phy[ord[1:k]]                            # data on tips so far
    nNodeNow <- 2L * (k - 1L) - 1L                     # nodes in current (k-1)-tip tree
    best <- NULL; bestLen <- Inf
    for (w in 0:nNodeNow) {                            # 0 = above root
      cand <- tryCatch(AddTip(tr, where = w, label = tip),
                       error = function(e) NULL)
      if (is.null(cand)) next
      L <- TreeLength(cand, phySub)
      if (L < bestLen) { bestLen <- L; best <- cand }
    }
    tr <- best
  }
  TreeLength(tr, phy)
}

# Reference: our fast-formula kernel Wagner, same seeds
at <- attributes(phy); contrast <- at$contrast
tipData <- matrix(unlist(phy, use.names = FALSE), nrow = n, byrow = TRUE)
weight  <- TreeSearch:::.ScaleWeight(at$weight); levels <- at$levels
FastWagner <- function(seed) {
  set.seed(seed)
  TreeSearch:::ts_random_wagner_tree(contrast, tipData, weight, levels)$score
}

cat(sprintf("== %s | n=%d | exact-insertion vs fast-formula RAS Wagner ==\n", nm, n))
cat("(TNT no-swap RAS Wagner ~1283-1325; optimum ~1261)\n\n")
for (s in seq_len(NSEED)) {
  t0 <- proc.time()["elapsed"]
  ex <- ExactWagner(s)
  el <- proc.time()["elapsed"] - t0
  fa <- FastWagner(s)
  cat(sprintf("  seed %d: exact=%.0f  fast=%.0f  gap(fast-exact)=%+.0f  [%.0fs]\n",
              s, ex, fa, fa - ex, el))
}
