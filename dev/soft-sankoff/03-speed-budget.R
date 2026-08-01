#!/usr/bin/env Rscript
# GATE B — what does a soft score cost against SIMD Fitch?
#
# Two distinct costs, and the second is likely the larger:
#   1. Loss of bit-parallelism: O(k^2) doubles plus exp/log per character per
#      node, against word-packed bitwise Fitch (src/ts_simd.h).
#   2. Loss of incremental rescoring.  src/ts_sankoff.h is full-rescore only;
#      TBR draws much of its throughput from rescoring only the affected
#      subtree.  That changes the complexity of the inner loop, not a constant.
#
# HONESTY NOTE: the wall-clock ratio below compares *pure R* soft-Sankoff with
# *compiled SIMD* Fitch, so it is an upper bound on the real cost and not the
# Gate B number.  The implementation-independent quantity is the operation
# count, reported alongside it.  The real number needs a C++ prototype; this
# script exists to bound the problem and to size that prototype.
#
# Usage:
#   Rscript dev/soft-sankoff/03-speed-budget.R

source("tests/testthat/helper-soft-sankoff.R")
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = ".agent-softsankoff")
})

set.seed(2026)
LADDER <- c(16, 32, 64, 128)
N_CHARS <- 100L
N_STATES <- 2L
REPS <- 5L

MakeData <- function(nTip, nChars, nStates) {
  mat <- matrix(sample(as.character(seq_len(nStates) - 1L),
                       nTip * nChars, replace = TRUE), nrow = nTip)
  rownames(mat) <- paste0("t", seq_len(nTip))
  phangorn::phyDat(mat, type = "USER",
                   levels = as.character(seq_len(nStates) - 1L))
}

cost <- matrix(1, N_STATES, N_STATES)
diag(cost) <- 0

rows <- list()
for (nTip in LADDER) {
  dataset <- MakeData(nTip, N_CHARS, N_STATES)
  tree <- TreeTools::Preorder(ape::rtree(nTip, rooted = TRUE, br = NULL))
  tree[["tip.label"]] <- paste0("t", seq_len(nTip))
  chars <- as.character(dataset)[tree[["tip.label"]], , drop = FALSE]
  coded <- matrix(match(chars, c("0", "1")), nrow = nTip)

  fitchTime <- system.time(
    for (i in seq_len(REPS * 20L)) TreeLength(tree, dataset, concavity = Inf)
  )[["elapsed"]] / (REPS * 20L)

  softTime <- system.time(
    for (i in seq_len(REPS)) {
      for (j in seq_len(ncol(coded))) {
        SoftSankoffScore(tree[["edge"]],
                         TipCosts(as.list(coded[, j]), N_STATES),
                         cost, temperature = 0.25)
      }
    }
  )[["elapsed"]] / REPS

  # Implementation-independent op counts, per whole-tree score.
  nNode <- nTip - 1L
  # Fitch: one bitwise AND/OR pass per node, over ceil(nChars/64) words.
  fitchOps <- nNode * ceiling(N_CHARS / 64)
  # Soft-Sankoff: k^2 adds + k exp + 1 log, per character, per node.
  softOps <- nNode * N_CHARS * (N_STATES^2 + N_STATES + 1)

  rows[[length(rows) + 1]] <- data.frame(
    nTip = nTip,
    fitch_s = fitchTime, soft_s = softTime,
    wall_ratio = softTime / fitchTime,
    fitch_ops = fitchOps, soft_ops = softOps,
    op_ratio = softOps / fitchOps
  )
  cat(sprintf(
    "nTip=%-4d  Fitch %8.2f us   softR %9.2f us   wall x%-8.0f  ops x%.0f\n",
    nTip, fitchTime * 1e6, softTime * 1e6, softTime / fitchTime,
    softOps / fitchOps))
}

result <- do.call(rbind, rows)
utils::write.csv(result, "dev/soft-sankoff/03-speed-budget.csv",
                 row.names = FALSE)

cat("\n=== Gate B ===\n")
cat(sprintf("Operation-count ratio (implementation independent): x%.0f-%.0f\n",
            min(result[["op_ratio"]]), max(result[["op_ratio"]])))
cat(sprintf("Pure-R wall ratio (UPPER BOUND, not the gate): x%.0f-%.0f\n",
            min(result[["wall_ratio"]]), max(result[["wall_ratio"]])))
cat("\nOrientation threshold from the plan: if a compiled soft score costs\n")
cat("more than ~50x a Fitch score, an annealed search cannot pay for itself\n")
cat("against the ratchet and Step 4 is dead regardless of Gate A.\n")
cat("\nStill unmeasured, and probably dominant: the loss of incremental\n")
cat("rescoring under TBR.  That needs the C++ prototype.\n")
