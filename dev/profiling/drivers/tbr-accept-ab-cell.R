# One A/B cell for the TBR-rerooting accept path — issue #38 (T-300)
#
# Prints ONE line so a shell loop can interleave the arms (this machine is
# shared, so alternating arms is what stops a drifting background load being
# read as an effect).  The reported numbers come from the DLL's own
# steady_clock counters around the accept branch, not from R-level timing.
#
# Usage: TS_NA_TIMING=1 Rscript dev/profiling/drivers/tbr-accept-ab-cell.R \
#          <lib> <arm> <dataset> <mode> <seed> [nCycles]

args <- commandArgs(trailingOnly = TRUE)
libDir <- args[[1]]
arm <- args[[2]]
dsName <- args[[3]]
mode <- args[[4]]
seed <- as.integer(args[[5]])
nCycles <- if (length(args) >= 6L) as.integer(args[[6]]) else 8L

suppressMessages(library(TreeSearch, lib.loc = libDir))
stopifnot(nzchar(Sys.getenv("TS_NA_TIMING")))

dataset <- TreeSearch::inapplicable.phyData[[dsName]]
at <- attributes(dataset)
tipData <- matrix(unlist(dataset, use.names = FALSE),
                  nrow = length(dataset), byrow = TRUE)
weight <- TreeSearch:::.ScaleWeight(at$weight)
concavity <- if (mode == "EW") -1 else 10

set.seed(seed)
tr <- ape::rtree(length(dataset), tip.label = names(dataset), rooted = FALSE)
startEdge <- ape::root(tr, 1L, resolve.root = TRUE)$edge

set.seed(seed)
res <- TreeSearch:::ts_ratchet_search(
  edge = startEdge, contrast = at$contrast, tip_data = tipData,
  weight = weight, levels = at$levels,
  nCycles = nCycles, perturbProb = 0.04, maxHits = 1L, concavity = concavity)

cat(sprintf("%s,%s,%s,%d,%.4f,%.2f,%.2f,%.0f\n",
            arm, dsName, mode, seed, res$score,
            res$na_t_total_ms, res$na_t_accept_ms, res$na_n_accept))
