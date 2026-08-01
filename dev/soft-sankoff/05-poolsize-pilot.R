#!/usr/bin/env Rscript
# PILOT — is `cidRankInMpt` measuring the CRITERION, or the SEARCH's retention?
#
# `cidRankInMpt` (the winner's quantile rank among its own MPT set's
# distance-to-truth distribution) is the statistic carrying the Step 3a headline:
# on congreveLamsdellMatrices it came out at 0.33 against a null of 0.5,
# p ~ 1e-6.  There, MPT sets were mostly complete (median 11 trees).
#
# On the O'Reilly matrices MaximizeParsimony() returns exactly 100 trees, which
# is `SearchControl()$poolMaxSize` -- a HARD CAP, not an exhausted set.  So the
# 100 trees are whichever 100 the search happened to retain, and a rank computed
# within them is a statement about the criterion AND about retention.  If
# retention correlates with anything CID-related, the rank statistic inherits it.
#
# This pilot asks whether the statistic is stable against pool size.  Same
# matrices, two caps.  If mean rank moves with the cap, the statistic is partly
# measuring retention and must not be reported from a capped set without that
# caveat; if it holds, the O'Reilly sweep is worth running.
#
# Cheap by design: a handful of matrices, one temperature, no Mk oracle.
#
# Usage:
#   Rscript dev/soft-sankoff/05-poolsize-pilot.R [nMatrices] [caps...]

args <- commandArgs(trailingOnly = TRUE)
N_MATRICES <- if (length(args) >= 1) as.integer(args[1]) else 10L
CAPS <- if (length(args) >= 2) as.integer(args[-1]) else c(100L, 300L)

TS_LIB <- Sys.getenv("SOFT_SANKOFF_LIB", ".agent-softsankoff")
OUT_DIR <- Sys.getenv("SOFT_SANKOFF_OUT", "dev/soft-sankoff")
OR_ROOT <- Sys.getenv("OR_ROOT", "C:/Users/pjjg18/GitHub/OReillyEtAl2016")

if (nzchar(TS_LIB)) .libPaths(c(TS_LIB, .libPaths()))
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = if (nzchar(TS_LIB)) TS_LIB else NULL)
  library("TreeTools")
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
source(file.path(OR_ROOT, "data", "orReferenceTree.R"))

# A grid, not one T.  The first pass used T = 0.1 alone and answered only the
# low-temperature question; T = 0.5 is where the congreveLamsdellMatrices result
# was strongest, so the generalisation test has to cover it.  Pool construction
# dominates the cost and is shared across temperatures, so the grid is nearly
# free.
TEMPERATURES <- c(0.02, 0.1, 0.25, 0.5, 1)
SEED <- 20260801L
MATRIX_DIR <- file.path(OR_ROOT, "data-raw", "Matrices", "100_char_matrices")

# Per-pattern tip costs, in `tipLabels` order, plus the shared cost matrices.
PatternCosts <- function(dataset, tipLabels) {
  at <- attributes(dataset)
  nStates <- length(at[["levels"]])
  tipData <- matrix(unlist(dataset, use.names = FALSE),
                    nrow = length(dataset), byrow = TRUE,
                    dimnames = list(names(dataset), NULL))
  tipData <- tipData[tipLabels, , drop = FALSE]
  contrast <- at[["contrast"]]
  costs <- lapply(seq_len(ncol(tipData)), function(p) {
    out <- matrix(Inf, length(tipLabels), nStates)
    for (i in seq_along(tipLabels)) {
      out[i, which(contrast[tipData[i, p], ] > 0)] <- 0
    }
    out
  })
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0
  list(costs = costs, weight = as.numeric(at[["weight"]]),
       costMatrices = rep(list(cost), length(costs)))
}

SoftScore <- function(edge, patterns, temperature) {
  r <- TreeSearch:::ts_soft_sankoff_test(
    edge = edge, n_tip = nrow(patterns[["costs"]][[1]]),
    tip_costs = patterns[["costs"]],
    cost_matrices = patterns[["costMatrices"]],
    temperature = temperature)
  sum(r[["per_char"]] * patterns[["weight"]])
}

tipLabels <- orReferenceTree[["tip.label"]]
files <- sort(list.files(MATRIX_DIR, pattern = "\\.NEX$", full.names = TRUE))
files <- files[seq_len(min(N_MATRICES, length(files)))]

rows <- list()
for (fi in seq_along(files)) {
  dataset <- ReadAsPhyDat(files[fi])
  patterns <- PatternCosts(dataset, tipLabels)

  for (cap in CAPS) {
    set.seed(SEED)   # same search seed at both caps, so only the cap differs
    elapsed <- system.time(
      mpt <- MaximizeParsimony(dataset, verbosity = 0, collapse = FALSE,
                               poolMaxSize = cap)
    )[["elapsed"]]
    if (inherits(mpt, "phylo")) mpt <- structure(list(mpt), class = "multiPhylo")
    mpt <- lapply(mpt, function(tr) {
      Preorder(RenumberTips(if (ape::is.binary(tr)) tr else ape::multi2di(tr),
                            tipLabels))
    })

    scores <- vapply(mpt, function(tr) {
      as.numeric(TreeLength(tr, dataset, concavity = Inf))
    }, numeric(1))
    keep <- which(scores <= min(scores) + 1e-9)     # strict MPTs only
    if (length(keep) < 2) next

    cid <- as.numeric(TreeDist::ClusteringInfoDist(
      structure(mpt[keep], class = "multiPhylo"), orReferenceTree,
      normalize = TRUE))
    randomMpt <- cid[sample.int(length(keep), 1)]
    edges <- lapply(mpt[keep], function(tr) {
      e <- tr[["edge"]]; storage.mode(e) <- "integer"; e
    })

    ranks <- numeric(0)
    for (temperature in TEMPERATURES) {
      soft <- vapply(edges, function(e) SoftScore(e, patterns, temperature),
                     numeric(1))
      chosen <- which.min(soft)
      # Rank of the winner among the OTHER MPTs, so a set of size n gives a
      # rank in [0, 1] that is uniform under the null.
      rank <- mean(cid[-chosen] < cid[chosen])
      ranks <- c(ranks, rank)
      rows[[length(rows) + 1]] <- data.frame(
        matrix = basename(files[fi]), cap = cap, temperature = temperature,
        seconds = elapsed, nReturned = length(mpt), nMpt = length(keep),
        optimum = min(scores), cidRankInMpt = rank,
        cidChosen = cid[chosen], cidMptMean = mean(cid),
        cidRandomMpt = randomMpt)
    }
    cat(sprintf("%s cap=%-4d %5.1fs  MPTs %3d  ranks %s\n",
                basename(files[fi]), cap, elapsed, length(keep),
                paste(sprintf("%.2f", ranks), collapse = " ")))
    utils::flush.console()
  }
}

result <- do.call(rbind, rows)
utils::write.csv(result, file.path(OUT_DIR, "05-poolsize-pilot.csv"),
                 row.names = FALSE)

# QUESTION 1: does the rank statistic replicate at 75 tips at all?
cat("\n=== Does cidRankInMpt replicate at 75 tips? ===\n")
cat("congreveLamsdellMatrices (22 tips, 54 patterns): mean rank 0.328,\n")
cat("Wilcoxon vs 0.5 p ~ 1e-6, n = 92.  Below, per temperature, at each cap:\n\n")
byT <- do.call(rbind, lapply(CAPS, function(cp) {
  do.call(rbind, lapply(TEMPERATURES, function(tt) {
    s <- result[result[["cap"]] == cp & result[["temperature"]] == tt, ]
    rk <- s[["cidRankInMpt"]]
    data.frame(cap = cp, temperature = tt, n = length(rk),
               meanRank = mean(rk), medianRank = stats::median(rk),
               rankP = if (length(rk) > 2) {
                 stats::wilcox.test(rk - 0.5)[["p.value"]]
               } else NA_real_,
               chosenBeatsMean = sum(s[["cidChosen"]] < s[["cidMptMean"]]),
               randomBeatsMean = sum(s[["cidRandomMpt"]] < s[["cidMptMean"]]))
  }))
}))
print(byT, row.names = FALSE)
cat("\nRank distribution at the lowest cap (0-.2 .2-.4 .4-.6 .6-.8 .8-1):\n")
for (tt in TEMPERATURES) {
  rk <- result[result[["cap"]] == CAPS[1] & result[["temperature"]] == tt,
               "cidRankInMpt"]
  cat(sprintf("  T=%-5g %s\n", tt,
              paste(table(cut(rk, seq(0, 1, 0.2))), collapse = " ")))
}
cat("A U-shaped distribution means the criterion is choosing DECISIVELY but not\n")
cat("in a truth-correlated way -- near-best about as often as near-worst.  That\n")
cat("is a different failure from choosing at random, and it is what a\n")
cat("density-tracking criterion would look like where density and truth have\n")
cat("come apart.\n")

# QUESTION 2: is the statistic an artifact of the truncated pool?
cat("\n=== Is cidRankInMpt stable against the pool cap? ===\n")
byCap <- do.call(rbind, lapply(CAPS, function(cp) {
  s <- result[result[["cap"]] == cp, ]
  data.frame(cap = cp, n = nrow(s),
             medianReturned = stats::median(s[["nReturned"]]),
             medianMpt = stats::median(s[["nMpt"]]),
             meanRank = mean(s[["cidRankInMpt"]]),
             medianSeconds = stats::median(s[["seconds"]]))
}))
print(byCap, row.names = FALSE)

# Paired across caps on the same matrices, at one temperature so the pairing is
# one-to-one.
PAIR_T <- TEMPERATURES[which.min(abs(TEMPERATURES - 0.1))]
wide <- merge(result[result[["cap"]] == CAPS[1] &
                       result[["temperature"]] == PAIR_T,
                     c("matrix", "cidRankInMpt")],
              result[result[["cap"]] == CAPS[length(CAPS)] &
                       result[["temperature"]] == PAIR_T,
                     c("matrix", "cidRankInMpt")],
              by = "matrix", suffixes = c(".lo", ".hi"))
if (nrow(wide) > 2) {
  cat(sprintf("\nPaired on %d matrices: mean rank %.3f at cap %d vs %.3f at cap %d\n",
              nrow(wide), mean(wide[["cidRankInMpt.lo"]]), CAPS[1],
              mean(wide[["cidRankInMpt.hi"]]), CAPS[length(CAPS)]))
  cat(sprintf("Paired Wilcoxon on the difference: p = %.4g\n",
              stats::wilcox.test(wide[["cidRankInMpt.lo"]],
                                 wide[["cidRankInMpt.hi"]],
                                 paired = TRUE)[["p.value"]]))
  cat("A large p (no shift) means the rank statistic is not an artifact of where\n")
  cat("the search stopped retaining, and the O'Reilly sweep is worth running.\n")
  cat("A shift means it partly measures retention and must not be reported from\n")
  cat("a capped pool without saying so.\n")
}
cat(sprintf("\nAlso note whether nReturned == cap exactly: that is the hard-cap\n"))
cat("signature, and it tells you the set is truncated rather than exhausted.\n")
