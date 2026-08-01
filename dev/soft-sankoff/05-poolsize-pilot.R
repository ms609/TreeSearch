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

TEMPERATURE <- 0.1        # inside the MPT set on C-L: 0/100 selections escaped
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
    soft <- vapply(mpt[keep], function(tr) {
      e <- tr[["edge"]]; storage.mode(e) <- "integer"
      SoftScore(e, patterns, TEMPERATURE)
    }, numeric(1))

    chosen <- which.min(soft)
    rows[[length(rows) + 1]] <- data.frame(
      matrix = basename(files[fi]), cap = cap, seconds = elapsed,
      nReturned = length(mpt), nMpt = length(keep), optimum = min(scores),
      cidRankInMpt = mean(cid[keep != keep[chosen]] < cid[chosen]),
      cidChosen = cid[chosen], cidMptMean = mean(cid),
      cidRandomMpt = cid[sample.int(length(keep), 1)])
    cat(sprintf("%s cap=%-4d %5.1fs  returned %3d  MPTs %3d  rank %.3f\n",
                basename(files[fi]), cap, elapsed, length(mpt), length(keep),
                rows[[length(rows)]][["cidRankInMpt"]]))
    utils::flush.console()
  }
}

result <- do.call(rbind, rows)
utils::write.csv(result, file.path(OUT_DIR, "05-poolsize-pilot.csv"),
                 row.names = FALSE)

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

# Paired across caps on the same matrices: the honest comparison.
wide <- merge(result[result[["cap"]] == CAPS[1], c("matrix", "cidRankInMpt")],
              result[result[["cap"]] == CAPS[length(CAPS)],
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
