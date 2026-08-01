#!/usr/bin/env Rscript
# GATE A — does the soft-score tilt point toward the truth, or toward density?
#
# softmin rewards topologies admitting MORE near-optimal reconstructions, and
# reconstruction ambiguity is exactly what inflates MPT sets.  So smoothing may
# walk the search *toward* the broad plateaux it was meant to escape.  This
# script measures the sign.
#
# Design: collect an MPT set on a Congreve & Lamsdell matrix (generating tree
# known), rank the MPTs by soft score at a grid of T, and correlate that rank
# with distance to the generating tree and with Mk log likelihood.
#
#   PASS: correlation materially NEGATIVE at some T (better soft score =>
#         closer to truth), consistently across matrices.
#   FAIL: null or positive.  Steps 4 and 5 of the plan die; Steps 2 and 3 live.
#
# Costs are equal and symmetric throughout, so the DP is rooting-invariant and
# the T-374/T-385 rooting defects do not apply.
#
# Usage:
#   Rscript dev/soft-sankoff/02-tilt-direction.R [nMatrices] [maxHits]

source("tests/testthat/helper-soft-sankoff.R")
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = ".agent-softsankoff")
  library("TreeDist")
  # phangorn must be attached, not just namespaced: `as.character.phyDat` is
  # registered but not exported, so `as.character(dataset)` silently returns a
  # bare vector unless the package is on the search path.
  library("phangorn")
})

args <- commandArgs(trailingOnly = TRUE)
nMatrices <- if (length(args) >= 1) as.integer(args[1]) else 5L
# v2.0.0 replaced maxHits/ratchIter with maxReplicates/targetHits/maxSeconds.
maxReplicates <- if (length(args) >= 2) as.integer(args[2]) else 48L

data("congreveLamsdellMatrices", package = "TreeSearch")
data("referenceTree", package = "TreeSearch")

TEMPERATURES <- c(0.02, 0.05, 0.1, 0.25, 0.5, 1)

# --- adapter: phyDat -> pattern-weighted tip-state list -------------------
PatternStates <- function(dataset, tipLabels) {
  chars <- as.character(dataset)[tipLabels, , drop = FALSE]
  levels <- sort(unique(as.vector(chars[!is.na(chars)])))
  coded <- matrix(match(chars, levels), nrow = nrow(chars))
  # Collapse duplicate site patterns; ambiguity/missing becomes "any state".
  key <- apply(coded, 2, paste, collapse = ",")
  keep <- !duplicated(key)
  weight <- as.vector(table(key)[key[keep]])
  list(coded = coded[, keep, drop = FALSE], weight = weight,
       nStates = length(levels))
}

# Weighted soft score of one tree over all patterns.
SoftTreeScore <- function(edge, pat, cost, temperature) {
  total <- 0
  for (j in seq_len(ncol(pat[["coded"]]))) {
    states <- lapply(pat[["coded"]][, j], function(s) {
      if (is.na(s)) seq_len(pat[["nStates"]]) else s
    })
    total <- total + pat[["weight"]][j] *
      SoftSankoffScore(edge, TipCosts(states, pat[["nStates"]]), cost,
                       temperature)
  }
  total
}

rows <- list()
for (m in seq_len(nMatrices)) {
  dataset <- congreveLamsdellMatrices[[m]]
  cat(sprintf("\n--- matrix %d ---\n", m))

  mpts <- MaximizeParsimony(dataset, concavity = Inf, verbosity = 0,
                            maxReplicates = maxReplicates, nThreads = 1L)
  mpts <- unique(mpts)
  if (inherits(mpts, "phylo")) mpts <- structure(list(mpts), class = "multiPhylo")
  cat(sprintf("  %d distinct MPTs, score %s\n", length(mpts),
              format(TreeLength(mpts[[1]], dataset, concavity = Inf))))
  if (length(mpts) < 8) {
    cat("  too few MPTs to correlate; skipping\n")
    next
  }

  tipLabels <- mpts[[1]][["tip.label"]]
  pat <- PatternStates(dataset, tipLabels)
  cost <- matrix(1, pat[["nStates"]], pat[["nStates"]])
  diag(cost) <- 0

  rooted <- lapply(mpts, function(tr) {
    TreeTools::Preorder(ape::multi2di(
      ape::root(tr, outgroup = tipLabels[1], resolve.root = TRUE)))
  })

  truth <- list(
    rf = as.numeric(TreeDist::RobinsonFoulds(mpts, referenceTree)),
    cid = as.numeric(TreeDist::ClusteringInfoDistance(mpts, referenceTree))
  )
  mk <- vapply(mpts, function(tr) {
    fit <- phangorn::pml(ape::compute.brlen(tr, 0.1), dataset)
    as.numeric(phangorn::optim.pml(fit, optEdge = TRUE,
                                   control = phangorn::pml.control(trace = 0)
                                   )[["logLik"]])
  }, numeric(1))

  for (temperature in TEMPERATURES) {
    soft <- vapply(rooted, function(tr) {
      SoftTreeScore(tr[["edge"]], pat, cost, temperature)
    }, numeric(1))
    if (stats::sd(soft) < 1e-12) {
      cat(sprintf("  T = %-6g soft score constant across MPTs (no tilt)\n",
                  temperature))
      next
    }
    rows[[length(rows) + 1]] <- data.frame(
      matrix = m, temperature = temperature, nMPT = length(mpts),
      # Negative => a better (lower) soft score means a closer tree: PASS.
      rho_rf = stats::cor(soft, truth[["rf"]], method = "spearman"),
      rho_cid = stats::cor(soft, truth[["cid"]], method = "spearman"),
      # Negative => a better soft score means a higher Mk likelihood: PASS.
      rho_mk = stats::cor(soft, -mk, method = "spearman")
    )
    cat(sprintf("  T = %-6g rho(RF) = %+.3f  rho(CID) = %+.3f  rho(Mk) = %+.3f\n",
                temperature, rows[[length(rows)]][["rho_rf"]],
                rows[[length(rows)]][["rho_cid"]],
                rows[[length(rows)]][["rho_mk"]]))
  }
}

if (!length(rows)) {
  cat("\nNo matrix produced a usable MPT set. Increase maxHits.\n")
  quit(status = 1)
}

result <- do.call(rbind, rows)
cat("\n=== Gate A summary: mean Spearman rho by temperature ===\n")
cat("(negative = soft score points toward the truth = PASS)\n\n")
summary <- aggregate(cbind(rho_rf, rho_cid, rho_mk) ~ temperature,
                     data = result, FUN = mean)
print(summary, row.names = FALSE, digits = 3)

outFile <- "dev/soft-sankoff/02-tilt-direction.csv"
utils::write.csv(result, outFile, row.names = FALSE)
cat(sprintf("\nPer-matrix rows written to %s\n", outFile))

nUsable <- length(unique(result[["matrix"]]))
verdict <- min(summary[["rho_cid"]])
mkAt <- summary[["rho_mk"]][which.min(summary[["rho_cid"]])]

cat(sprintf("\nUsable matrices (>= 8 distinct MPTs): %d of %d\n",
            nUsable, nMatrices))
cat(sprintf("Best (most negative) mean rho(CID): %+.3f   rho(Mk) there: %+.3f\n",
            verdict, mkAt))

# A single matrix is an anecdote.  Do not let the script declare a gate passed
# on one data point, and do not let a topology signal that disagrees with the
# likelihood signal read as agreement.
cat(sprintf("\nVerdict: %s\n",
            if (nUsable < 5) {
              paste0("INSUFFICIENT DATA - only ", nUsable,
                     " matrix/matrices had a usable MPT set. Raise ",
                     "`nMatrices`, or widen the design to near-optimal trees ",
                     "(see the plan's Gate A note).")
            } else if (verdict < -0.2 && mkAt < 0) {
              "GATE A PASS - tilt points toward the truth on both measures"
            } else if (verdict < -0.2) {
              "GATE A SPLIT - topology measure favourable, Mk measure is not"
            } else if (verdict > 0.2) {
              "GATE A FAIL - tilts toward density"
            } else {
              "GATE A INCONCLUSIVE - no material tilt"
            }))
