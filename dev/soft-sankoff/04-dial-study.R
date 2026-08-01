#!/usr/bin/env Rscript
# STEP 3a — the parsimony/likelihood dial study.
#
# The standing methodological weakness of the parsimony-vs-likelihood
# literature is that every comparison confounds CRITERION with IMPLEMENTATION:
# different programs, different search intensity, different rooting
# conventions, different stopping rules.  A temperature parameter dissolves
# that confound.  Same dynamic program, same data, same candidate trees; one
# scalar varied.
#
# Question: what T maximises recovery of the generating tree, and does the
# answer track homoplasy?
#
# WHY THIS IS NOT GATE A.  Gate A asked whether a soft-guided SEARCH would walk
# toward truth or toward MPT density.  That question existed to protect Step 4
# (annealed search), which 03-speed-budget.R killed on cost: a compiled soft
# score costs a median x1147 of a Fitch score for binary characters, against a
# x50 orientation threshold.  Gate A is therefore retired rather than answered,
# and 02-tilt-direction.R is left as the record of what it said.  This script
# asks the criterion question instead, which survives Gate B because it never
# needed the score to be cheap.
#
# WHAT THIS DESIGN CAN AND CANNOT CLAIM.  Search is held constant on purpose:
# the candidate pool comes from ordinary hard-parsimony search, and T only
# chooses among its members.  So this measures which tree the criterion at T
# PREFERS out of a common pool.  It does not measure what a soft-objective
# search would find, and the pool inherits hard parsimony's bias about which
# trees are worth having in it.  That is the price of dissolving the
# implementation confound, and it is the right price here, but it is not
# nothing: a criterion can only be credited with recovering a tree the pool
# contains.
#
# Distance to the generating tree is ClusteringInfoDist, normalised.  Not
# Robinson-Foulds, which inflates and is dominated by rogue tips.
#
# Usage:
#   Rscript dev/soft-sankoff/04-dial-study.R [nMatrices] [poolSize] [mkOracle]
#     nMatrices  how many of the 100 matrices to use (default all)
#     poolSize   target candidate trees per matrix (default 200)
#     mkOracle   "yes" to compute the Mk likelihood of each selected tree
#                (default yes; needs phangorn, and it is the slow part)

args <- commandArgs(trailingOnly = TRUE)
N_MATRICES <- if (length(args) >= 1) as.integer(args[1]) else 100L
POOL_SIZE <- if (length(args) >= 2) as.integer(args[2]) else 200L
MK_ORACLE <- if (length(args) >= 3) !identical(args[3], "no") else TRUE

# Library and output locations are overridable so the same script runs against
# the local isolated build and against a Hamilton project library, with no
# cluster-specific fork to drift out of step.
TS_LIB <- Sys.getenv("SOFT_SANKOFF_LIB", ".agent-softsankoff")
OUT_DIR <- Sys.getenv("SOFT_SANKOFF_OUT", "dev/soft-sankoff")

if (nzchar(TS_LIB)) .libPaths(c(TS_LIB, .libPaths()))
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = if (nzchar(TS_LIB)) TS_LIB else NULL)
  library("TreeTools")
})
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TEMPERATURES <- c(0, 0.02, 0.05, 0.1, 0.25, 0.5, 1, 2)
MAX_EXTRA_STEPS <- 5      # pool admits trees within this many steps of optimal
PERTURB_ATTEMPTS <- 40L   # rearrangement attempts per pool slot before giving up
SEED <- 20260801L

data("congreveLamsdellMatrices", package = "TreeSearch")
data("referenceTree", package = "TreeSearch")

# ---------------------------------------------------------------------------
# Soft scoring of a whole matrix
# ---------------------------------------------------------------------------

# Per-pattern tip-cost matrices for a dataset, in the tip order `tipLabels`.
# Returns the cost matrices plus the pattern weights, since phyDat compresses
# to unique patterns and the score is the weighted sum over them.
PatternCosts <- function(dataset, tipLabels) {
  at <- attributes(dataset)
  nStates <- length(at[["levels"]])
  # tip_data: one row per tip, one column per pattern, entries indexing levels.
  tipData <- matrix(unlist(dataset, use.names = FALSE),
                    nrow = length(dataset), byrow = TRUE,
                    dimnames = list(names(dataset), NULL))
  tipData <- tipData[tipLabels, , drop = FALSE]
  contrast <- at[["contrast"]]

  costs <- lapply(seq_len(ncol(tipData)), function(p) {
    out <- matrix(Inf, length(tipLabels), nStates)
    for (i in seq_along(tipLabels)) {
      admissible <- which(contrast[tipData[i, p], ] > 0)
      out[i, admissible] <- 0
    }
    out
  })
  list(costs = costs, weight = as.numeric(at[["weight"]]), nStates = nStates)
}

# Weighted soft score of one tree.  One binding call per tree per temperature;
# per_char is dotted with the pattern weights.
SoftScore <- function(edge, patterns, temperature) {
  nStates <- patterns[["nStates"]]
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0
  result <- TreeSearch:::ts_soft_sankoff_test(
    edge = edge,
    n_tip = nrow(patterns[["costs"]][[1]]),
    tip_costs = patterns[["costs"]],
    cost_matrices = rep(list(cost), length(patterns[["costs"]])),
    temperature = temperature
  )
  sum(result[["per_char"]] * patterns[["weight"]])
}

# ---------------------------------------------------------------------------
# Candidate pool
# ---------------------------------------------------------------------------

# Optimal trees plus near-optimal neighbours, all renumbered to one tip order
# so a single set of tip-cost matrices serves every member.
BuildPool <- function(dataset, tipLabels, poolSize, maxExtra) {
  # collapse = FALSE is REQUIRED, not a preference.  The default contracts
  # zero-length (unsupported) branches before returning, so MaximizeParsimony
  # hands back polytomies on any matrix with an unsupported node -- and both
  # TreeLength() and the soft kernel refuse a non-binary tree.  With the default
  # this script dies partway through the sweep ("`tree` must be binary"), which
  # on 100 matrices means discovering it a quarter of an hour in.
  optimal <- MaximizeParsimony(dataset, verbosity = 0, collapse = FALSE)
  if (inherits(optimal, "phylo")) optimal <- structure(list(optimal),
                                                       class = "multiPhylo")

  # Belt and braces: resolve anything still polytomous rather than trusting the
  # argument to keep meaning what it means.  An arbitrary resolution of an
  # unsupported node is exactly what collapse = FALSE returns anyway, so this
  # changes nothing when the argument works -- it just refuses to fail late.
  Binary <- function(tr) if (ape::is.binary(tr)) tr else ape::multi2di(tr)
  optimal <- lapply(optimal, Binary)

  optimum <- min(vapply(optimal, function(tr) {
    as.numeric(TreeLength(tr, dataset, concavity = Inf))
  }, numeric(1)))

  Canonical <- function(tr) Preorder(RenumberTips(Binary(tr), tipLabels))
  pool <- lapply(optimal, Canonical)
  scores <- rep(optimum, length(pool))
  seen <- vapply(pool, function(tr) paste(as.character(TreeTools::as.Splits(tr)),
                                          collapse = "|"), character(1))

  attempts <- 0L
  maxAttempts <- poolSize * PERTURB_ATTEMPTS
  while (length(pool) < poolSize && attempts < maxAttempts) {
    attempts <- attempts + 1L
    parent <- pool[[sample.int(length(pool), 1)]]
    candidate <- try(Canonical(TBR(parent)), silent = TRUE)
    if (inherits(candidate, "try-error")) next
    score <- as.numeric(TreeLength(candidate, dataset, concavity = Inf))
    if (score > optimum + maxExtra) next
    key <- paste(as.character(TreeTools::as.Splits(candidate)), collapse = "|")
    if (key %in% seen) next
    seen <- c(seen, key)
    pool[[length(pool) + 1]] <- candidate
    scores <- c(scores, score)
  }

  list(trees = pool, scores = scores, optimum = optimum,
       nOptimal = length(optimal))
}

# ---------------------------------------------------------------------------
# Mk likelihood oracle (secondary, descriptive)
# ---------------------------------------------------------------------------

# Parsimony trees carry no branch lengths, and phangorn::pml() refuses a tree
# without them ("tree must have edge weights").  Seed them, then let optim.pml
# fit them.
#
# The first failure is reported rather than swallowed.  An earlier version wrapped
# this in a bare try() and the whole mkLogLik column came back NA on a run that
# otherwise looked successful; an oracle that silently degrades to NA is worse
# than one that is switched off.
mkWarned <- FALSE
MkLogLik <- function(tree, dataset) {
  if (!requireNamespace("phangorn", quietly = TRUE)) {
    if (!mkWarned) {
      warning("phangorn unavailable; Mk oracle skipped", call. = FALSE)
      mkWarned <<- TRUE
    }
    return(NA_real_)
  }
  tree[["edge.length"]] <- rep(0.1, nrow(tree[["edge"]]))
  fit <- tryCatch(
    suppressWarnings(phangorn::optim.pml(
      phangorn::pml(tree, dataset), optEdge = TRUE,
      control = phangorn::pml.control(trace = 0))),
    error = function(e) e)
  if (inherits(fit, "error")) {
    if (!mkWarned) {
      warning("Mk oracle failed: ", conditionMessage(fit), call. = FALSE)
      mkWarned <<- TRUE
    }
    return(NA_real_)
  }
  as.numeric(fit[["logLik"]])
}

# ---------------------------------------------------------------------------
# Homoplasy
# ---------------------------------------------------------------------------

# Consistency index: minimum possible steps over observed steps at the optimum.
# 1 = no homoplasy; lower = more.  Computed from the contrast matrix, so it does
# not depend on the search.
ConsistencyIndex <- function(dataset, optimum) {
  at <- attributes(dataset)
  contrast <- at[["contrast"]]
  minSteps <- pmax(rowSums(contrast > 0) - 1L, 0L)
  tipData <- matrix(unlist(dataset, use.names = FALSE),
                    nrow = length(dataset), byrow = TRUE)
  perPattern <- vapply(seq_len(ncol(tipData)), function(p) {
    max(minSteps[tipData[, p]], length(unique(tipData[, p])) - 1L)
  }, numeric(1))
  sum(perPattern * as.numeric(at[["weight"]])) / optimum
}

# ---------------------------------------------------------------------------
# Main sweep
# ---------------------------------------------------------------------------

set.seed(SEED)
tipLabels <- referenceTree[["tip.label"]]
nUse <- min(N_MATRICES, length(congreveLamsdellMatrices))
rows <- list()

for (m in seq_len(nUse)) {
  dataset <- congreveLamsdellMatrices[[m]]
  patterns <- PatternCosts(dataset, tipLabels)
  pool <- BuildPool(dataset, tipLabels, POOL_SIZE, MAX_EXTRA_STEPS)
  ci <- ConsistencyIndex(dataset, pool[["optimum"]])

  edges <- lapply(pool[["trees"]], function(tr) {
    e <- tr[["edge"]]
    storage.mode(e) <- "integer"
    e
  })
  # Distance of every pool member to the generating tree, computed once.
  cid <- as.numeric(TreeDist::ClusteringInfoDist(
    structure(pool[["trees"]], class = "multiPhylo"), referenceTree,
    normalize = TRUE))

  for (temperature in TEMPERATURES) {
    soft <- vapply(edges, function(e) SoftScore(e, patterns, temperature),
                   numeric(1))
    # Ties are real and common at T = 0 (the MPT set); average the distance
    # over every tree the criterion cannot separate, rather than letting pool
    # order decide the winner.
    best <- which(soft <= min(soft) + 1e-9)
    chosen <- best[1]

    mkLogLik <- NA_real_
    if (MK_ORACLE) {
      mkLogLik <- MkLogLik(pool[["trees"]][[chosen]], dataset)
    }

    rows[[length(rows) + 1]] <- data.frame(
      matrix = m,
      temperature = temperature,
      optimum = pool[["optimum"]],
      consistencyIndex = ci,
      poolSize = length(pool[["trees"]]),
      nOptimal = pool[["nOptimal"]],
      nTied = length(best),
      cidChosen = cid[chosen],
      cidTiedMean = mean(cid[best]),
      cidPoolBest = min(cid),
      cidPoolMean = mean(cid),
      hardScoreChosen = pool[["scores"]][chosen],
      mkLogLik = mkLogLik
    )
  }
  cat(sprintf("matrix %3d/%3d  optimum %5.0f  CI %.3f  pool %3d (%2d optimal)\n",
              m, nUse, pool[["optimum"]], ci, length(pool[["trees"]]),
              pool[["nOptimal"]]))
  utils::flush.console()
}

result <- do.call(rbind, rows)
utils::write.csv(result, file.path(OUT_DIR, "04-dial-study.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

cat("\n=== Recovery of the generating tree by temperature ===\n")
cat("(cidTiedMean: distance to truth, averaged over trees the criterion at T\n")
cat(" cannot separate.  Lower is better.  T = 0 is hard parsimony.)\n\n")
byT <- do.call(rbind, lapply(TEMPERATURES, function(tt) {
  sub <- result[result[["temperature"]] == tt, ]
  data.frame(temperature = tt,
             medianCid = stats::median(sub[["cidTiedMean"]]),
             meanCid = mean(sub[["cidTiedMean"]]),
             medianTied = stats::median(sub[["nTied"]]),
             medianMkLogLik = stats::median(sub[["mkLogLik"]], na.rm = TRUE))
}))
print(byT, row.names = FALSE)

nMk <- sum(!is.na(result[["mkLogLik"]]))
if (MK_ORACLE) {
  cat(sprintf("\nMk oracle: %d of %d selections fitted (%s)\n",
              nMk, nrow(result),
              if (nMk == 0) "COLUMN IS EMPTY -- treat medianMkLogLik as absent"
              else "secondary descriptive only"))
}

baseline <- byT[byT[["temperature"]] == 0, "medianCid"]
cat(sprintf("\nHard parsimony (T = 0) median distance to truth: %.4f\n", baseline))
better <- byT[byT[["medianCid"]] < baseline & byT[["temperature"]] > 0, ]
if (nrow(better)) {
  cat("Temperatures beating it on the median:\n")
  print(better, row.names = FALSE)
} else {
  cat("No T > 0 beats hard parsimony on the median distance to truth.\n")
}

# Per-matrix sign test: does the best T > 0 beat T = 0 more often than chance?
perMatrix <- do.call(rbind, lapply(unique(result[["matrix"]]), function(m) {
  sub <- result[result[["matrix"]] == m, ]
  hard <- sub[sub[["temperature"]] == 0, "cidTiedMean"]
  warm <- sub[sub[["temperature"]] > 0, ]
  bestRow <- warm[which.min(warm[["cidTiedMean"]]), ]
  data.frame(matrix = m, consistencyIndex = sub[["consistencyIndex"]][1],
             hardCid = hard, bestT = bestRow[["temperature"]],
             bestCid = bestRow[["cidTiedMean"]])
}))
utils::write.csv(perMatrix, file.path(OUT_DIR, "04-dial-study-per-matrix.csv"),
                 row.names = FALSE)

wins <- sum(perMatrix[["bestCid"]] < perMatrix[["hardCid"]] - 1e-9)
ties <- sum(abs(perMatrix[["bestCid"]] - perMatrix[["hardCid"]]) <= 1e-9)
losses <- nrow(perMatrix) - wins - ties
cat(sprintf("\nPer-matrix, best T > 0 vs hard parsimony: %d better, %d tied, %d worse\n",
            wins, ties, losses))
cat("NOTE: bestT is selected per matrix using the answer, so this is an\n")
cat("optimistic ceiling on what a fixed T could achieve, not an estimate of it.\n")
if (wins + losses > 0) {
  cat(sprintf("Sign test p = %.4g\n",
              stats::binom.test(wins, wins + losses)[["p.value"]]))
}

cat("\n=== Does the best temperature track homoplasy? ===\n")
warmOnly <- perMatrix[perMatrix[["bestCid"]] < perMatrix[["hardCid"]] - 1e-9, ]
if (nrow(warmOnly) >= 5) {
  rho <- stats::cor(warmOnly[["consistencyIndex"]], warmOnly[["bestT"]],
                    method = "spearman")
  cat(sprintf("Spearman(consistency index, best T) over %d matrices where a\n",
              nrow(warmOnly)))
  cat(sprintf("warm T helped: rho = %+.3f\n", rho))
  cat("Consistency index falls as homoplasy rises, so a NEGATIVE rho means\n")
  cat("more homoplasy favours a HIGHER temperature -- the plan's prediction.\n")
} else {
  cat(sprintf("Only %d matrices had any T > 0 beat hard parsimony; too few to\n",
              nrow(warmOnly)))
  cat("correlate best T against homoplasy.\n")
}
