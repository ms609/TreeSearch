#!/usr/bin/env Rscript
# Step 0 evidence: the semiring identity, checked numerically.
#
#   T -> 0  :  soft-Sankoff == hard weighted parsimony
#   T  = 1  :  soft-Sankoff with cost = -log(P) == Felsenstein pruning
#
# Requires no build and no TreeSearch install; `ape` is used only to generate
# random trees, `phangorn` only as a third-party likelihood oracle (optional).
#
# Run from the repository root:  Rscript dev/soft-sankoff/01-identity-check.R

source("tests/testthat/helper-soft-sankoff.R")
suppressPackageStartupMessages(library("ape"))

set.seed(1)

Report <- function(label, value, tolerance = 1e-10) {
  ok <- is.finite(value) && abs(value) < tolerance
  cat(sprintf("%-58s %12.3e  %s\n", label, value, if (ok) "PASS" else "FAIL"))
  invisible(ok)
}

RandomTree <- function(nTip) {
  tr <- ape::rtree(nTip, rooted = TRUE, br = NULL)
  tr[["edge"]]
}

results <- logical(0)

# ---------------------------------------------------------------------------
cat("\n== 1. softmin -> min as T -> 0 ==\n")
x <- c(3.2, 1.7, 1.7, 9.0)
for (temperature in c(1, 1e-1, 1e-2, 1e-3, 1e-4)) {
  cat(sprintf("  T = %-8g softmin = %.10f   (min = %.10f)\n",
              temperature, SoftMin(x, temperature), min(x)))
}
results <- c(results, Report("softmin(T=0) - min", SoftMin(x, 0) - min(x)))
cat(sprintf("  softmin <= min at every T: %s\n",
            all(vapply(10^seq(0, -6), function(tt) SoftMin(x, tt) <= min(x) + 1e-12,
                       logical(1)))))

# softmin does NOT converge to min at rate T; it converges to
# min - T * log(multiplicity of the minimum).  That gap is the whole point:
# it is the entropy of the near-optimal set, i.e. the degree to which the
# criterion integrates over reconstructions instead of optimising them.
# `x` has two tied minima, so (min - softmin) / T -> log(2).
cat("\n   multiplicity law: (min - softmin) / T -> log(#minima)\n")
for (temperature in c(1e-3, 1e-5, 1e-7)) {
  cat(sprintf("  T = %-8g (min - softmin)/T = %.10f   (log 2 = %.10f)\n",
              temperature, (min(x) - SoftMin(x, temperature)) / temperature,
              log(2)))
}
results <- c(results, Report("(min - softmin)/T - log(2)  at T = 1e-7",
                             (min(x) - SoftMin(x, 1e-7)) / 1e-7 - log(2), 1e-6))

# ---------------------------------------------------------------------------
cat("\n== 2. T -> 0 recovers hard Sankoff (random cost matrices) ==\n")
worstHard <- 0
for (rep in 1:20) {
  nTip <- sample(5:12, 1)
  nStates <- sample(2:5, 1)
  edge <- RandomTree(nTip)
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  cost <- matrix(round(runif(nStates^2, 0.5, 4), 2), nStates, nStates)
  diag(cost) <- 0
  soft <- SoftSankoffScore(edge, tc, cost, temperature = 1e-7)
  hard <- HardSankoffScore(edge, tc, cost)
  worstHard <- max(worstHard, abs(soft - hard))
}
results <- c(results, Report("max |soft(T=1e-7) - hardSankoff| over 20 trees",
                             worstHard, 1e-4))

# ---------------------------------------------------------------------------
cat("\n== 3. T -> 0 with equal costs recovers Fitch ==\n")
worstFitch <- 0
for (rep in 1:20) {
  nTip <- sample(5:12, 1)
  nStates <- sample(2:4, 1)
  edge <- RandomTree(nTip)
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0
  soft <- SoftSankoffScore(edge, tc, cost, temperature = 1e-7)
  fitch <- FitchScore(edge, states, nStates)
  worstFitch <- max(worstFitch, abs(soft - fitch))
}
results <- c(results, Report("max |soft(T=1e-7) - Fitch| over 20 trees",
                             worstFitch, 1e-4))

# ---------------------------------------------------------------------------
cat("\n== 4. T = 1 with cost = -log(P) recovers Felsenstein pruning ==\n")
worstFels <- 0
for (rep in 1:20) {
  nTip <- sample(5:12, 1)
  nStates <- sample(2:4, 1)
  tr <- ape::rtree(nTip, rooted = TRUE)
  edge <- tr[["edge"]]
  # Per-branch transition matrices, indexed by child node.
  P <- vector("list", max(edge))
  cost <- vector("list", max(edge))
  for (e in seq_len(nrow(edge))) {
    child <- edge[e, 2]
    P[[child]] <- MkTransition(nStates, tr[["edge.length"]][e])
    cost[[child]] <- -log(P[[child]])
  }
  states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
  tc <- TipCosts(states, nStates)
  tipLik <- exp(-tc)                      # 1 where observed, 0 elsewhere
  rootFreq <- rep(1 / nStates, nStates)

  soft <- SoftSankoffScore(edge, tc, cost, temperature = 1,
                           rootCost = -log(rootFreq))
  fels <- FelsensteinLogLik(edge, tipLik, P, rootFreq)
  worstFels <- max(worstFels, abs(soft - (-fels)))
}
results <- c(results, Report("max |soft(T=1) - (-logLik)| over 20 trees",
                             worstFels, 1e-10))

# ---------------------------------------------------------------------------
cat("\n== 5. Third-party oracle: phangorn::pml ==\n")
if (requireNamespace("phangorn", quietly = TRUE)) {
  nTip <- 8
  nStates <- 2
  tr <- ape::rtree(nTip, rooted = FALSE)
  tr[["edge.length"]] <- round(runif(nrow(tr[["edge"]]), 0.02, 0.4), 3)
  states <- sample(c("0", "1"), nTip, replace = TRUE)
  names(states) <- tr[["tip.label"]]
  dat <- phangorn::phyDat(as.matrix(states), type = "USER", levels = c("0", "1"))
  pmlFit <- phangorn::pml(tr, dat)

  # Match phangorn's convention: unrooted tree scored by rooting at a tip's
  # neighbour is equivalent for a reversible model; use ape's rooted form.
  trR <- ape::multi2di(ape::root(tr, outgroup = tr[["tip.label"]][1],
                                 resolve.root = TRUE))
  edge <- trR[["edge"]]
  P <- vector("list", max(edge))
  cost <- vector("list", max(edge))
  for (e in seq_len(nrow(edge))) {
    child <- edge[e, 2]
    bl <- trR[["edge.length"]][e]
    if (is.na(bl)) bl <- 0
    P[[child]] <- MkTransition(2, bl)
    cost[[child]] <- -log(P[[child]])
  }
  lookup <- match(trR[["tip.label"]], names(states))
  tc <- TipCosts(lapply(states[lookup], function(s) if (s == "0") 1L else 2L), 2)
  rootFreq <- c(0.5, 0.5)
  soft <- SoftSankoffScore(edge, tc, cost, temperature = 1,
                           rootCost = -log(rootFreq))
  cat(sprintf("  phangorn logLik = %.10f\n", as.numeric(pmlFit[["logLik"]])))
  cat(sprintf("  soft-Sankoff    = %.10f  (negated: %.10f)\n",
              soft, -soft))
  results <- c(results, Report("|(-soft) - phangorn logLik|",
                               -soft - as.numeric(pmlFit[["logLik"]]), 1e-8))
} else {
  cat("  phangorn not installed - skipped\n")
}

# ---------------------------------------------------------------------------
cat("\n== 6. Score is monotone non-increasing in T ==\n")
nTip <- 10
nStates <- 3
edge <- RandomTree(nTip)
states <- replicate(nTip, sample(seq_len(nStates), 1), simplify = FALSE)
tc <- TipCosts(states, nStates)
cost <- matrix(1, nStates, nStates)
diag(cost) <- 0
grid <- c(0, 10^seq(-4, 0.5, length.out = 12))
scores <- vapply(grid, function(tt) {
  SoftSankoffScore(edge, tc, cost, temperature = tt)
}, numeric(1))
for (i in seq_along(grid)) {
  cat(sprintf("  T = %-10.4g score = %.6f\n", grid[i], scores[i]))
}
results <- c(results, Report("max increase between consecutive T",
                             max(c(0, diff(scores))), 1e-9))

# ---------------------------------------------------------------------------
cat("\n== 7. Up-pass marginals ==\n")
marg <- SoftSankoffMarginals(edge, tc, cost, temperature = 0.5)
cat("  row sums (internal nodes):",
    sprintf("%.6f", rowSums(marg)[(nTip + 1):nrow(marg)]), "\n")
results <- c(results, Report("max |rowSum - 1| over all nodes",
                             max(abs(rowSums(marg) - 1)), 1e-10))

# ---------------------------------------------------------------------------
cat(sprintf("\n%d/%d checks passed\n", sum(results), length(results)))
if (!all(results)) quit(status = 1)
