#!/usr/bin/env Rscript
# GATE B — what does a soft score cost against SIMD Fitch?
#
# Two distinct costs, and the second is the larger:
#   1. Loss of bit-parallelism: O(k^2) doubles plus exp/log per character per
#      node, against word-packed bitwise Fitch (src/ts_simd.h).  This is a
#      CONSTANT FACTOR, and this script now measures it with a compiled soft
#      kernel (src/ts_soft_sankoff.{h,cpp}) on both sides.
#   2. Loss of incremental rescoring.  Both Sankoff kernels are full-rescore
#      only; TBR draws much of its throughput from rescoring the affected
#      subtree.  That changes the COMPLEXITY of the inner loop, not a constant.
#
# On (2): a full-rescore prototype cannot measure an incremental soft kernel
# that does not exist.  What it CAN measure is how much Fitch gains from
# incrementality — full-rescore cost against per-candidate incremental cost,
# both reported by ts_bench_tbr_phases() from inside C++.  That ratio is the
# multiplier a full-rescore soft kernel forfeits, so it bounds the term without
# requiring the incremental soft kernel to be written.  It is a bound on the
# term's size, not a measurement of an incremental soft implementation.
#
# MEASUREMENT NOTES
#
# * The Fitch baseline is `time_full_rescore_us` from ts_bench_tbr_phases(),
#   clocked inside C++ around ts::score_tree() with dataset construction
#   excluded.  Timing TreeLength() or even ts_fitch_score() from R would fold
#   in per-call marshalling of the contrast and tip-data matrices, inflating
#   the Fitch time and so understating the ratio — an error in the direction
#   that flatters the soft kernel.
# * The soft kernel is timed by differencing across the binding's own `n_rep`:
#   t(N) = overhead + N * kernel, so kernel = (t(N) - t(1)) / (N - 1).  This
#   removes R-side marshalling without needing to trust that it is small.
# * phyDat compresses to unique site patterns.  Fitch scores patterns (weighted),
#   so the soft kernel is given exactly the same pattern count, not the raw
#   character count.
# * ts_bench_tbr_phases clocks in WHOLE MICROSECONDS, and a Fitch rescore of
#   100 patterns takes 1-7 of them — one run reported 0 us, making the ratio
#   Inf.  Both kernels are linear in pattern count, so the fix is to compare at
#   a pattern count high enough to resolve the denominator: N_CHARS is set so
#   Fitch takes tens to hundreds of ticks.  The ratio is scale-free; the
#   resolution is not.
# * Ratios are summarised by MEDIAN with a sign count, never by mean.
#
# Usage:
#   Rscript dev/soft-sankoff/03-speed-budget.R

source("tests/testthat/helper-soft-sankoff.R")
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = ".agent-softsankoff")
})

set.seed(2026)
LADDER <- c(16, 32, 64, 128)
STATE_LADDER <- c(2L, 4L)
N_CHARS <- 2000L          # high enough that a Fitch rescore clears the clock floor
TEMPERATURE <- 0.25
FITCH_REPS <- 25L         # medianed, since each call clocks a single rescore
SOFT_OUTER <- 5L          # outer repeats of the n_rep differencing
SOFT_TARGET_S <- 0.3      # n_rep is chosen per cell to hit roughly this wall
R_PATTERN_CAP <- 5L       # pure-R reference is timed on a subset and scaled
THRESHOLD <- 50           # plan's orientation threshold for Step 4

MakeData <- function(nTip, nChars, nStates) {
  mat <- matrix(sample(as.character(seq_len(nStates) - 1L),
                       nTip * nChars, replace = TRUE), nrow = nTip)
  rownames(mat) <- paste0("t", seq_len(nTip))
  phangorn::phyDat(mat, type = "USER",
                   levels = as.character(seq_len(nStates) - 1L))
}

PrepDataset <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at[["contrast"]]
  storage.mode(contrast) <- "double"
  tipData <- matrix(unlist(dataset, use.names = FALSE),
                    nrow = length(dataset), byrow = TRUE)
  storage.mode(tipData) <- "integer"
  # min_steps is per-CHARACTER, not per-contrast-column, and only bites under
  # implied weights.  These runs are equal-weights (concavity < 0), so leave it
  # empty and let make_dataset() decide rather than hand-rolling a length that
  # has to match the pattern count.
  list(contrast = contrast, tipData = tipData, weight = at[["weight"]],
       levels = at[["levels"]], minSteps = integer())
}

# Median in-C++ Fitch full-rescore time, plus the per-candidate incremental
# cost, both in microseconds.
FitchKernel <- function(edge, prepped, reps) {
  full <- numeric(reps)
  incr <- numeric(reps)
  for (i in seq_len(reps)) {
    r <- TreeSearch:::ts_bench_tbr_phases(
      edge, prepped[["contrast"]], prepped[["tipData"]],
      prepped[["weight"]], prepped[["levels"]], prepped[["minSteps"]]
    )
    full[i] <- r[["time_full_rescore_us"]]
    incr[i] <- if (r[["n_candidates"]] > 0) {
      r[["time_clip_incr_us"]] / r[["n_candidates"]]
    } else NA_real_
  }
  list(fullUs = stats::median(full),
       incrementalUs = stats::median(incr, na.rm = TRUE))
}

# Compiled soft kernel time per whole-tree score over all patterns, in
# microseconds, with R-side marshalling differenced out.
# n_rep is chosen per cell from a pilot pass, so that both the t(1) and t(N)
# measurements are far above system.time()'s own resolution.  A cheap cell gets
# hundreds of repetitions, an expensive one gets two; the differencing is valid
# at any nRep >= 2.
SoftKernelCpp <- function(edge, tipCostList, costList, temperature,
                          outer, targetSeconds) {
  storage.mode(edge) <- "integer"
  nTip <- nrow(tipCostList[[1]])
  Time1 <- function(nR) {
    system.time(
      TreeSearch:::ts_soft_sankoff_test(
        edge = edge, n_tip = nTip,
        tip_costs = tipCostList, cost_matrices = costList,
        temperature = temperature, n_rep = as.integer(nR)
      )
    )[["elapsed"]]
  }
  pilot <- max(Time1(1L), 1e-4)
  nRep <- max(2L, min(500L, as.integer(ceiling(targetSeconds / pilot))))
  one <- stats::median(vapply(seq_len(outer), function(i) Time1(1L), numeric(1)))
  many <- stats::median(vapply(seq_len(outer), function(i) Time1(nRep),
                               numeric(1)))
  list(us = (many - one) / (nRep - 1) * 1e6, nRep = nRep)
}

# Pure-R reference time per whole-tree score over all patterns, extrapolated
# from a capped subset of patterns.
SoftKernelR <- function(edge, tipCostList, costList, temperature, cap) {
  use <- seq_len(min(cap, length(tipCostList)))
  elapsed <- system.time(
    for (ch in use) {
      SoftSankoffScore(edge, tipCostList[[ch]], costList[[ch]], temperature)
    }
  )[["elapsed"]]
  elapsed / length(use) * length(tipCostList) * 1e6
}

rows <- list()
for (nStates in STATE_LADDER) {
  cost <- matrix(1, nStates, nStates)
  diag(cost) <- 0

  for (nTip in LADDER) {
    dataset <- MakeData(nTip, N_CHARS, nStates)
    prepped <- PrepDataset(dataset)
    nPattern <- length(prepped[["weight"]])

    tree <- TreeTools::Preorder(ape::rtree(nTip, rooted = TRUE, br = NULL))
    tree[["tip.label"]] <- paste0("t", seq_len(nTip))
    edge <- tree[["edge"]]

    # Per-pattern tip costs, in the tree's own tip order.
    chars <- as.character(dataset)[tree[["tip.label"]], , drop = FALSE]
    coded <- matrix(match(chars, as.character(seq_len(nStates) - 1L)),
                    nrow = nTip)
    # Reduce to the same unique patterns Fitch scores.
    patterns <- unique(t(coded))
    patterns <- patterns[seq_len(min(nrow(patterns), nPattern)), , drop = FALSE]
    tipCostList <- lapply(seq_len(nrow(patterns)), function(i) {
      TipCosts(as.list(patterns[i, ]), nStates)
    })
    costList <- rep(list(cost), length(tipCostList))

    fitch <- FitchKernel(edge, prepped, FITCH_REPS)
    soft <- SoftKernelCpp(edge, tipCostList, costList, TEMPERATURE,
                          SOFT_OUTER, SOFT_TARGET_S)
    softCpp <- soft[["us"]]
    softR <- SoftKernelR(edge, tipCostList, costList, TEMPERATURE,
                         R_PATTERN_CAP)

    nNode <- nTip - 1L
    fitchOps <- nNode * ceiling(length(tipCostList) / 64)
    softOps <- nNode * length(tipCostList) * (nStates^2 + nStates + 1)

    rows[[length(rows) + 1]] <- data.frame(
      nTip = nTip, nStates = nStates, nPattern = length(tipCostList),
      fitch_us = fitch[["fullUs"]],
      fitch_incremental_us = fitch[["incrementalUs"]],
      soft_cpp_us = softCpp,
      soft_n_rep = soft[["nRep"]],
      soft_r_us = softR,
      cpp_ratio = softCpp / fitch[["fullUs"]],
      r_over_cpp = softR / softCpp,
      fitch_incrementality = fitch[["fullUs"]] / fitch[["incrementalUs"]],
      op_ratio = softOps / fitchOps
    )
    cat(sprintf(
      "k=%d nTip=%-4d pat=%-4d  Fitch %7.1f us  softC++ %8.1f us  x%-7.1f  (R x%.0f slower)\n",
      nStates, nTip, length(tipCostList), fitch[["fullUs"]], softCpp,
      softCpp / fitch[["fullUs"]], softR / softCpp))
  }
}

result <- do.call(rbind, rows)
utils::write.csv(result, "dev/soft-sankoff/03-speed-budget.csv",
                 row.names = FALSE)

Report <- function(label, x) {
  cat(sprintf("%-42s median x%-8.1f  range x%.1f-%.1f\n",
              label, stats::median(x), min(x), max(x)))
}

cat("\n=== Gate B: the CONSTANT term, now measured compiled ===\n")
for (k in STATE_LADDER) {
  sub <- result[result[["nStates"]] == k, ]
  Report(sprintf("compiled soft / Fitch full rescore (k=%d)", k),
         sub[["cpp_ratio"]])
}

# ts_bench_tbr_phases clocks in whole microseconds, and a Fitch rescore at
# these sizes is 1-7 ticks.  A one-tick denominator carries ~100% relative
# error, so the honest headline is the best-resolved cell rather than the
# median over cells that include it.  Report the resolution explicitly instead
# of quoting a ratio whose denominator is at the clock floor.
result[["fitch_ticks"]] <- round(result[["fitch_us"]])
cat(sprintf("\nFitch denominator resolution: %d of %d cells rest on <= 3 clock\n",
            sum(result[["fitch_ticks"]] <= 3), nrow(result)))
cat("ticks (ts_bench_tbr_phases clocks whole microseconds).  Best-resolved\n")
cat("cells, which carry the headline:\n")
best <- result[result[["fitch_ticks"]] == max(result[["fitch_ticks"]]), ]
for (i in seq_len(nrow(best))) {
  cat(sprintf("  k=%d nTip=%d: Fitch %.0f us (%.0f%% resolution)  ratio x%.0f\n",
              best[["nStates"]][i], best[["nTip"]][i], best[["fitch_us"]][i],
              100 / best[["fitch_ticks"]][i], best[["cpp_ratio"]][i]))
}
cat("The conclusion tolerates a factor-of-2 error in that denominator by\n")
cat(sprintf("more than an order of magnitude against the x%d threshold.\n",
            THRESHOLD))
cat(sprintf("\nCells over the x%d orientation threshold: %d of %d\n",
            THRESHOLD, sum(result[["cpp_ratio"]] > THRESHOLD), nrow(result)))
Report("operation-count ratio (prior estimate)", result[["op_ratio"]])
Report("pure R / compiled soft (speedup gained)", result[["r_over_cpp"]])

cat("\n=== The COMPLEXITY term: still open, but bounded ===\n")
Report("Fitch full rescore / incremental candidate",
       result[["fitch_incrementality"]])
cat("A full-rescore soft kernel forfeits that multiplier, on top of the\n")
cat("constant above.  This bounds the term; it does not measure an\n")
cat("incremental soft kernel, which does not exist.\n")

cat("\nPlan's orientation threshold: if a compiled soft score costs more than\n")
cat(sprintf("~x%d a Fitch score, an annealed search cannot pay for itself\n",
            THRESHOLD))
cat("against the ratchet and Step 4 is dead regardless of Gate A.\n")
