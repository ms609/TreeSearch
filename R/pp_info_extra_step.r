#' Information content of a character known to contain _e_ steps
#'
#' `StepInformation()` calculates the phylogenetic information content of a
#' character `char` when _e_ extra steps are present, for all possible 
#' values of _e_.
#' 
#' Calculates the number of trees consistent with the character having 
#' _e_ extra steps, where _e_ ranges from its minimum possible value
#' (i.e. number of different tokens minus one) to its maximum.
#'
#' For characters with 2 informative tokens, uses the exact formula of
#' Carter _et al._ (1990) via [LogCarter1()].
#' For characters with 3--5 informative tokens, uses the recursive algorithm
#' of Maddison & Slatkin (1991) via [MaddisonSlatkin()].
#' Characters with more than 5 informative tokens are reduced to their 5 most
#' frequent tokens (with a warning).
#'
#' When the Maddison & Slatkin computation would be infeasible (exponential
#' in the number of tips for a given number of tokens), behaviour depends on
#' the `approx` argument.  With `"auto"` (default), infeasible characters are
#' reduced to their 2 most frequent tokens and scored using the faster Carter
#' formula.  This occurs for 3-state characters with more than 27 informative
#' tips, 4-state characters with more than 19 tips, and 5-state characters with
#' more than 13 tips.  With `"mc"`, the MC-calibrated approximation is used
#' instead: the step-count distribution is estimated from a random sample of
#' trees, anchored at the exact minimum-steps probability.
#'
#' @param char Vector of tokens listing states for the character in question.
#' @param ambiguousTokens Vector specifying which tokens, if any, correspond to
#' the ambiguous token (`?`).
#' @param approx Character string controlling the fallback for infeasible
#'   multi-state characters: `"auto"` (default) reduces to the 2 most frequent
#'   tokens; `"mc"` uses a Monte Carlo approximation that preserves the full
#'   multi-state character (see Details); `"exact"` forces exact computation
#'   regardless of cost (may be very slow for large characters).
#' @param n_mc Integer.  Number of random trees used by the `"mc"` approximation.
#'   Larger values improve accuracy but increase computation time.
#'   Default: 5000.
#' 
#' @return `StepInformation()` returns a numeric vector detailing the amount
#' of phylogenetic information (in bits) associated with the character when
#' 0, 1, 2&hellip; extra steps are present.  The vector is named with the
#' total number of steps associated with each entry in the vector: for example,
#' a character with three observed tokens must exhibit two steps, so the first
#' entry (zero extra steps) is named `2` (two steps observed).
#' 
#' @examples
#' character <- rep(c(0:3, "?", "-"), c(8, 5, 1, 1, 2, 2))
#' StepInformation(character)
#' @template MRS
#' @importFrom fastmatch %fin%
#' @importFrom stats setNames dnorm sd
#' @importFrom TreeTools Log2Unrooted LnUnrooted MatrixToPhyDat NUnrooted
#'   NUnrootedMult RandomTree
#' @family profile parsimony functions
#' @export
StepInformation <- function (char, ambiguousTokens = c("-", "?"),
                             approx = "auto", n_mc = 5000L) {
  NIL <- c("0" = 0)
  char <- char[!char %fin% ambiguousTokens]
  if (length(char) == 0) {
    return(NIL)
  }
  split <- sort(as.integer(table(char)), decreasing = TRUE)
  singletons <- split == 1L
  nSingletons <- sum(singletons)
  nInformative <- sum(split) - nSingletons
  minSteps <- length(split) - 1L
  if (minSteps == 0L) {
    return(NIL)
  }
  
  split <- split[!singletons]
  if (length(split) < 2L) {
    return(setNames(0, minSteps))
  }
  
  k <- length(split)
  
  if (k > 5L) {
    warning("More than 5 informative tokens; keeping the 5 most frequent.")
    split <- split[seq_len(5L)]
    k <- 5L
  }
  
  nTips <- sum(split)
  
  # Feasibility check: use partition-aware split_count rather than n alone.
  # Balanced partitions blowup at lower n than skewed ones; the split_count
  # (coeff of x^floor(n/2) in prod_i(1+x+...+x^{a_i})) captures this.
  # Thresholds from .MS_SC_THRESHOLD in data_manipulation.R.
  infeasible <- k >= 3L &&
    .MSSplitCount(split) > .MS_SC_THRESHOLD[k]
  
  if (infeasible && identical(approx, "mc")) {
    # MC-calibrated normal approximation: preserves multi-state character
    return(.ApproxStepInformation(split, n_mc = n_mc, nSingletons = nSingletons))
  }
  
  if (infeasible && !identical(approx, "exact")) {
    warning(
      "Exact multi-state profile computation infeasible for ",
      k, "-state character with ", nTips, " tips; ",
      "reducing to 2 most frequent tokens."
    )
    split <- split[seq_len(2L)]
    k <- 2L
    nTips <- sum(split)
  }
  
  if (k == 2L) {
    # Binary: use Carter (fast, exact)
    logProfile <- vapply(seq_len(split[2]), LogCarter1, double(1),
                         split[1], split[2])
    # Convert log-count to log-probability
    logP <- logProfile - LnUnrooted(nTips)
    reducedMinSteps <- 1L
  } else {
    # Multi-state (3-5): use MaddisonSlatkin
    nStates <- 2L^k - 1L
    states <- integer(nStates)
    for (i in seq_along(split)) {
      states[2L^(i - 1L)] <- split[i]
    }
    reducedMinSteps <- k - 1L
    maxSteps <- nTips - 1L
    logP <- MaddisonSlatkin(reducedMinSteps:maxSteps, states)
  }
  
  # Trim trailing -Inf entries (impossible step counts)
  finite_idx <- which(is.finite(logP))
  if (length(finite_idx) == 0L) {
    return(setNames(0, minSteps))
  }
  logP <- logP[seq_len(max(finite_idx))]
  
  # Cumulative information: -log2(cumsum(P))
  ret <- -.LogCumSumExp(logP) / log(2)
  
  # Name with total step counts (reduced steps + singleton offset)
  names(ret) <- seq.int(reducedMinSteps,
                         reducedMinSteps + length(ret) - 1L) + nSingletons
  
  ret[ret < sqrt(.Machine[["double.eps"]])] <- 0 # Floating point error inevitable
  
  # Return:
  ret
}

# MC-calibrated normal approximation for infeasible multi-state characters.
# Returns a named IC vector matching the format of StepInformation().
#
# @param split Integer vector of informative token frequencies (sorted
#   decreasing, singletons removed).
# @param n_mc Integer. Number of Monte Carlo trees to score.
# @param nSingletons Integer. Number of singleton tokens (for step offset).
# @return Named numeric vector of IC (bits) by step count.
# @keywords internal
.ApproxStepInformation <- function(split, n_mc = 5000L, nSingletons = 0L) {
  k <- length(split)
  n <- sum(split)
  s_min <- k - 1L
  s_max <- n - 1L
  
  # 1. Exact P(s_min) via product-of-double-factorials formula O(k)
  log_p_min <- log(NUnrootedMult(split)) - log(NUnrooted(n))
  
  # 2. MC: sample random trees and score this character with Fitch parsimony
  labels  <- paste0("t", seq_len(n))
  char_vec <- rep(seq_len(k) - 1L, split)
  names(char_vec) <- labels
  dat <- MatrixToPhyDat(matrix(char_vec, ncol = 1L,
                                dimnames = list(labels, "c1")))
  
  trees <- lapply(seq_len(n_mc), function(i) {
    TreeTools::RootTree(RandomTree(labels), 1L)
  })
  class(trees) <- "multiPhylo"
  mc_scores <- TreeLength(trees, dat)
  
  # 3. Fit normal to MC body (for interpolating the near-s_min gap)
  mu_hat <- mean(mc_scores)
  sd_hat <- sd(mc_scores)
  
  # 4. Build log-probability vector for s_min..s_max
  steps  <- s_min:s_max
  log_p  <- numeric(length(steps))
  
  for (i in seq_along(steps)) {
    s <- steps[[i]]
    if (s == s_min) {
      log_p[[i]] <- log_p_min
    } else {
      mc_count <- sum(mc_scores == s)
      if (mc_count > 0L) {
        log_p[[i]] <- log(mc_count / n_mc)
      } else {
        # Normal extrapolation for step counts not reached by MC
        log_p[[i]] <- dnorm(s, mu_hat, sd_hat, log = TRUE)
      }
    }
  }
  
  # 5. Trim trailing negligible entries
  finite_idx <- which(is.finite(log_p) & log_p > -700)
  if (length(finite_idx) == 0L) {
    return(setNames(0, s_min + nSingletons))
  }
  log_p <- log_p[seq_len(max(finite_idx))]
  steps  <- steps[seq_len(max(finite_idx))]
  
  # 6. Cumulative IC
  ret <- -.LogCumSumExp(log_p) / log(2)
  names(ret) <- steps + nSingletons
  ret[ret < sqrt(.Machine[["double.eps"]])] <- 0
  
  ret
}

# Adapted from https://rpubs.com/FJRubio/LSE
.LogCumSumExp <- function (x) { 
  n <- length(x)
  Lk <- c(x[1], double(n - 1L))
  for (k in 1L + seq_len(n - 1L)) {
    Lk[k] <- Lk[k - 1]
    Lk[k] <- max(x[k], Lk[k]) + log1p(exp(-abs(x[k] - Lk[k])))
  }
  
  # Return:
  Lk
}

#' Number of trees with _m_ steps
#' 
#' Calculate the number of unrooted binary trees on which Fitch parsimony
#' reconstructs exactly _m_ steps for a character.
#' 
#' `Carter1()` (and its logarithmic variants `Log2Carter1()`, `LogCarter1()`)
#' implement theorem 1 of \insertCite{Carter1990;textual}{TreeTools} for
#' **binary** characters, where _a_ leaves bear one state and _b_ bear the
#' other.
#' 
#' `MaddisonSlatkin()` generalises this result to characters with multiple
#' states using the recursive approach of
#' \insertCite{Maddison1991;textual}{TreeSearch}.
#' It returns the **log-probability** (i.e. log of the fraction of unrooted
#' binary trees) for each requested step count.  Currently supports 2--5
#' character tokens.
#' 
#' @param m,steps Number of steps.
#' @param a,b Number of leaves labelled `0` and `1`.
#' @param states Integer vector giving the number of leaves bearing each
#'   possible combination of states, laid out in binary fashion.
#'   Entry 1 = state `1` (binary `001`), entry 2 = state `2` (binary `010`),
#'   entry 3 = ambiguous state `{1,2}` (binary `011`), and so on.
#'   Only observed singleton states need non-zero counts; polymorphic entries
#'   are typically zero.
#' 
#' @references 
#' \insertAllCited{}
#' 
#' See also:
#' 
#' \insertRef{Steel1993}{TreeSearch}
#' 
#' \insertRef{Steel1995}{TreeSearch}
#' 
#' (\insertRef{Steel1996}{TreeSearch})
#' @importFrom TreeTools LogDoubleFactorial
#' @examples 
#' # The character `0 0 0 1 1 1`
#' Carter1(1, 3, 3) # Exactly one step
#' Carter1(2, 3, 3) # Two steps (one extra step)
#' 
#' # Number of trees that the character can map onto with exactly _m_ steps
#' # if non-parsimonious reconstructions are permitted:
#' cumsum(sapply(1:3, Carter1, 3, 3))
#' 
#' # Three steps allow the character to map onto any of the 105 six-leaf trees.
#' 
#' @template MRS
#' @family profile parsimony functions
#' @export
Carter1 <- function (m, a, b) {
  n <- a + b
  twoN <- n + n
  twoM <- m + m
  N <- function (n, m) {
    if (n < m) {
      0
    } else {
      nMinusM <- n - m
      factorial(n + nMinusM - 1L) / prod(
        factorial(nMinusM),
        factorial(m - 1L),
        2 ^ (nMinusM))
    }
  }
  prod(
    (twoN - twoM - m),
    N(a, m),
    N(b, m),
    exp(lfactorial(m - 1L) + 
          LogDoubleFactorial(twoN - 5L) -
          LogDoubleFactorial(twoN - twoM - 1L))
    )
}

#' @rdname Carter1
#' @export
#' @importFrom TreeTools Log2DoubleFactorial
Log2Carter1 <- function (m, a, b) {
  n <- a + b
  twoN <- n + n
  twoM <- m + m
  Log2N <- function (n, m) {
    if (n < m) -Inf else {
      nMinusM <- n - m
      (lfactorial(n + nMinusM - 1L) - 
       lfactorial(nMinusM) -
       lfactorial(m - 1L)) / log(2) - nMinusM
    }
  }
  sum(
    log2(twoN - twoM - m),
    (lfactorial(m - 1L) / log(2)),
    Log2DoubleFactorial(twoN - 5L),
    Log2N(a, m),
    Log2N(b, m)
  ) - Log2DoubleFactorial(twoN - twoM - 1L)
}

#' @rdname Carter1
#' @export
#' @importFrom TreeTools LogDoubleFactorial
LogCarter1 <- function (m, a, b) {
  n <- a + b
  twoN <- n + n
  twoM <- m + m
  LogN <- function (n, m) {
    if (n < m) -Inf else {
      nMinusM <- n - m
      (lfactorial(n + nMinusM - 1L) - 
       lfactorial(nMinusM) -
       lfactorial(m - 1L)) - (nMinusM * log(2))
    }
  }
  sum(
    log(twoN - twoM - m),
    lfactorial(m - 1L),
    LogDoubleFactorial(twoN - 5L),
    LogN(a, m),
    LogN(b, m)
  ) - LogDoubleFactorial(twoN - twoM - 1L)
}

# TODO: Replace the below with an advanced version of Maddison & Slakey 1991, 
# or use the results of Carter et al. 1990; Steel 1993 to estimate +0 & +1 steps,
# and approximate the rest.

## @importFrom TreeTools Log2UnrootedMult
# Old_IC_Approx <- function() {
#   
#   nIter <- min(maxIter, round(iter))
#   if (nIter == maxIter) {
#     warning ("Will truncate number of iterations at maxIter = ", maxIter)
#   }
#   n01ExtraSteps <- nOneExtraStep + nNoExtraSteps
#   analyticIC <- Log2Unrooted(sum(split)) -  setNames(c(
#     Log2UnrootedMult(split), log2(n01ExtraSteps)),
#     minSteps + 0:1)
#   analyticP <- 2 ^ -analyticIC[2]
#   
#   if (warn) {
#     message("  Token count ", split, " = ",
#             signif(analyticIc0, ceiling(log10(maxIter))),
#             " bits @ 0 extra steps. \n  Simulating ", nIter, 
#             " trees to estimate cost of further steps.")
#     # message(c(round(analyticIc0, 3), "bits @ 0 extra steps;", round(analyticIc1, 3),
#     #    "@ 1; attempting", nIter, "iterations.\n"))
#   }
#   
#   morphyObj <- SingleCharMorphy(rep(seq_along(split) - 1L, split))
#   on.exit(morphyObj <- UnloadMorphy(morphyObj))
#   steps <- vapply(rep(nInformative, nIter), RandomTreeScore,
#                   integer(1), morphyObj) + nSingletons
#   
#   tabSteps <- table(steps[steps > (minSteps - nSingletons + 1)]) # Quicker than table(steps)[-1]
#   
#   approxP <- tabSteps / sum(tabSteps) * (1 - analyticP)
#   approxSE <- sqrt(approxP * (1 - approxP) / nIter)
#   cumP <- cumsum(c(analyticP, approxP))[-1]
# 
#   approxIC <- -log2(cumP)
#   icLB <- -log2(cumP - approxSE)
#   icError <- icLB - approxIC
#   if (warn || max(icError) > tolerance) {
#     message("  Approx. std. error < ", signif(max(icError) * 1.01, 2))
#   }
# }


#' Number of trees with one extra step
#' @param \dots Vector or series of integers specifying the number of leaves
#' bearing each distinct non-ambiguous token.
#' @importFrom TreeTools NRooted NUnrooted
#' @examples
#' WithOneExtraStep(1, 2, 3)
#' @importFrom TreeTools NUnrootedMult DoubleFactorial
#' @family profile parsimony functions
#' @export
WithOneExtraStep <- function (...) {
  splits <- c(...)
  # Ignore singletons, which can be added at the end...
  singletonSplits <- splits == 1
  splits <- sort(splits[!singletonSplits], decreasing = TRUE)
  nSplits <- length(splits)
  if (nSplits < 2) return (0)
  if (nSplits == 2) {
    prod(
      # Zone 0; Zone 1 will be added at Root
      NRooted(splits[1]),
      sum(
        vapply(seq_len(splits[2] - 1L), function (beforeStep) {
          NRooted(beforeStep) * # Zone 1 will sit at root of Zone 0
            sum(
             # Case 1: Zone 1 & 2 not adjacent
             (splits[1] + splits[1] - 4L) * # Edges not touching Zone 1
               DoubleFactorial(splits[2] + splits[2] - 4L) /
               DoubleFactorial(beforeStep + beforeStep - 2L),
             # Case 2: Zone 1 & Zone 2 adjacent
             2 * # Two edges adjacent to Zone 1
               DoubleFactorial(splits[2] + splits[2] - 1L) /
               DoubleFactorial(beforeStep + beforeStep + 1L)
            )
        }, double(1))
        
      ),
      # Add singleton splits
      (2 * (sum(splits) + seq_len(sum(singletonSplits))) - 5)
    )
    
  } else {
    
    stop("Not implemented.")
                                                                                # nocov start
    sum(vapply(seq_along(splits), function (omit) {
      backboneSplits <- splits[-omit]
      omitted.tips <- splits[omit]
      backbone.tips <- sum(backboneSplits)
      backbones <- NUnrootedMult(backboneSplits)
      backbone.edges <- max(0L, 2L * backbone.tips - 3L)
      attachTwoRegions <- backbone.edges * (backbone.edges - 1L)
      
      
      prod( # omitted tips form two separate regions
        backbones,
        attachTwoRegions,
        sum(
        vapply(seq_len(omitted.tips - 1), function (first.group) { 
          # For each way of splitsting up the omitted tips, e.g. 1|16, 2|15, 3|14, etc
          choose(omitted.tips, first.group) * 
            NRooted(first.group) * NRooted(omitted.tips - first.group)
        }, double(1))
      ) / 2) +
        prod(
          # paraphyletic.  Worry: This is equivalent to splitting gp. 0
          # Double count: (0, 0, (0, (1, (0, 1))
        )
  
    }, double(1))
    )
                                                                                # nocov end
  }
}

#' Clear `MaddisonSlatkin()` cache
#'
#' Releases the internal C++ cache used by `MaddisonSlatkin()`.
#' Needed only in testing or if memory pressure is a concern.
#'
#' @name MaddisonSlatkin_clear_cache
#' @keywords internal
#' @export
NULL
