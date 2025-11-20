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
#' @param char Vector of tokens listing states for the character in question.
#' @param ambiguousTokens Vector specifying which tokens, if any, correspond to
#' the ambiguous token (`?`).
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
#' @importFrom stats setNames
#' @importFrom TreeTools Log2Unrooted
#' @family profile parsimony functions
#' @export
StepInformation <- function (char, ambiguousTokens = c("-", "?")) {
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
  
  if (length(split) > 2L) {
    warning("Ignored least informative tokens where more than two informative ",
            "tokens present.")
    ranked <- order(order(split, decreasing = TRUE))
    split <- split[ranked < 3]
  }
  
  logProfile <- vapply(seq_len(split[2]), LogCarter1, double(1),
                       split[1], split[2])
  ret <- setNames(Log2Unrooted(sum(split[1:2]))
                  - (.LogCumSumExp(logProfile) / log(2)),
                  seq_len(split[2]) + sum(singletons))
  ret[ret < sqrt(.Machine[["double.eps"]])] <- 0 # Floating point error inevitable
  
  # Return:
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
#' Functions to compute the number of trees with a given parsimony score.
#' 
#' `Carter1()` calculates the number of trees in which Fitch parsimony will
#' reconstruct  _m_ steps, where _a_ leaves are labelled with one state,
#' and _b_ leaves are labelled with a second state, using theorem 1 of
#' \insertCite{Carter1990;textual}{TreeTools}
#' 
#' `MaddisonSlatkin()` generalises this result to characters with multiple
#' steps using the recursive approach of
#' \insertCite{Maddison1991;textual}{TreeSearch}.
#' 
#' @param m,steps Number of steps.
#' @param a,b Number of leaves labelled `0` and `1`.
#' @param states Number of leaves labelled with each possible state.
#' States are presented in binary fashion.  The first entry of the vector
#' corresponds to state `1` (binary `001`),
#' the second to state `2` (binary `010`),
#' and the third to the ambiguous state `01` (binary `011`).
#' 
#' 
#' @seealso [TreeTools::NUnrooted()]: number of unrooted trees with _n_ leaves.
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
#' Carter1(3, 3, 3)
#' NUnrooted(3 + 3)
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
  ret <- sum(
    log2(twoN - twoM - m),
    (lfactorial(m - 1L) / log(2)),
    Log2DoubleFactorial(twoN - 5L),
    Log2N(a, m),
    Log2N(b, m)
  ) - Log2DoubleFactorial(twoN - twoM - 1L)
  
  # Return:
  if (abs(ret) < sqrt(.Machine$double.eps)) 0 else ret
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
  ret <- sum(
    log(twoN - twoM - m),
    lfactorial(m - 1L),
    LogDoubleFactorial(twoN - 5L),
    LogN(a, m),
    LogN(b, m)
  ) - LogDoubleFactorial(twoN - twoM - 1L)
  
  # Return:
  if (abs(ret) < sqrt(.Machine$double.eps)) 0 else ret
}



# R is the probability that in a randomly selected tree of _n_ taxa, the
# smaller of the two basal subclades will have _m_ taxa.
# 
# The formula in Maddison & Slatkin 1991 overcounts by a factor of 2; rather
# than doubling the count where n != 2m (why?), we should halve the count
# where n == 2m: because choosing A B : C D gives the same trees as C D : A B,
# choose(2m, m) counts each tree twice.
# 
# Without this correction, we can see that R can exceed 1,
# e.g. R(1 | 3) = 2 * 4 * 1 * 3 / 15, using T = NRooted.
# (Note the misplaced brackets in Maddison & Slatkin 1991's definition of T(n))
# 
# TODO simplify by cancelling the choose(n,m) term with the one in .LogD
.LogR <- function(m, n) {
  log(if (n == m + m) 1 / 2 else 1) + 
    lchoose(n, m) +
    LnRooted(m) + LnRooted(n - m) - LnRooted(n)
}

# D is the probability that, in a randomly selected tree on `leaves`, the
# smaller subclade of taxa will receive taxa with labels `drawn`.
# If both subclades are the same size, the 'smaller' subclade is drawn at
# random
.LogD <- function(drawn, leaves) {
  sum(lchoose(leaves, drawn)) - lchoose(sum(leaves), sum(drawn))
}

#' @importFrom fastmap fastmap
.LogB_cache <- new.env(parent = emptyenv())

#' @importFrom stringi stri_paste
.HashLeaves <- function(v) {
  stri_paste(v, collapse = ",")
}

# B(b | tokens) is the probability that state b is reconstructed at the base
# of a clade with leaves labelled `tokens`
.LogB <- function(token, leaves, dp) {
  if (token < 1 || token > length(leaves)) {
    stop("`token` must be 1 .. length(leaves)")
  }
  
  leafHash <- .HashLeaves(leaves)
  sub <- .LogB_cache[[leafHash]]
  if (is.null(sub)) {
    sub <- new.env(parent = emptyenv())
    .LogB_cache[[leafHash]] <- sub
  }
  val <- sub[[as.character(token)]]
  if (!is.null(val)) {
    return(val)
  }
  
  if (missing(dp)) {
    nLevels <- floor(log2(length(leaves))) + 1
    nStates <- 2 ^ nLevels - 1
    dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  }
  n <- sum(leaves)
  result <- if (n == 1) {
    # If the one leaf we're looking at bears the token, p = 1; else p = 0
    log(leaves[[token]] == 1)
  } else if (n == 2) {
    result <- if (any(leaves == 2)) {
      which.max(leaves)
    } else {
      dp[rbind(which(leaves > 0))]
    }
    log(if (token == result) 1 else 0)
  } else {
    LogSumExp(apply(.ValidDraws(leaves), 1, function(drawn) {
      m <- sum(drawn)
      undrawn <- leaves - drawn
      # if both subtrees are equal size, D will give half the probability we want.
      balancedCorrection <- if (m + m == n) log(2) else log(1)
      if (all(drawn == undrawn)) {
        # Double counting in symmetrical case
        balancedCorrection <- balancedCorrection - log(2)
      }
      
      # Return:
      balancedCorrection + 
        .LogR(m, n) +
        .LogD(drawn, leaves) +
        LogSumExp(apply(
          which(dp == token, arr.ind = TRUE), 1, function(pair) {
            LogProdExp(list(.LogB(pair[[1]], drawn, dp),
                            .LogB(pair[[2]], undrawn, dp)))
          }))
    }))
  }
  sub[[as.character(token)]] <- result
  result
}


.ValidDraws <- function(leaves) {
  grid <- do.call(expand.grid, lapply(leaves, function(mc) seq.int(0, mc))) |>
    as.matrix()
  nDrawn <- rowSums(grid)
  n <- sum(leaves)
  evenSplits <- which(nDrawn == n / 2)
  duplicated <- evenSplits[
    apply(grid[evenSplits, , drop = FALSE], 1,
          function(drawn) {
            undrawn <- leaves - drawn
            cmp <- NA_integer_
            for (i in seq_along(drawn)) {
              if (drawn[[i]] < undrawn[[i]]) {cmp <- -1L; break }
              if (drawn[[i]] > undrawn[[i]]) {cmp <- +1L; break }
            }
            # Return:
            !is.na(cmp) && cmp == +1L
          })]
  validDraws <- setdiff(which(nDrawn > 0 & nDrawn <= floor(n / 2)), duplicated)
  grid[validDraws, , drop = FALSE]
}

#' @rdname Carter1
#' @examples
#' 
#' LogCarter1(1, 8, 24)
#' LogCarter1(2, 8, 24)
#' LogCarter1(3, 8, 24)
#' 
#' MaddisonSlatkin(1, c(2, 2))
#' MaddisonSlatkin(1, c(1, 1))
#' LogCarter1(1,2,2)
#' 
#' @importFrom TreeTools LnRooted
#' @export
MaddisonSlatkin <- function(steps, states) {
  nTaxa <- sum(states)
  nLevels <- floor(log2(length(states))) + 1
  nStates <- 2 ^ nLevels - 1
  length(states) <- nStates
  states[is.na(states)] <- 0
  dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  
  # Probability of seeing s steps given tokens at leaves, and `token` at root
  .LogP <- function(s, leaves, token) {
    n <- sum(leaves)
    if (n == 1) {
      if(leaves[[token]] != 1) {
        # Impossible to have `token` at base given `leaves`
        return(log(0))
      }
      return(log(s == 0))
    }
    if (n == 2) {
      return(log(if (any(leaves == 2)) {
        s == 0
      } else {
        s == attr(dp, "step")[rbind(which(leaves > 0))]
      }))
    }
    
    denominator <- .LogB(token, leaves, dp)
    if (is.finite(denominator)) {
      
      apply(.ValidDraws(leaves), 1, function(drawn) {
        m <- sum(drawn)
        undrawn <- leaves - drawn
        
        # Return:
        .LogR(m, n) +
          .LogD(drawn, leaves) +
          LogSumExp(
            
            
            # 0:s where the conjunction doesn't add a step
            vapply(0:s, function(r) {
              apply(which(dp == token & !attr(dp, "step"), arr.ind = TRUE), 1, function(pair)
                LogProdExp(list(
                  .LogP(r, drawn, pair[[1]]),
                  .LogB(pair[[1]], drawn, dp),
                  .LogP(s - r, undrawn, pair[[2]]),
                  .LogB(pair[[2]], undrawn, dp)
                ))
              ) |> LogSumExp()
            }, double(1)),
            
            # 0:(s - 1) where the conjunction adds a step
            vapply(seq_len(s) - 1, function(r) {
              which(dp == token & attr(dp, "step"), arr.ind = TRUE) |>
                # Empty call to which will deliver pair = c(0, 0)
                apply(1, function(pair) if (pair[[1]] == 0) -Inf else {
                  LogProdExp(list(
                    .LogP(r, drawn, pair[[1]]),
                    .LogB(pair[[1]], drawn, dp),
                    # s - r - 1 because we're adding a step here
                    .LogP(s - r - 1, undrawn, pair[[2]]),
                    .LogB(pair[[2]], undrawn, dp)
                  ))
                }) |> LogSumExp()
              }, double(1))
          )
      }) |> LogSumExp() - denominator
    } else {
      denominator
    }
  }
  
  
  p <- log(0)
  for (state in seq_len(nStates)) {
    p <- LogSumExp(p, LogProdExp(list(
      .LogB(state, states, dp),
      .LogP(steps, states, state))
      ))
  }
  p
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
    # TODO test splits <- 2 2 4
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
        # TODO would be quicker to calculate just first half; special case:
        #  omitted.tips %% 2
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
