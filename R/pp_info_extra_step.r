#' Information Content Steps
#'
#' `ICSteps()` estimates the phylogenetic information content of a character 
#' `char` when _e_ extra steps are present, for all possible values of _e_.
#' 
#' Calculates the number of trees consistent with the character having 
#' _e_ extra steps, where _e_ ranges from its minimum possible value
#' (i.e. number of different tokens minus one) to its maximum.
#' The number of trees with no extra steps can be calculated exactly; the 
#' number of trees with more additional steps must be approximated.
#' The function samples `maxIter` trees, or enough trees that the trees with 
#' the minimum number of steps will be recovered at least `expectedMinima`
#' times, in order to obtain precise results.
#'
#' @param char Vector of tokens listing states for the character in question.
#' @param ambiguousToken Which token, if any, corresponds to the ambiguous token
#' (?) (not yet fully implemented).
#' @param tolerance Numeric specifying an anticipated upper bound on the
#' standard error of approximate profile scores.
#' @param maxIter Integer specifying maximum number of trees to score when
#' approximating profile scores.  Small values may result in errors
#' exceeding `tolerance`.
#' @template warnParam
#' 
#' @template MRS
#' @references
#'  \insertRef{Faith2001}{TreeSearch}
#' 
#' @keywords tree
#' 
#' @examples
#' character <- rep(c(0:3, '?'), c(8, 5, 1, 1, 2))
#' ICSteps(character, warn = FALSE)
#' 
#' @importFrom TreeTools NUnrooted Log2Unrooted Log2UnrootedMult NUnrootedMult
#' @export
ICSteps <- function (char, ambiguousToken = '?', tolerance = 0.05,
                     maxIter = 40000L, warn = TRUE) {
  #char <- matrix(2L ^ char[char != ambiguousToken], ncol = 1L)
  char <- char[char != ambiguousToken]
  names(char) <- paste0('t', seq_along(char))
  
  split <- sort(as.integer(table(char)))
  singletons <- split == 1L
  nSingletons <- sum(singletons)
  nInformative <- sum(split) - nSingletons
  minSteps <- length(split) - 1L
  if (minSteps == 0L) return (c('0' = 1L))
  
  split <- split[!singletons]
  if (length(split) < 2L) return (1L)
  
  
  nNoExtraSteps <- NUnrootedMult(split)
  log2ExtraSteps <- LnUnrootedMult(split)
  nOneExtraStep <- WithOneExtraStep(split)
  
  p1 <- nOneExtraStep / NUnrooted(sum(split))
  ic1 <- -log2(nOneExtraStep / NUnrooted(sum(split)))
  icTol <- tolerance - ic1
  pTol <- (2 ^ icTol) - p1
  iter <- p1 * (1 - p1) / (pTol ^ 2)
  
  nIter <- min(maxIter, round(iter))
  if (nIter == maxIter && warn) {
    warning ("Will truncate number of iterations at maxIter = ", maxIter)
  }
  n01ExtraSteps <- nOneExtraStep + nNoExtraSteps
  analyticIC <- Log2Unrooted(sum(split)) -  setNames(c(
    Log2UnrootedMult(split), log2(n01ExtraSteps)),
    minSteps + 0:1)
  analyticP <- 2 ^ -analyticIC[2]
  
  if (warn) {
    message('  Token count ', split, " = ",
            signif(analyticIc0, ceiling(log10(maxIter))),
            ' bits @ 0 extra steps. \n  Simulating ', nIter, 
            ' trees to estimate cost of further steps.')
    # message(c(round(analyticIc0, 3), 'bits @ 0 extra steps;', round(analyticIc1, 3),
    #    '@ 1; attempting', nIter, 'iterations.\n'))
  }
  
  morphyObj <- SingleCharMorphy(rep(seq_along(split) - 1L, split))
  on.exit(morphyObj <- UnloadMorphy(morphyObj))
  steps <- vapply(rep(nInformative, maxIter), RandomTreeScore,
                  integer(1), morphyObj)
  
  tabSteps <- table(steps[steps > (minSteps - nSingletons + 1)]) # Quicker than table(steps)[-1]
  
  approxP <- tabSteps / sum(tabSteps) * (1 - analyticP)
  approxSE <- sqrt(approxP * (1 - approxP) / nIter)
  cumP <- cumsum(c(analyticP, approxP))[-1]
                 
  approxIC <- -log2(cumP)
  icLB <- -log2(cumP - approxSE)
  icError <- icLB - approxIC
  if (warn || max(icError) > tolerance) {
    message("  Approx. std. error < ", signif(max(icError) * 1.01, 2))
  }
  ret <- c(analyticIC, approxIC)
  ret[length(ret)] <- 0 # Floating point error inevitable
  
  # Return:
  ret
}

#' @describeIn ICPerStep Memoized calculating function
#'
#' @param a \code{min(splits)}.
#' @param b \code{max(splits)}.
#' @param m \code{maxIter}.
#' @template warnParam
#'
#' @importFrom R.cache addMemoization
#' @keywords internal
#' @export
ICS <- addMemoization(function(a, b, m, warn = TRUE) {
  ICSteps(c(rep(1, a), rep(2, b)), maxIter = m, warn = warn)
})

#' Information content per step
#' @template splitsParam
#' @param maxIter number of iterations to use when estimating concavity constant
#' @template warnParam
#' @export
ICPerStep <- function(splits, maxIter, warn = TRUE) {
  ICS(min(splits), max(splits), maxIter, warn)
}

#' Number of trees with _m_ additional steps
#' 
#' Calculate the number of trees with _m_ extra steps under Fitch parsimony
#' where _a_ leaves are labelled with one state, and _b_ leaves labelled with
#' a second state.
#' 
#' Implementation of theorem 1 from Carter _et al._ (1990)
#' 
#' @param m Number of steps
#' @param a,b Number of leaves labelled `0` and `1`.
#' 
#' @references 
#' \insertRef{Carter1990}{TreeTools}
#' See also:
#' 
#' \insertRef{Steel1993}{TreeSearch}
#' 
#' \insertRef{Steel1995}{TreeSearch}
#' @importFrom TreeTools DoubleFactorial
#' @export
Carter1 <- function (m, a, b) {
  n <- a + b
  twoN <- n + n
  twoM <- m + m
  N <- function (n, m) {
    if (n < m) 0 else {
      nMinusM <- n - m
      factorial(n + nMinusM - 1L) / prod(
        factorial(nMinusM),
        factorial(m - 1L),
        2 ^ (nMinusM))
    }
  }
  prod(
    factorial(m - 1L),
    (twoN - twoM - m),
    DoubleFactorial(twoN - 5L),
    N(a, m),
    N(b, m)
  ) / DoubleFactorial(twoN - twoM - 1L)
}

#' @rdname Carter1
#' @export
Log2Carter1 <- function (m, a, b) {
  n <- a + b
  twoN <- n + n
  twoM <- m + m
  Log2N <- function (n, m) {
    if (n < m) 0 else {
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


#' @param n Number of leaves in tree
#' @param x Vector specifying number of leaves with state 1, 2, ...
#' @param b State reconstructed at base of tree
#' States are converted to binary, so 
#' `1` denotes 'state 1',
#' `2` denotes 'state 2', and
#' `3` denotes 'state 1 OR state 2'.
#' 
#' @references
#' \insertRef{Maddison1991}{TreeSearch}
#' @examples
#' @template MRS
#' @importFrom TreeTools NRooted
#' @export
MaddisonSlatkin <- function (a, b) {
  n <- sum(a, b)
  i <- min(a, b)
  
  # D: The probability that, in a randomly selected tree of _n_ taxa, _i_ of
  # which have state 1, the smaller subclade of _m_ taxa will receive _j_ taxa
  # with state 1
  .D <- function (j, i, m, n) {
    choose(i, j) * choose(n - i, m - j) / choose(n, m)
  }
  # R: The probability that in a randomly selected tree of _n_ leaves, the
  # smaller of the two basal subclades will have _m_ leaves.
  # My function gives half the value given in Maddison & Slatkin 1991 [!]
  .R <- function (m, n) {
    prod(ifelse(n == m + m, 1L, 2L) / 2L, 
         choose(n, m),
         NRooted(m),
         NRooted(n - m)
    ) / NRooted(n)
  } #TODO UNTESTED: Is it quicker to use LogNRooted?
  
  # B{b}: The probability that state _b_ is reconstructed at the base of a clade
  # with _n_ leaves and _i_ leaves with state 1.
  # My function runs from j = 0 to i, not j = 1 to i as printed in M&S91 [!]
  .B0 <- function (n, i) {
    if (n == 2L && i == 1L) {
      0
    } else if (n == i) {
      0
    } else if (i == 0L) {
      1
    } else {
      sum(vapply(seq_len(n / 2L), function (m) {
        .R(m, n) * sum(vapply(0:(min(i, m)), function (j) {
          .D(j, i, m, n) * sum(
            .B0(m, j) * .B0(n - m, i - j),
            .B01(m, j) * .B0(n - m, i - j),
            .B0(m, j) * .B01(n - m, i - j)
          )
        }, double(1)))
      }, double(1)))
    }
  }
  .B1 <- function (n, i) {
    if (n == 2L && i == 1L) {
      0
    } else if (n == i) {
      1
    } else if (i == 0L) {
      0
    } else {
      sum(vapply(seq_len(n / 2L), function (m) {
        .R(m, n) * sum(vapply(0:min(i, m), function (j) {
          .D(j, i, m, n) * sum(
            .B1(m, j) * .B1(n - m, i - j),
            .B01(m, j) * .B1(n - m, i - j),
            .B1(m, j) * .B01(n - m, i - j)
          )
        }, double(1)))
      }, double(1)))
    }
  }
  .B01 <- function (n, i) {
    if (n == 2L && i == 1L) {
      1
    } else if (n == i) {
      0
    } else if (i == 0L) {
      0
    } else {
      sum(vapply(seq_len(n / 2L), function (m) {
        .R(m, n) * sum(vapply(0:min(i, m), function (j) {
          .D(j, i, m, n) * sum(
            .B01(m, j) * .B01(n - m, i - j),
            .B0(m, j) * .B1(n - m, i - j),
            .B1(m, j) * .B0(n - m, i - j)
          )
        }, double(1)))
      }, double(1)))
    }
  }
  .PRow <- function (P1, B1, P2, B2, i, j, m, n, s, rMax = s) {
    sum(vapply(0:rMax, function (r) prod(
      P1(r, m, j),
      B1(m, j),
      P2(s - r, n - m, i - j),
      B2(n - m, i - j)), double(1)))
  }
  
  # P: The probability of having _s_ steps in a clade of _n_ leaves given that
  # _i_ of the leaves have state 1.
  # P{b}: P, given state _b_ is reconstructed at the base of the clade.
  .P0 <- function (s, n, i) {
    if (n == 2L && i == 1L) {
      if (s == 1L) 1 else 0
    } else if (n == i) {
      0 # Undefined?
    } else if (i == 0L) {
      if (s == 0L) 1 else 0
    } else {
      b0 <- .B0(n, i)
      if (b0 == 0) 0 else {
        sum(vapply(seq_len(n / 2), function (m) {
          .R(m, n) / b0 * sum(vapply(seq_len(min(i, m)), function (j) {
            .D(j, i, m, n) * sum(
              .PRow(.P0, .B0, .P0, .B0, i, j, m, n, s),
              .PRow(.P0, .B0, .P01, .B01, i, j, m, n, s),
              .PRow(.P01, .B01, .P0, .B0, i, j, m, n, s)
              )
          }, double(1)))
        }, double(1)))
      }
    }
  }
  .P1 <- function (s, n, i) {
    if (n == 2L && i == 1L) {
      0
    } else if (n == i) {
      1
    } else if (i == 0L) {
      0
    } else {
      b1 <- .B1(n, i)
      if (b1 == 0) 0 else {
        sum(vapply(seq_len(n / 2), function (m) {
          .R(m, n) / b1 * sum(vapply(seq_len(min(i, m)), function (j) {
            .D(j, i, m, n) * sum(
              .PRow(.P1, .B1, .P1, .B1, i, j, m, n, s),
              .PRow(.P1, .B1, .P01, .B01, i, j, m, n, s),
              .PRow(.P01, .B01, .P1, .B1, i, j, m, n, s)
            )
          }, double(1)))
        }, double(1)))
      }
    }
  }
  .P01 <- function (s, n, i) {
    if (n == 2L && i == 1L) {
      1
    } else if (n == i) {
      0
    } else if (i == 0L) {
      0
    } else {
      b01 <- .B01(n, i)
      if (b01 == 0) 0 else {
        sum(vapply(seq_len(n / 2), function (m) {
          .R(m, n) / b01 * sum(vapply(seq_len(min(i, m)), function (j) {
            .D(j, i, m, n) * sum(
              .PRow(.P1, .B1, .P1, .B1, i, j, m, n, s),
              .PRow(.P1, .B1, .P01, .B01, i, j, m, n, s, s - 1),
              .PRow(.P01, .B01, .P1, .B1, i, j, m, n, s, s - 1)
            )
          }, double(1)))
        }, double(1)))
      }
    }
  }
  .P0(0, 3, 2)
  .P1(0, 3, 2)
  .P01(0, 3, 2)
  
  .P0(0, 3, 1)
  .P1(0, 3, 1)
  .P01(0, 3, 1)
  
  .P0(1, 3, 1)
  .P1(1, 3, 1)
  .P01(1, 3, 1)
  
  P1(1, 3, 2)
  .P01(1, 3, 2)
  
  # Return:
  sum(
    .P0(s, n, i) * .B0(n, i),
    .P1(s, n, i) * .B1(n, i),
    .P01(s, n, i) * .B01(n, i))
}

# TODO: Remove the below, or if necessary replace with results of Carter et al. 1990; Steel 1993

#' Number of trees with one extra step
#' @param \dots Vector or series of integers specifying the number of leaves
#' bearing each distinct non-ambiguous token.
#' @importFrom TreeTools NRooted NUnrooted
#' @examples
#' WithOneExtraStep(1, 2, 3)
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
        attachTwoClades,
        sum(
        # TODO would be quicker to calculate just first half; special case: omitted.tips %% 2
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
  }
}

#' @importFrom TreeTools NUnrooted
#' @export
ExtraSteps <- function (a, b) {
  NRooted(a) * .ExtraSteps(a, b - 1L, 1L, 1L, 0L, 1L)
}

.ExtraSteps <- function (on1, off2, on2 = 0L,
                         clades2 = 0L, grades2 = 0L,
                         steps = 0L) {
  if (off2 == 0L) {
    return (c(1L, 0L, 0L, 0L, 0L))
  }
  counts <- integer(5L)
  # Don't change shape
  if (on2) {
    counts <- counts + (
                  (on2 + on2 + grades2 - clades2) *
                    .ExtraSteps(on1, off2 - 1L, on2 + 1L, clades2, grades2,
                                steps)
                  )
  }
  
  if (steps < 4L) {
    newCladePoints <- on1 + on1 - 3L - clades2 - grades2
    if (newCladePoints > 0L) {
      # Add a new clade
      counts[-1] <- counts[-1] + (
        newCladePoints *
        .ExtraSteps(on1, off2 - 1L, on2 + 1L, clades2 + 1L, grades2,
                    steps + 1L)[-5]
      )
    }
    # Convert clade to grade
    if (clades2 != 0) {
      counts[-1] <- counts[-1] + (
      (clades2 + clades2) *
      .ExtraSteps(on1, off2 - 1L, on2 + 1L, clades2 - 1L, grades2 + 1L,
                  steps + 1L)[-5]
      )
    }
  }
  
  # Return:
  counts
}

#' Logistic Points
#' Extract points from a fitted model
#'
#' @param x an integer vector giving x co-ordinates.
#' @param fittedModel a fitted model, perhaps generated using 
#' `nls(cumP ~ SSlogis(nSteps, Asym, xmid, scal))`.
#'
#' @return values of y co-ordinates corresponding to the x co-ordinates provided
#' @author Martin R. Smith
#' @export
LogisticPoints <- function (x, fittedModel) {
  coefL <- summary(fittedModel)$coef[, 1]
  # Return:
  coefL[1] / (1 + exp((coefL[2] - x) / coefL[3]))
}

#' Evaluate tree
#' @template treeParam
#' @template datasetParam
#' @template warnParam
#' @importFrom stats nls
#' @export
Evaluate <- function (tree, dataset, warn=TRUE) {
  totalSteps <- TreeLength(tree, dataset)
  chars <- matrix(unlist(dataset), attr(dataset, 'nr'))
  ambiguousToken <- which(attr(dataset, 'allLevels') == "?")
  asSplits <- apply(chars, 1, function (x) {
    ret <- table(x)
    ret[names(ret) != ambiguousToken] 
  })
  if (class(asSplits) == 'matrix') asSplits <- lapply(seq_len(ncol(asSplits)), function(i) asSplits[, i])
  ic.max <- round(vapply(asSplits,
                         function (split) -log(NUnrootedMult(split)/
                                                 NUnrooted(sum(split)))/log(2),
                         double(1)), 12)
  infoLosses <- apply(chars, 1, ICSteps, ambiguousToken=ambiguousToken, maxIter=1000, warn=warn)
  infoAmounts <- lapply(infoLosses, function(p) {
    #message(length(p))
    cumP <- cumsum(p)
    nSteps <- as.integer(names(p))
    infer <- min(nSteps):max(nSteps)
    infer <- infer[!(infer %in% nSteps)]
    calculatedP <- double(max(nSteps))
    calculatedP[nSteps] <- cumP
    if (length(infer)) {
      fitL <- nls(cumP ~ SSlogis(nSteps, Asym, xmid, scal))
      calculatedP[infer] <- LogisticPoints(infer, fitL)
    }
    calcIC <- -log(calculatedP) / log(2)
    calcIC
  })
  
  info.used <- double(length(totalSteps))
  for (i in seq_along(totalSteps)) {
    if (totalSteps[i] > 0) info.used[i] <- infoAmounts[[i]][totalSteps[i]]
  }
  info.lost <- round(ic.max - info.used, 13)
  index <- attr(dataset, 'index')
  total.info <- sum(ic.max[index])
  info.misleading <- sum(info.lost[index])
  proportion.lost <- info.misleading / total.info
  info.needed <- -log(1 / NUnrooted(length(dataset))) / log(2)
  info.overkill <- total.info / info.needed
  info.retained <- sum(info.used[index])
  signal.noise <- info.retained / info.misleading
  message(total.info, ' bits, of which ', round(info.retained, 2), ' kept, ',
          round(total.info - info.retained, 2), ' lost,',
          round(info.needed, 2), ' needed.  SNR = ', signal.noise, "\n")
  # Return:
  c(signal.noise, info.retained/info.needed)
}

#' Amount of information in each character
#'
#' As presently implemented, this function requires that there be no ambiguous
#' tokens and two applicable tokens, '1' and '2'.
#'
#' @param tokenTable A matrix, where each row corresponds to a character, each
#' column to a leaf, and each entry to the value (1 or 2) of the character at 
#' that leaf.
#' @param precision Integer specifyign how many random trees to generate when 
#' calculating profile curves.
#' @template warnParam
#'
#' @return information content of each extra step, in bits
#' 
#' @author Martin R. Smith
#'
#' @export
InfoAmounts <- function (tokenTable, precision = 1e+06, warn = TRUE) {
  # The below is simplified from info_extra_step.r::evaluate
  if (length(unique(as.integer(tokenTable))) > 2L) {
    stop ("Cannot calculate information amouts for",
          "characters unless only tokens are 1 and 2. See ?InfoAmounts().")
  }
  splits <- apply(tokenTable, 1, table)
  infoLosses <- apply(splits, 2, ICPerStep, maxIter = precision, warn = warn)
  
  blankReturn <- double(max(colSums(splits)))
  ret <- vapply(infoLosses, function(p) {
    calcIC <- blankReturn
    cumP <- cumsum(p)
    nSteps <- as.integer(names(p))
    infer <- min(nSteps):max(nSteps)
    infer <- infer[!(infer %in% nSteps)]
    calculatedP <- double(max(nSteps))
    calculatedP[nSteps] <- cumP
    if (length(infer)) {
      if (cumP[infer[1] - 1L] > 0.999998) {
        # We can cope with rounding at the sixth decimal place, I think
        calculatedP[infer] <- 1
      } else {
        fitL <- nls(cumP ~ SSlogis(nSteps, Asym, xmid, scal))
        calculatedP[infer] <- LogisticPoints(infer, fitL)
        ##if (warn) warning('Concavity function generated by approximation')
      }
    }
    calcIC[seq_along(calculatedP)] <- -log(calculatedP) / log(2)
    calcIC
  }, blankReturn)
  ret[rowSums(ret) > 0, ]
}
