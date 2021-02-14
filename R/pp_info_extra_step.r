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
#' @param char The character in question.
#' @param ambiguousToken Which token, if any, corresponds to the ambiguous token
#' (?) (not yet fully implemented).
#' @param expectedMinima sample enough trees that the rarest step counts is
#' expected to be seen at least this many times.
#' @param maxIter Maximum iterations to conduct.
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
#' suppressWarnings(ICSteps(character))
#' 
#' # MS Temp Example testing
#' char <- character
#' ambiguousToken = '?'
#' expectedMinima = 25
#' maxIter = 100
#' warn = TRUE
#' 
#' @importFrom TreeTools NUnrooted Log2Unrooted Log2UnrootedMult NUnrootedMult
#' @export
ICSteps <- function (char, ambiguousToken = '?', expectedMinima = 25L,
                     maxIter = 10000L, warn = TRUE) {
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
  #nOneExtraStep <- WithOneExtraStep(split)
  
  nIter <- min(maxIter, round(expectedMinima * NUnrooted(nInformative) / nNoExtraSteps))
  if (nIter == maxIter && warn) warning ("Will truncate number of iterations at maxIter = ", maxIter)
  analyticIc0 <- Log2Unrooted(sum(split)) - Log2UnrootedMult(split)
  #analyticIc1<- -log(nOneExtraStep/NUnrooted(sum(split))) / log(2)
  #analyticIc1<- -log((nNoExtraSteps + nOneExtraStep)/NUnrooted(sum(split))) / log(2)
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
  
  tabSteps <- table(steps[steps > (minSteps - nSingletons + 0)]) # Quicker than table(steps)[-1]
  #tabSteps <- table(steps[steps > (minSteps + 1)])
  
  tabSteps <- c(analyticSteps, tabSteps * (nIter - sum(analyticSteps)) / sum(tabSteps))
  pSteps <- tabSteps / sum(tabSteps)
  
  # Return:
  pSteps
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


#' Number of trees with one extra step
#' @template splitsParam
#' @importFrom TreeTools NRooted
#' @export
WithOneExtraStep <- function (splits) {
  # Ignore singletons, which can be added at the end...
  splits.withSplitstables <- splits[splits > 1]
  if (length(splits.withSplitstables) < 2) return (0)
  
  # TODO test splits: 1 1 2 4, splits: 2 2 4
  sum(vapply(seq_along(splits), function (omit) {
    backboneSplits <- splits[-omit]
    omitted.tips <- splits[omit]
    if (omitted.tips < 2) return (0)
    backbone.tips <- sum(backboneSplits)
    backbones <- NUnrootedMult(backboneSplits)
    backbone.edges <- max(0, 2 * backbone.tips - 3)
    backbone.attachments <- backbone.edges * (backbone.edges - 1)
    prod(sum( # Branch unambiguously splits along first group
      vapply(1:(omitted.tips - 1), function (first.group) { # For each way of splitsting up the omitted tips, e.g. 1|16, 2|15, 3|14, etc
        choose(omitted.tips, first.group) * 
          NRooted(first.group) * NRooted(omitted.tips - first.group)
      }, double(1))
    ) / 2, backbone.attachments, backbones) + prod(
      # Second group added adjacent to first group, thus new edge could belong to the backbone or the omitted tip group
      sum(vapply(1:(omitted.tips - 1), function (first.group) { # For each way of splitsting up the omitted tips, e.g. 1|16, 2|15, 3|14, etc
        choose(omitted.tips, first.group) * 
          NRooted(first.group) * NRooted(omitted.tips - first.group) # backbone tips have already been splits - when we selected a branch
      }, double(1))) / 2,
      backbones,
      backbone.edges,
      2 # left or right of group addition location
      / 2 # Will be counted again when 'added group' becomes the 'backbone group'
    )
    
  }, double(1))
  )
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
  if (length(unique(as.integer(tokenTable))) > 2) {
    stop ("Cannot calculate information amouts for",
          "characters unless only tokens are 1 and 2. See ?InfoAmounts.")
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
