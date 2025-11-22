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
                  - (LogCumSumExp(logProfile) / log(2)),
                  seq_len(split[2]) + sum(singletons))
  ret[ret < sqrt(.Machine[["double.eps"]])] <- 0 # Floating point error inevitable
  
  # Return:
  ret
}

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
