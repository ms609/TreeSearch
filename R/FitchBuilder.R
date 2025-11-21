#' @rdname Carter1
#' @examples
#' # Number of trees with 2 steps for character 0011122
#' FitchBuilder(2, c(2, 3, 2)) * NUnrooted(7)
#' 
#' @importFrom TreeTools LnRooted
#' @export
FitchBuilder <- function(steps, tokenCount) {
  nTaxa <- sum(tokenCount)
  nLevels <- length(tokenCount)
  nStates <- 2 ^ nLevels - 1
  states <- integer(nStates)
  states[2 ^ (seq_len(nLevels) - 1)] <- tokenCount
  dp <- `mode<-`(.DownpassOutcome(nLevels), "integer")
  dpStep <- attr(dp, "step")
  diag(dp) <- NA_integer_
  dp[lower.tri(dp, diag = TRUE)] <- NA_integer_
  
  .BR <- function(x) .BinaryRepresentation(x, digits = nLevels)
  
  .Recurse <- function(rootState, tokensAvailable, stepsAvailable) {
    
    if (stepsAvailable < 1) return(rootState)
    
    apply(which(dp == rootState, arr.ind = TRUE), 1, function(childStates) {
      lState <- childStates[[1]]
      rState <- childStates[[2]]
      
      if (lState == 0 && stepsAvailable == 0) {
        # We've hit the end of the line
        return(rootState)
      }
      
      addsStep <- dpStep[rbind(childStates)]
      lTokens <- .BR(lState)
      rTokens <- .BR(rState)
      
      tokensNeeded <- lTokens + rTokens
      if (any(tokensNeeded > tokensAvailable)) return(NULL)
      
      lStepsNeeded <- sum(lTokens) - 1L
      rStepsNeeded <- sum(rTokens) - 1L
      if (addsStep) {
        stepsAvailable <- stepsAvailable - 1L
      }
      if (lStepsNeeded + rStepsNeeded > stepsAvailable) return(NULL)
      
      lapply(lStepsNeeded:(stepsAvailable - rStepsNeeded), function(lSteps) {
        rSteps <- stepsAvailable - lSteps
        list("L" = .Recurse(lState, tokensAvailable - rTokens, lSteps),
             "R" = .Recurse(rState, tokensAvailable - lTokens, rSteps))
      })
    })
  }
  
  lapply(seq_len(nStates), .Recurse,
         tokensAvailable = tokenCount,
         stepsAvailable = steps)
}

.AsNewick <- function(x) {
  if (!is.list(x)) {
    # atomic: just return as character
    return(as.character(x))
  }
  # recursive case
  children <- vapply(x, to_newick, FUN.VALUE = character(1))
  paste0("(", paste(children, collapse = ", "), ")")
}

.BinaryRepresentation <- function(x, digits = 32) {
  intToBits(x)[1:digits] != as.raw(0)
}

# TODO Delete: # For reference, here's the definition of:
# .DownpassOutcome <- function(nStates) {
#   # Assume we have 3 states.  Number them from zero.
#   # Encode each state as 2^(0, 1, 2)
#   # Each node must exhibit a state: cannot be all FALSE
#   stateSpace <- as.raw(seq_len(2 ^ nStates - 1))
#   intersect <- outer(stateSpace, stateSpace, `&`)
#   union <- outer(stateSpace, stateSpace, `|`)
#   ret <- intersect
#   newStep <- intersect == as.raw(0)
#   ret[newStep] <- union[newStep]
#   structure(ret, step = newStep)
# }

