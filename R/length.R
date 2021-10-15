#' Minimum length
#' 
#' The smallest length that a character can obtain on any tree.
#' 
#' 
#' @param x An object of class `phyDat`,
#' or an integer vector listing the tokens that may be present at each 
#' tip along a single character, with each token represented as a binary digit;
#' e.g. a value of 11 ( = 2^0 + 2^1 + 2^3) means that
#' the tip may have tokens 0, 1 or 3.
#' 
#' Inapplicable tokens should be denoted with the integer `0` (not 2^0).
#' 
#' Tokens that are ambiguous for an inapplicable and an applicable
#' state are not presently supported; for an approximate value, denote such
#' ambiguity with the integer `0`.
#' @template compressParam
#' 
#' @return `MinimumLength()` returns a vector of integers specifying the 
#' minimum number of steps that each character must contain.
#'
#' @examples
#' data('inapplicable.datasets')
#' myPhyDat <- inapplicable.phyData[[4]] 
#' MinimumLength(myPhyDat)
#' MinimumLength(myPhyDat, compress = TRUE)
#' 
#' 
#' class(myPhyDat) # phyDat object
#' # load your own data with
#' # my.PhyDat <- as.phyDat(read.nexus.data('filepath'))
#' # or Windows users can select a file interactively using:
#' # my.PhyDat <- as.phyDat(read.nexus.data(choose.files()))
#' 
#' # Convert list of character codings to an array
#' myData <- vapply(myPhyDat, I, myPhyDat[[1]])
#' 
#' # Convert phyDat's representation of states to binary
#' myContrast <- attr(myPhyDat, 'contrast')
#' tokens <- colnames(myContrast)
#' binaryContrast <- integer(length(tokens))
#' tokenApplicable <- tokens != '-'
#' binaryContrast[tokenApplicable] <- 2 ^ (seq_len(sum(tokenApplicable)) - 1)
#' binaryValues <- apply(myContrast, 1, 
#'   function (row) sum(binaryContrast[as.logical(row)]))
#' myStates <- matrix(binaryValues[myData], nrow = nrow(myData),
#'                    ncol = ncol(myData), dimnames = dimnames(myData))
#'
#' # Finally, work out minimum steps 
#' apply(myStates, 1, MinimumLength)
#' @template MRS
#' @family tree scoring
#' @export
MinimumLength <- function (x, compress = FALSE) UseMethod('MinimumLength')

#' @rdname MinimumLength
#' @export
MinimumLength.phyDat <- function (x, compress = FALSE) {
  
  at <- attributes(x)
  nLevel <- length(at$level)
  nChar <- at$nr
  nTip <- length(x)
  cont <- at$contrast
  if (is.null(colnames(cont))) {
    colnames(cont) <- as.character(at$levels)
  }
  simpleCont <- ifelse(rowSums(cont) == 1,
                       apply(cont != 0, 1, function (x) colnames(cont)[x][1]),
                       '?')
  inappLevel <- at$levels == '-'
  
  if (any(inappLevel)) {
    # TODO this is a workaround until MinimumLength.numeric can handle {-, 1}
    cont[cont[, inappLevel] > 0, ] <- 0
    ambiguousToken <- at$allLevels == '?'
    cont[ambiguousToken, ] <- colSums(cont[!ambiguousToken, ]) > 0
  }
  
  powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- as.integer(cont %*% powersOf2)
  unlisted <- unlist(x, use.names = FALSE)
  binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
  
  ret <- apply(binaryMatrix, 1, MinimumLength)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(x, 'index')]
  }
}

#' @rdname MinimumLength
#' @export
MinimumLength.numeric <- function (x, compress = NA) {
  
  uniqueStates <- unique(x[x > 0])
  if (length(uniqueStates) < 2) return (0)
  tokens <- vapply(uniqueStates, intToBits, raw(32)) != 00
  tokens <- tokens[apply(tokens, 1, any), ]
  
  lastDim <- dim(tokens)
  tokensUsed <- 0
  
  repeat {
    tokens <- tokens[!duplicated(tokens), , drop = FALSE]
    unambiguous <- colSums(tokens) == 1
    tokenNecessary <- apply(tokens[, unambiguous, drop = FALSE], 1, any)
    statesRemaining <- !unambiguous
    statesRemaining[statesRemaining] <- colSums(
      tokens[tokenNecessary, statesRemaining, drop = FALSE]) == 0
    tokensUsed <- tokensUsed + sum(tokenNecessary)
    
    if (!any(statesRemaining)) {
      # Return:
      return (tokensUsed - 1L)
    }
    
    tokens <- tokens[!tokenNecessary, statesRemaining, drop = FALSE]
    if (identical(dim(tokens), lastDim)) {
      occurrences <- rowSums(tokens)
      unnecessary <- occurrences == 1
      if (any(unnecessary)) {
        tokens <- tokens[!unnecessary, , drop = FALSE]
      } else {
        squish <- which.max(occurrences)
        tokensUsed <- tokensUsed + 1L
        tokens <- tokens[, tokens[!squish], drop = FALSE]
      }
    }
    lastDim <- dim(tokens)
  }
}

#' @rdname MinimumLength
MinimumSteps <- function(x) {
  .Deprecated("MinimumLength",
              msg = 'Renamed to `MinimumLength()` and recoded to better support inapplicable tokens')
  MinimumLength(x, compress = TRUE)
}
