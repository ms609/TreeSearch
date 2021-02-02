#' Prepare data for Profile Parsimony
#' 
#' @param dataset dataset of class \code{phyDat}
#' @param precision number of random trees to generate when calculating Profile curves. 
#'                  With 22 tokens (taxa):
#'                  - Increasing precision from 4e+05 to 4e+06 reduces error by a mean of 
#'                  0.005 bits for each step after the first (max = 0.11 bits, sd=0.017 bits)
#'                  - Increasing precision from 1e+06 to 4e+06 reduces error by a mean of 
#'                  0.0003 bits for each step after the first (max = 0.046 bits, sd=0.01 bits)
#'                  
#' @template warnParam
#'
#' @return An object of class `phyDat`, with additional attributes.
#' `PrepareDataProfile` adds the attributes:
#' 
#'   - \code{info.amounts}: details the information represented by each 
#'     character when subject to N additional steps.
#'   
#'   - \code{split.sizes}: The size of the splits implied by each character
#'   
#'   - \code{bootstrap}: The character vector 
#'     \code{c('info.amounts', 'split.sizes')}, indicating attributes to sample
#'      when bootstrapping the dataset (e.g. in Ratchet searches).
#'
#' `PrepareDataIW` adds the attribute:
#' 
#'  - \code{min.length}: The minimum number of steps that must be present in each
#'    transformation series.
#'
#' @author Martin R. Smith; written with reference to 
#' `phangorn:::prepareDataFitch()`
#' @export
PrepareDataProfile <- function (dataset, precision = 1e+06, warn = TRUE) {
  at <- attributes(dataset)
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  powers.of.2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  unlisted <- unlist(dataset, recursive=FALSE, use.names=FALSE)
  binaryMatrix <- tmp[unlisted]
  attr(binaryMatrix, 'dim') <- c(nChar, nTip)
  
  attr(dataset, 'info.amounts') <- InfoAmounts(binaryMatrix, precision, warn=warn)
  if (!any(attr(dataset, 'bootstrap') == 'info.amounts')) {
    attr(dataset, 'bootstrap') <- c(attr(dataset, 'bootstrap'), 'info.amounts')
  }
  
  ####  inappLevel <- which(at$levels == "-")
  ####  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  ####  
  ####  attr(dataset, 'split.sizes') <- apply(binaryMatrix, 1, function(x) {
  ####      vapply(applicableTokens, function (y) sum(x == y), integer(1))
  ####    })
  
  dataset
}


#' @describeIn PrepareDataProfile Prepare data for implied weighting
#' @export
PrepareDataIW <- function (dataset) {
  at <- attributes(dataset)
  nLevel <- length(at$level)
  nChar <- at$nr
  nTip <- length(dataset)
  cont <- at$contrast
  inappLevel <- at$levels == '-'
  
  if (any(inappLevel)) {
    # TODO this is a workaround until MinimumLength can handle {-, 1}
    cont[cont[, inappLevel] > 0, ] <- 0
    ambiguousToken <- at$allLevels == '?'
    cont[ambiguousToken, ] <- colSums(cont[!ambiguousToken, ]) > 0
  }
  
  # Perhaps replace with previous code:
  # inappLevel <- which(at$levels == "-")
  # cont[, inappLevel] <- 0
  
  
  powersOf2 <- 2L ^ c(0L, seq_len(nLevel - 1L))
  tmp <- as.integer(cont %*% powersOf2)
  unlisted <- unlist(dataset, use.names = FALSE)
  binaryMatrix <- matrix(tmp[unlisted], nChar, nTip, byrow = FALSE)
  
  attr(dataset, 'min.length') <- apply(binaryMatrix, 1, MinimumLength)
  
  # Return:
  dataset
}

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
#' 
#' @return A vector of integers specifying the minimum number of steps that 
#' each character must contain.
#'
#' @examples
#' data('inapplicable.datasets')
#' myPhyDat <- inapplicable.phyData[[4]] 
#' MinimumLength(myPhyDat)
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
#'   
#'
#' @author Martin R. Smith
#' @export
MinimumLength <- function (x) UseMethod('MinimumLength')

#' @rdname MinimumLength
#' @export
MinimumLength.phyDat <- function (x) {
  
  at <- attributes(x)
  nLevel <- length(at$level)
  nChar <- at$nr
  nTip <- length(x)
  cont <- at$contrast
  if (is.null(colnames(cont))) colnames(cont) <- as.character(at$levels)
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
  minLength <- apply(binaryMatrix, 1, MinimumLength)
  
}

#' @rdname MinimumLength
#' @export
MinimumLength.numeric <- function (x) {
  
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
    
    if (!any(statesRemaining)) return (tokensUsed - 1L)
    
    tokens <- tokens[!tokenNecessary, statesRemaining, drop = FALSE]
    if (identical(dim(tokens), lastDim)) {
      unnecessary <- rowSums(tokens) == 1
      if (any(unnecessary)) {
        tokens <- tokens[!unnecessary, , drop = FALSE]
      } else {
        stop("The token configuration [", paste(states, collapse=" "), 
             "] is not correctly handled by MinimumLength()\n",
             "Please report this bug at ",
             "https://github.com/ms609/TreeSearch/issues/new")
      }
    }
    lastDim <- dim(tokens)
  }
  
}

#' @rdname MinimumLength
MinimumSteps <- function(states) {
  .Deprecated(MinimumLength, msg='Renamed and recoded to better support
              inapplicable tokens')
  MinimumLength(states)
}
