#' Prepare data for Profile Parsimony
#' 
#' @param dataset dataset of class \code{phyDat}
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
#' @examples 
#' data('congreveLamsdellMatrices')
#' dataset <- congreveLamsdellMatrices[[42]]
#' PrepareDataProfile(dataset)
#' @author Martin R. Smith; written with reference to 
#' `phangorn:::prepareDataFitch()`
#' @export
PrepareDataProfile <- function (dataset) {
  at <- attributes(dataset)
  nLevel <- length(at$level)
  nChar <- at$nr
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  index <- at$index
  allLevels <- as.character(at$allLevels)
  
  contSums <- rowSums(cont)
  qmLevel <- which(contSums == ncol(cont))
  ambigs <- which(contSums > 1L & contSums < ncol(cont))
  inappLevel <- which(colnames(cont) == '-')
  if (length(inappLevel) != 0L) {
    message("Inapplicable tokens treated as ambiguous for profile parsimony")
    inappLevel <- which(apply(unname(cont), 1, identical,
                              as.double(colnames(cont) == '-')))
    dataset[] <- lapply(dataset, function (i) {
      i[i %in% inappLevel] <- qmLevel
      i
    })
  }
  
  if (length(ambigs) != 0L) {
    message("Ambiguous tokens ", paste(at$allLevels[ambigs], collapse = ', '),
            " converted to '?'")
    dataset[] <- lapply(dataset, function (i) {
        i[i %in% ambigs] <- qmLevel
        i
      })
  }
  # TODO could improve efficiency by re-compressing, as some characters may
  # now be identical that previously differed in their ambiguous tokens
  
  mataset <- matrix(unlist(dataset, recursive = FALSE, use.names = FALSE),
                    max(index))
  #TODO when require R4.1: replace with
  # info <- apply(mataset, 1, StepInformation, 
  #               ambiguousTokens = c(qmLevel, inappLevel),
  #               simplify = FALSE)
  info <- lapply(seq_along(mataset[, 1]), function (i) 
    StepInformation(mataset[i, ], ambiguousTokens = c(qmLevel, inappLevel)))
  
  
  maxSteps <- max(vapply(info, function (i) max(as.integer(names(i))), integer(1)))
  info <- vapply(info,
                 function (x) {
                    ret <- double(maxSteps)
                    ret[seq_along(x)] <- max(x) - x
                    ret
                  }, double(maxSteps))
  if (is.null(dim(info))) {
    dim(info) <- c(1L, length(info))
  }
  attr(dataset, 'info.amounts') <- info
  
  if (!any(attr(dataset, 'bootstrap') == 'info.amounts')) {
    attr(dataset, 'bootstrap') <- c(attr(dataset, 'bootstrap'), 'info.amounts')
  }
  
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
  
  # Return:
  apply(binaryMatrix, 1, MinimumLength)
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
  MinimumLength(x)
}
