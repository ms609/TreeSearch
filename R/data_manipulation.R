#' Prepare data for Profile Parsimony
#' 
#' Calculates profiles for each character in a dataset.  Will also simplify
#' characters, with a warning, where they are too complex for the present
#' implementation of profile parsimony: 
#' - inapplicable tokens will be replaced with the ambiguous token
#'    (i.e. `-` &rarr; `?`);
#' - Ambiguous tokens will be treated as fully ambiguous
#'   (i.e. `{02}` &rarr; `?`)
#' - Where more than two states are informative (i.e. unambiguously present in
#'   more than two taxa), states beyond the two most informative will be
#'   ignored.
#TODO can do something more complex like first two to one TS, second two to another   
#' 
#' @param dataset dataset of class \code{phyDat}
#'
#' @return An object of class `phyDat`, with additional attributes.
#' `PrepareDataProfile` adds the attributes:
#' 
#'   - `info.amounts`: details the information represented by each 
#'     character when subject to N additional steps.
#'   
#'   - `informative`: logical specifying which characters contain any
#'     phylogenetic information.
#'   
#'   - `bootstrap`: The character vector 
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
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  index <- at$index
  allLevels <- as.character(at$allLevels)
  
  contSums <- rowSums(cont)
  qmLevel <- which(contSums == ncol(cont))
  
  if (length(qmLevel) == 0) {
    attr(dataset, "contrast") <- rbind(attr(dataset, "contrast"), 1)
    attr(dataset, "allLevels") <- c(attr(dataset, "allLevels"), '{?}')
    qmLevel <- length(allLevels) + 1L
  }
  
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
    # Message unnecessary until multiple informative states are supported
    # message("Ambiguous tokens ", paste(at$allLevels[ambigs], collapse = ', '),
    #         " converted to '?'")
    dataset[] <- lapply(dataset, function (i) {
        i[i %in% ambigs] <- qmLevel
        i
      })
  }
  
  mataset <- matrix(unlist(dataset, recursive = FALSE, use.names = FALSE),
                    max(index))
  
  .RemoveExtraTokens <- function (char, ambiguousTokens) {
    unambig <- char[!char %in% ambiguousTokens]
    if (length(unambig) == 0) {
      return(char)
    }
    split <- table(unambig)
    ranking <- order(order(split, decreasing = TRUE))
    ignored <- ranking > 2L
    if (any(split[ignored] > 1L)) {
      warning("Can handle max. 2 informative tokens. Dropping others.")
    }
    if (length(ambiguousTokens) == 0) {
      stop("No ambiguous token available for replacement")
    }
    tokens <- names(split)
    most <- tokens[which.min(ranking)]
    vapply(setdiff(names(split)[split > 1], most), function (kept) {
           simplified <- char
           simplified[!simplified %in% c(most, kept)] <- ambiguousTokens[1]
           simplified
    }, char)
  }
  
  decomposed <- lapply(seq_along(mataset[, 1]), function (i) 
    .RemoveExtraTokens(mataset[i, ], ambiguousTokens = qmLevel))
  nChar <- vapply(decomposed, dim, c(0, 0))[2, ]
  attr(dataset, 'nr') <- sum(nChar)
  attr(dataset, 'weight') <- rep.int(at$weight, nChar)
  newIndex <- seq_len(sum(nChar))
  oldIndex <- rep.int(seq_along(nChar), nChar)
  index <- unlist(lapply(index, function (i) {
    newIndex[oldIndex == i]
  }))
  
  mataset <- unname(do.call(cbind, decomposed))
  
  .Recompress <- function (char, ambiguousTokens) {
    tokens <- unique(char)
    nonAmbig <- setdiff(tokens, ambiguousTokens)
    stopifnot(length(nonAmbig) == 2L)
    #available <- setdiff(seq_along(c(nonAmbig, ambiguousTokens)), ambiguousTokens)
    
    cipher <- seq_len(max(tokens))
    cipher[nonAmbig] <- 1:2 # available[seq_along(nonAmbig)]
    cipher[ambiguousTokens] <- 3
    
    # Return:
    cipher[char]
  }
  mataset <- apply(mataset, 2, .Recompress, qmLevel)
  if (length(mataset) == 0) {
    stop("No informative characters in `dataset`.")
  }
  dupCols <- duplicated(t(mataset))
  copies <- lapply(which(!dupCols), function (i) {
    i + which(apply(mataset[, -seq_len(i), drop = FALSE], 2, identical, mataset[, i]))
  })
  firstOccurrence <- seq_len(dim(mataset)[2])
  # This could be faster:
  for (i in seq_along(copies)) {
    firstOccurrence[copies[[i]]] <- i
  }
  kept <- unique(firstOccurrence)
  cipher <- seq_len(max(kept))
  cipher[kept] <- order(kept)
  index <- cipher[firstOccurrence][index]
  
  mataset <- mataset[, !dupCols, drop = FALSE]
  dataset[] <- lapply(seq_len(length(dataset)), function (i) mataset[i, ])
  
  
  #TODO when require R4.1: replace with
  # info <- apply(mataset, 1, StepInformation, 
  #               ambiguousTokens = c(qmLevel, inappLevel),
  #               simplify = FALSE)
  info <- lapply(seq_along(mataset[1, ]), function (i) 
    StepInformation(mataset[, i], ambiguousTokens = qmLevel))
  
  
  maxSteps <- max(vapply(info, function (i) max(as.integer(names(i))), integer(1)))
  info <- vapply(info,
                 function (x) {
                    ret <- setNames(double(maxSteps), seq_len(maxSteps))
                    x <- x[setdiff(names(x), '0')]
                    if (length(x)) {
                      ret[names(x)] <- max(x) - x
                    }
                    ret
                  }, double(maxSteps))
  if (is.null(dim(info))) {
    dim(info) <- c(1L, length(info))
  }
  attr(dataset, 'index') <- index
  weight <- as.integer(table(index))
  attr(dataset, 'weight') <- weight
  attr(dataset, 'nr') <- length(weight)
  attr(dataset, 'info.amounts') <- info
  attr(dataset, 'informative') <- colSums(info) > 0
  lvls <- c('0', '1')
  attr(dataset, 'levels') <- lvls
  attr(dataset, 'allLevels') <- c(lvls, '?')
  attr(dataset, 'contrast') <- matrix(c(1,0,1,0,1,1), length(lvls) + 1L, length(lvls), 
                                      dimnames = list(NULL, lvls))
  attr(dataset, 'nc') <- length(lvls)
  
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
