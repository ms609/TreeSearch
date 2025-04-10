#' Minimum and Maximum lengths possible for a character
#' 
#' The smallest and largest length that a phylogenetic character can attain on
#' any tree.
#' 
#' Ambiguous inapplicable states (e.g. `{0, -}`) are currently replaced with the
#' plain inapplicable token `-`, reflecting the current behaviour of Morphy.
#' 
#' @param x An object of class `phyDat`;
#' or a string to be coerced to a `phyDat` object via 
#' [`TreeTools::StringToPhyDat()`];
#' or an integer vector listing the tokens that may be present at each 
#' tip along a single character, with each token represented as a binary digit;
#' e.g. a value of 11 ( = 2^0 + 2^1 + 2^3) means that
#' the tip may have tokens 0, 1 or 3.
#' 
#' Inapplicable tokens should be denoted with the integer `0` (not 2^0).
#' 
#' @inheritParams CharacterLength
#' 
#' @return `MinimumLength()` returns a vector of integers specifying the 
#' minimum number of steps that each character must contain.
#'
#' @examples
#' data("inapplicable.datasets")
#' myPhyDat <- inapplicable.phyData[[4]]
#' 
#' # load your own data with
#' # my.PhyDat <- as.phyDat(read.nexus.data("filepath"))
#' # or Windows users can select a file interactively using:
#' # my.PhyDat <- as.phyDat(read.nexus.data(choose.files()))
#' 
#' class(myPhyDat) # phyDat object
#' 
#' # Minimum length of each character in turn
#' MinimumLength(myPhyDat)
#' 
#' # Collapse duplicate characters, per phyDat compression
#' MinimumLength(myPhyDat, compress = TRUE)
#' 
#' # Calculate length of a single character from its textual representation
#' MinimumLength("-{-1}{-2}{-3}2233")
#' MaximumLength("----0011")
#' @template MRS
#' @family tree scoring
#' @export
MinimumLength <- function (x, compress = FALSE) UseMethod("MinimumLength")

#' @rdname MinimumLength
#' @export
MinimumLength.phyDat <- function (x, compress = FALSE) {
  
  at <- attributes(x)
  levels <- at[["levels"]]
  nLevel <- length(levels)
  nChar <- at[["nr"]]
  nTip <- length(x)
  cont <- at[["contrast"]]
  if (is.null(colnames(cont))) {
    colnames(cont) <- as.character(levels)
  }
  
  inappLevel <- levels == "-"
  powersOf2 <- 2L ^ (seq_len(nLevel - sum(inappLevel)) - 1L)
  
  # Treat {-, 1} as {1}
  unlisted <- unlist(x, use.names = FALSE)
  tmp <- as.integer(cont[, colnames(cont) != "-", drop = FALSE] %*% powersOf2)
  ambigIsApp <- matrix(tmp[unlisted], nChar, nTip)
  
  if (any(inappLevel)) {
    # Treat {-, 1} as {-}
    tmp[cont[, "-"] == 1] <- 0
    ambigIsInapp <- matrix(tmp[unlisted], nChar, nTip)
    
    inappCount <- rowSums(matrix(unlisted %in% which(at[["allLevels"]] == "-"),
                                 nChar, nTip))
    binaryMatrix <- ambigIsApp
    binaryMatrix[inappCount > 1, ] <- ambigIsInapp[inappCount > 1, ]
  } else {
    binaryMatrix <- ambigIsApp
  }
  
  ret <- apply(binaryMatrix, 1, MinimumLength.numeric)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(x, "index")]
  }
}

#' @rdname MinimumLength
#' @export
MinimumLength.numeric <- function (x, compress = NA) {
  
  uniqueStates <- unique(x[x > 0])
  if (length(uniqueStates) < 2) {
    return (0)
  }
  tokens <- vapply(uniqueStates, intToBits, raw(32)) != 00
  tokens <- tokens[rowSums(tokens) > 0, ]
  
  lastDim <- dim(tokens)
  tokensUsed <- 0
  
  repeat {
    tokens <- tokens[!duplicated(tokens), , drop = FALSE]
    unambiguous <- colSums(tokens) == 1
    tokenNecessary <- rowSums(tokens[, unambiguous, drop = FALSE]) > 0
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

#' @importFrom TreeTools NexusTokens StringToPhyDat
.MumLength.character <- function (x, Func) {
  nTip <- length(NexusTokens(x[[1]]))
  vapply(x, function (x) Func(StringToPhyDat(x, nTip)), 1, USE.NAMES = FALSE)
}

#' @rdname MinimumLength
#' @export
MinimumLength.character <- function (x, compress = TRUE) {
  .MumLength.character(x, MinimumLength.phyDat)
}

#' @rdname MinimumLength
#' @export
MaximumLength.character <- function (x, compress = TRUE) {
  .MumLength.character(x, MaximumLength.phyDat)
}


#' @rdname MinimumLength
MinimumSteps <- function(x) {
  .Deprecated("MinimumLength",
              msg = "Renamed to `MinimumLength()` and recoded to better support inapplicable tokens")
  MinimumLength(x, compress = TRUE)
}

#' @rdname MinimumLength
#' @return `MaximumLength()` returns a vector of integers specifying the 
#' maximum number of steps that each character can attain in a parsimonious
#' reconstruction on a tree.  Inapplicable tokens are not yet supported.
#' @export
MaximumLength <- function(x, compress = TRUE) {
  UseMethod("MaximumLength")
}

#' @rdname MinimumLength
#' @export
MaximumLength.numeric <- function(x, compress = NA) {
  
  counts <- tabulate(x[x > 0])
  nInapp <- sum(x == 0)
  regions <- max(1, nInapp - 1)
  steps <- 0
  if (sum(counts > 0) > 1) {
    nState <- floor(log2(length(counts))) + 1L
    nToken <- 2 ^ nState - 1
    counts <- c(counts, double(nToken - length(counts)))
    
    rawTokens <- vapply(seq_len(nToken), intToBits, raw(32))[seq_len(nState), ]
    tokens <- rawTokens != 00
    tokenSums <- colSums(tokens)
    active <- c(rep(TRUE, nToken - 1), FALSE)
    
    bitTokens <- as.raw(seq_len(nToken))
    nonIntersect <- outer(bitTokens, bitTokens, `&`) == 00
    unions <- matrix(tokenSums[as.integer(outer(bitTokens, bitTokens, `|`))],
                     nToken, nToken)
    .Merge <- function(a, b) sum(2 ^ (which(tokens[, a] | tokens[, b]) - 1))
    loopCount <- 0
    
    
    # Think of this algorithm as building up a tree by joining regions in a
    # manner guaranteed to add a Fitch step, wherever possible.
    # 
    # The best we can do is include a leaf with each token 0..n to make a
    # region whose ancestral state is {012....n}
    # 
    # If we have inapplicable tokens, we can then plant these regions in a sea
    # of inapplicable tokens just large enough to make each region denote a
    # separate origin of each token, if we have enough (a) regions and (b)
    # inapplicable tokens.
    
    # Start with the token denoting the most ambiguous states
    repeat {
      amb <- max(-Inf, tokenSums[counts > 0 & active])
      if (amb < 1) {
        break
      }
      loopCount <- loopCount + 1
      if (loopCount > 1e4) {
        stop("MaximumLength() failed.",                                         # nocov
             " Please report this bug to TreeSearch maintainer.")               # nocov
      }
      escape <- FALSE
      
      #  We should optimally pair ...+++ with +++... to yield ++++++
      #  before considering ++.... for ++.+++
      for (unionSize in nState:(amb + 1)) {
        for (i in which(tokenSums == amb & counts > 0 & active)) {
          options <- counts > 0 & nonIntersect[i, ] & active
          candidates <- options & unions[i, ] == unionSize
          if (any(options)) {
            if (any(candidates)) {
              chosen <- which(candidates)[which.max(counts[candidates])]
              counts[c(i, chosen)] <- counts[c(i, chosen)] - 1
              product <- .Merge(i, chosen)
              counts[[product]] <- counts[[product]] + 1
              steps <- steps + 1
              escape <- TRUE
              break
            }
          } else {
            # There will never be a match for this ambiguity; disable 
            active[[i]] <- FALSE
          }
        }
        if (escape) {
          break
        }
      }
    }
  }
  # Return:
  steps + max(0, min(sum(counts), regions) - 1)
}

#' @export
MaximumLength.phyDat <- function (x, compress = FALSE) {
  
  at <- attributes(x)
  levels <- at[["levels"]]
  nLevel <- length(levels)
  nChar <- at[["nr"]]
  nTip <- length(x)
  cont <- at[["contrast"]]
  if (is.null(colnames(cont))) {
    colnames(cont) <- as.character(levels)
  }
  
  inappLevel <- levels == "-"
  powersOf2 <- 2L ^ (seq_len(nLevel - sum(inappLevel)) - 1L)
  
  # Treat {-, 1} as {1}
  unlisted <- unlist(x, use.names = FALSE)
  tmp <- as.integer(cont[, colnames(cont) != "-", drop = FALSE] %*% powersOf2)
  ambigIsApp <- matrix(tmp[unlisted], nChar, nTip)
  
  if (any(inappLevel)) {
    # Treat {-, 1} as {-}
    tmp[cont[, "-"] == 1] <- 0
    ambigIsInapp <- matrix(tmp[unlisted], nChar, nTip)
    
    inappCount <- rowSums(matrix(unlisted %in% which(at[["allLevels"]] == "-"),
                                 nChar, nTip))
    binaryMatrix <- ambigIsApp
    binaryMatrix[inappCount > 1, ] <- ambigIsInapp[inappCount > 1, ]
  } else {
    binaryMatrix <- ambigIsApp
  }
  
  ret <- apply(binaryMatrix, 1, MaximumLength.numeric)
  
  # Return:
  if (compress) {
    ret
  } else {
    ret[attr(x, "index")]
  }
}
