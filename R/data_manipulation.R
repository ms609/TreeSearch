#' Prepare data for Profile Parsimony
#' 
#' Calculates profiles for each character in a dataset.  Will also simplify
#' characters, with a warning, where they are too complex for the present
#' implementation of profile parsimony: 
#' - inapplicable tokens will be replaced with the ambiguous token
#'    (i.e. `-` \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} `?`);
#' - Ambiguous tokens will be treated as fully ambiguous
#'   (i.e. `{02}` \ifelse{html}{\out{&rarr;}}{\eqn{\rightarrow}{-->}} `?`)
#' - Where more than two states are informative (i.e. unambiguously present in
#'   more than one taxon), states beyond the two most informative will be
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
#'     \code{c("info.amounts", "split.sizes")}, indicating attributes to sample
#'      when bootstrapping the dataset (e.g. in Ratchet searches).
#'
#' `PrepareDataIW` adds the attribute:
#' 
#'  - \code{min.length}: The minimum number of steps that must be present in each
#'    transformation series.
#' @examples 
#' data("congreveLamsdellMatrices")
#' dataset <- congreveLamsdellMatrices[[42]]
#' PrepareDataProfile(dataset)
#' @author Martin R. Smith; written with reference to 
#' `phangorn:::prepareDataFitch()`
#' @importFrom cli cli_alert cli_alert_warning
#' @family profile parsimony functions
#' @encoding UTF-8
#' @export
PrepareDataProfile <- function (dataset) {
  if ("info.amounts" %fin% names(attributes(dataset))) {
    # Already prepared
    return(dataset)
  }
  at <- attributes(dataset)
  nLevel <- length(at[["levels"]])
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  index <- at[["index"]]
  allLevels <- as.character(at[["allLevels"]])
  
  contSums <- rowSums(cont)
  qmLevel <- which(contSums == ncol(cont))
  
  if (length(qmLevel) == 0) {
    attr(dataset, "contrast") <- rbind(attr(dataset, "contrast"), 1)
    attr(dataset, "allLevels") <- c(attr(dataset, "allLevels"), "{?}")
    qmLevel <- length(allLevels) + 1L
  }
  
  ambigs <- which(contSums > 1L & contSums < ncol(cont))
  inappLevel <- which(colnames(cont) == "-")
  if (length(inappLevel) != 0L) {
    cli_alert("Inapplicable tokens treated as ambiguous for profile parsimony")
    inappLevel <- which(apply(unname(cont), 1, identical,
                              as.double(colnames(cont) == "-")))
    dataset[] <- lapply(dataset, function (i) {
      i[i %fin% inappLevel] <- qmLevel
      i
    })
  }
  
  if (length(ambigs) != 0L) {
    # Message unnecessary until multiple informative states are supported
    # message("Ambiguous tokens ", paste(at[["allLevels"]][ambigs], collapse = ", "),
    #         " converted to "?"")
    dataset[] <- lapply(dataset, function (i) {
        i[i %fin% ambigs] <- qmLevel
        i
      })
  }
  
  mataset <- matrix(unlist(dataset, recursive = FALSE, use.names = FALSE),
                    max(index))
  
  .RemoveExtraTokens <- function (char, ambiguousTokens) {
    unambig <- char[!char %fin% ambiguousTokens]
    if (length(unambig) == 0) {
      return(matrix(nrow = length(char), ncol = 0))
    }
    split <- table(unambig)
    ranking <- order(order(split, decreasing = TRUE))
    ignored <- ranking > 2L
    if (any(split[ignored] > 1L)) {
      warningMsg <- "Can handle max. 2 informative tokens. Dropping others."
      if (interactive()) {
        cli_alert_warning(warningMsg)                                           # nocov
      } else {
        warning(warningMsg)
      }
    }
    if (length(ambiguousTokens) == 0) {
      stop("No ambiguous token available for replacement")
    }
    tokens <- names(split)
    most <- tokens[which.min(ranking)]
    vapply(setdiff(names(split)[split > 1], most), function (kept) {
           simplified <- char
           simplified[!simplified %fin% c(most, kept)] <- ambiguousTokens[1]
           simplified
    }, char)
  }
  
  decomposed <- lapply(seq_along(mataset[, 1]), function (i) 
    .RemoveExtraTokens(mataset[i, ], ambiguousTokens = qmLevel))
  nChar <- vapply(decomposed, dim, c(0, 0))[2, ]
  if (sum(nChar) == 0) {
    cli_alert("No informative characters in `dataset`.")
    attr(dataset, "info.amounts") <- double(0)
    return(dataset[0])
  }
  newIndex <- seq_len(sum(nChar))
  oldIndex <- rep.int(seq_along(nChar), nChar)
  index <- unlist(lapply(index, function (i) {
    newIndex[oldIndex == i]
  }))
  
  mataset <- unname(do.call(cbind, decomposed))
  
  NON_AMBIG <- 1:2
  AMBIG <- max(NON_AMBIG) + 1L
  .Recompress <- function (char, ambiguousTokens) {
    tokens <- unique(char)
    nonAmbig <- setdiff(tokens, ambiguousTokens)
    stopifnot(length(nonAmbig) == 2L)
    #available <- setdiff(seq_along(c(nonAmbig, ambiguousTokens)), ambiguousTokens)
    
    cipher <- seq_len(max(tokens))
    cipher[nonAmbig] <- NON_AMBIG # available[seq_along(nonAmbig)]
    cipher[ambiguousTokens] <- AMBIG
    
    # Return:
    cipher[char]
  }
  if (length(mataset) == 0) {
    cli_alert("No informative characters in `dataset`.")
    attr(dataset, "info.amounts") <- double(0)
    return(dataset[0])
  }
  mataset <- apply(mataset, 2, .Recompress, qmLevel)
  dupCols <- duplicated(t(mataset))
  kept <- which(!dupCols)
  copies <- lapply(kept, function (i) {
    i + which(apply(mataset[, -seq_len(i), drop = FALSE], 2, identical, mataset[, i]))
  })
  firstOccurrence <- seq_len(dim(mataset)[2])
  for (i in seq_along(copies)) {
    firstOccurrence[copies[[i]]] <- kept[i]
  }
  
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
    StepInformation(mataset[, i], ambiguousTokens = AMBIG))
  
  
  maxSteps <- max(vapply(info,
                         function (i) max(as.integer(names(i))),
                         integer(1)))
  info <- vapply(info,
                 function (x) {
                    ret <- setNames(double(maxSteps), seq_len(maxSteps))
                    x <- x[setdiff(names(x), "0")]
                    if (length(x)) {
                      ret[names(x)] <- max(x) - x
                    }
                    ret
                  }, double(maxSteps))
  if (is.null(dim(info))) {
    dim(info) <- c(1L, length(info))
  }
  attr(dataset, "index") <- index
  weight <- as.integer(table(index))
  attr(dataset, "weight") <- weight
  attr(dataset, "nr") <- length(weight)
  attr(dataset, "info.amounts") <- info
  attr(dataset, "informative") <- colSums(info) > 0
  lvls <- c("0", "1")
  attr(dataset, "levels") <- lvls
  attr(dataset, "allLevels") <- c(lvls, "?")
  attr(dataset, "contrast") <- matrix(c(1,0,1,0,1,1), length(lvls) + 1L, length(lvls), 
                                      dimnames = list(NULL, lvls))
  attr(dataset, "nc") <- length(lvls)
  
  if (!any(attr(dataset, "bootstrap") == "info.amounts")) {
    attr(dataset, "bootstrap") <- c(attr(dataset, "bootstrap"), "info.amounts")
  }
  
  dataset
}


#' @describeIn PrepareDataProfile Prepare data for implied weighting
#' @export
PrepareDataIW <- function (dataset) {
  attr(dataset, "min.length") <- MinimumLength(dataset, compress = TRUE)
  
  # Return:
  dataset
}
