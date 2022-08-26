#' Contribution of character to leaf instability
#' 
#' Would tree length change if a character was coded as ambiguous for each 
#' leaf \insertRef{Pol2009}{TreeSearch}?
#' 
#' When inapplicable tokens are present in a character, the applicability of
#' each coding is maintained: i.e. a leaf coded with an applicable token is
#' never allowed to take an inapplicable value; and an inapplicable token
#' remains inapplicable.
#' 
#' @param trees List of trees of class `phylo`, or `multiPhylo` object.
#' @param char `phyDat` object containing a single character.
#' @template concavityParam
#' 
#' @return `ExtraLength()` returns a named numeric vector listing the mean
#' absolute change to tree length resulting if the character were coded
#' ambiguous for each leaf in turn, under the specified concavity constant.
#' 
#' @references \insertAllCited{}
#' @template MRS
#' @examples
#' trees <- inapplicable.trees[["Vinther2008"]]
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' char <- dataset[, 11]
#' extra <- ExtraLength(trees, char)
#' 
#' PlotCharacter(
#'   tree = trees[[1]], 
#'   dataset = char,
#'   tip.color = 1 + extra[trees[[1]]$tip.label] # Colour by extra steps
#' ) -> XX # Suppress return value; display plot only
#' 
#' @export
ExtraLength <- function(trees, char, concavity = Inf, applicability = FALSE) {
  if(!inherits(char, "phyDat")) {
    stop("`char` must be a character of class `phyDat`.")
  }
  trees <- RootTree(trees, 1) # Avoid warnings in TreeLength()
  start <- TreeLength(trees, char, concavity)
  contApp <- cont[, setdiff(colnames(cont), "-")]
  
  # Define ambiguous state, depending on applicability
  qm <- which(rowSums(cont) == dim(cont)[2])
  if ("-" %in% colnames(cont)) {
    inapp <- as.logical(cont[, "-"])
    app <- as.logical(rowSums(contApp))
  } else {
    inapp <- logical(nrow(cont))
    app <- !inapp
  }
  inappLevel <- which.max(inapp)
  qmApp <- which(apply(contApp == 1, 1, all) & !inapp)
  if (length(qmApp) == 0) {
    attr(char, "contrast") <- rbind(cont, colnames(cont) != "-")
    qmApp <- 1 + nrow(cont)
  }
  
  QMScore <- function(leaf) {
    startToken <- char[[leaf]]
    if (!app[startToken]) {
      if (inapp[startToken]) {
        start
      } else {
        stop("Character has no valid levels for leaf ", leaf)
      }
    } else {
      charQm <- char
      charQm[[leaf]] <- if (inapp[startToken]) qm else qmApp
      TreeLength(trees, charQm, concavity)
    }
  }
  
  delta <- setNames(
    pmax(0, colSums(start - vapply(seq_along(char), QMScore, start))),
    # TODO pmax is a workaround for the Morphy bug that treats {+} as {?}.
    # TODO remove when that bug is fixed.
    names(char)
  )
  
  # Return:
  delta / length(trees)
}

#' @export
PolEscapa <- ExtraLength
