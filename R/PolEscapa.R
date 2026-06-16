#' Contribution of character to leaf instability
#' 
#' Would tree lengths change if a character was coded as ambiguous for each 
#' leaf \insertCite{Pol2009}{TreeSearch}?
#' 
#' High values for a leaf indicate that its coding contributes to instability
#' ("wildcard" or "roguish" behaviour; see
#' [\pkg{Rogue}](https://ms609.github.io/Rogue/) for further details).
#' The coding is in tension with other data, which may indicate that the
#' assumptions of homology that underlie the character's construction and
#' scoring require careful scrutiny &ndash; or that the taxon in question has
#' been subject to convergent evolution.
#' 
#' When inapplicable tokens are present in a character, the applicability of
#' each coding is maintained: i.e. a leaf coded with an applicable token is
#' never allowed to take an inapplicable value; and an inapplicable token
#' remains inapplicable.
#' 
#' @param trees List of trees of class `phylo`, or `multiPhylo` object.
#' @param char `phyDat` object containing a single character.
#' @inheritParams MaximizeParsimony
#' 
#' @return `LengthAdded()` returns a named numeric vector listing the mean
#' absolute change to tree length resulting if the character were coded
#' ambiguous for each leaf in turn, under the specified concavity constant.
#' 
#' @references \insertAllCited{}
#' @template MRS
#' @examples
#' trees <- inapplicable.trees[["Vinther2008"]]
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' char <- dataset[, 11]
#' added <- LengthAdded(trees, char)
#' 
#' PlotCharacter(
#'   tree = trees[[1]], 
#'   dataset = char,
#'   tip.color = 1 + added[trees[[1]]$tip.label] # Colour by added steps
#' ) -> XX # Suppress return value; display plot only
#' 
#' @family tree scoring
#' @export
LengthAdded <- function(trees, char, concavity = Inf) {
  if(!inherits(char, "phyDat")) {
    stop("`char` must be a character of class `phyDat`.")
  }
  if (attr(char, "nr") > 1L) {
    stop("`char` must comprise a single character; try char[, 1]")
  }
  cont <- attr(char, "contrast")
  zeroRows <- which(rowSums(cont) == 0)
  usedTokens <- unique(unlist(char, use.names = FALSE))
  if (any(zeroRows %in% usedTokens)) {
    stop("`char` contrast matrix lacks levels for token(s) ",
         paste(zeroRows[zeroRows %in% usedTokens], collapse = ", "))
  }
  if (inherits(trees, "phylo")) {
    trees <- c(trees)
  }
  
  trees <- RootTree(trees, 1) # Avoid warnings in TreeLength()
  start <- TreeLength(trees, char, concavity)
  contApp <- cont[, setdiff(colnames(cont), "-"), drop = FALSE]
  
  if (is.finite(concavity)) {
    # minLength attribute must be fixed.
    # Otherwise setting the only instance of a `1` to `?` will change the
    # calculated minimum length, potentially resulting in negative scores.
    char <- PrepareDataIW(char)
  } else if (.UseProfile(concavity)) {
    char <- PrepareDataProfile(char)
  }
  
  # Define ambiguous state, depending on applicability.
  # Take the first matching row when multiple rows are fully ambiguous, to
  # avoid silently assigning a vector to `charQm[[leaf]]` (analogous to the
  # T-302 fix for `qmApp`).
  qm <- which(rowSums(cont) == dim(cont)[2])
  if (length(qm) > 0L) {
    qm <- qm[[1L]]
  }
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
  } else {
    qmApp <- qmApp[[1L]]
  }
  # If no fully ambiguous (`?`) token row exists, `qm` is empty; a leaf whose
  # starting token is inapplicable would then be assigned `integer(0)` at
  # `charQm[[leaf]] <- qm`, silently corrupting the phyDat (dropping an
  # element).  Append an all-ones row (every state, applicable + inapplicable)
  # and point `qm` at it.  Done after the `qmApp` fallback so the row indices
  # of `cont`/`contApp`/`app`/`inapp` computed above remain consistent; this
  # only adds an extra row to the contrast that `qm` references.
  if (length(qm) == 0L) {
    newContrast <- rbind(attr(char, "contrast"), rep(1, dim(cont)[2]))
    attr(char, "contrast") <- newContrast
    qm <- nrow(newContrast)
  }
  
  QMScore <- function(leaf) {
    startToken <- char[[leaf]]
    if (!app[startToken]) {
      stopifnot(inapp[startToken])
      start
    } else {
      charQm <- char
      charQm[[leaf]] <- if (inapp[startToken]) qm else qmApp
      TreeLength(trees, charQm, concavity)
    }
  }
  
  deltas <- start - .vapply(seq_along(char), QMScore, start)
  delta <- setNames(colSums(deltas), names(char))
  
  # Return:
  delta / length(trees)
}

.vapply <- function(...) {
  ret <- vapply(...)
  if (is.null(dim(ret))) {
    ret <- matrix(ret, 1)
  }
  ret
}

#' @rdname LengthAdded
#' @export
PolEscapa <- LengthAdded
