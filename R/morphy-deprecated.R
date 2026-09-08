# Deprecated MorphyLib-era names, retained as thin shims over the native
# scoring layer (PrepareData / EdgeListScore / TreeScore / BootstrapTree /
# ReleaseData / is.ParsimonyData / RandomPostorderTree).  Scheduled for removal
# after 2.0.x.

#' Deprecated MorphyLib functions
#'
#' These MorphyLib-era functions are retained as deprecated aliases of their
#' native-scoring replacements and will be removed in a future release:
#' \itemize{
#'   \item `PhyDat2Morphy()` -> [`PrepareData()`]
#'   \item `SingleCharMorphy()` -> [`SingleCharData()`]
#'   \item `UnloadMorphy()` -> [`ReleaseData()`]
#'   \item `is.morphyPtr()` -> [`is.ParsimonyData()`]
#'   \item `MorphyLength()` -> [`EdgeListScore()`]
#'   \item `MorphyTreeLength()` -> [`TreeScore()`]
#'   \item `MorphyBootstrap()` -> [`BootstrapTree()`]
#'   \item `RandomMorphyTree()` -> [`RandomPostorderTree()`]
#' }
#'
#' @param phy,char,morphyObj,tree,parent,child,edgeList,nTip,gap,weight,inPostorder,nTaxa,\dots
#'  Passed to the replacement function.
#' @return The value of the corresponding replacement function.
#' @name TreeSearch-deprecated
#' @keywords internal
NULL

#' @rdname TreeSearch-deprecated
#' @export
PhyDat2Morphy <- function(phy, gap = "inapplicable",
                          weight = attr(phy, "weight")) {
  .Deprecated("PrepareData")
  if (!inherits(phy, "phyDat")) {
    stop("Invalid data type ", class(phy), "; should be phyDat.")
  }
  if (.GapHandler(gap) != "inapplicable") {
    stop("Only the `inapplicable` gap treatment is supported; recode your ",
         "data to treat gaps as missing or as an extra state.")
  }
  if (!is.null(weight)) {
    attr(phy, "weight") <- rep_len(weight, attr(phy, "nr"))
  }
  PrepareData(phy)
}

#' @rdname TreeSearch-deprecated
#' @export
SingleCharMorphy <- function(char, gap = "inapp") {
  .Deprecated("SingleCharData")
  if (.GapHandler(gap) != "inapplicable") {
    stop("Only the `inapplicable` gap treatment is supported.")
  }
  SingleCharData(char)
}

#' @rdname TreeSearch-deprecated
#' @export
UnloadMorphy <- function(morphyObj) {
  .Deprecated("ReleaseData")
  ReleaseData(morphyObj)
  invisible(0L)
}

#' @rdname TreeSearch-deprecated
#' @export
is.morphyPtr <- function(morphyObj) {
  .Deprecated("is.ParsimonyData")
  is.ParsimonyData(morphyObj)
}

#' @rdname TreeSearch-deprecated
#' @export
MorphyLength <- function(parent, child, morphyObj, inPostorder = FALSE,
                         nTaxa = NULL) {
  .Deprecated("EdgeListScore")
  EdgeListScore(parent, child, morphyObj, inPostorder = inPostorder)
}

#' @rdname TreeSearch-deprecated
#' @export
MorphyTreeLength <- function(tree, morphyObj) {
  .Deprecated("TreeScore")
  TreeScore(tree, morphyObj)
}

#' @rdname TreeSearch-deprecated
#' @export
MorphyBootstrap <- function(edgeList, morphyObj, ...) {
  .Deprecated("BootstrapTree")
  BootstrapTree(edgeList, morphyObj, ...)
}

#' @rdname TreeSearch-deprecated
#' @export
RandomMorphyTree <- function(nTip) {
  .Deprecated("RandomPostorderTree")
  RandomPostorderTree(nTip)
}
