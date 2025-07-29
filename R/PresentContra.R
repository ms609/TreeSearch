#' Group present or contradicted
#' 
#' Implements the Groups Present / Contradicted (GC) measure
#' \insertCite{Goloboff2003}{TreeSearch}.
#' 
#' @template sprint
#' 
#' @references \insertAllCited{}
#' @examples
#' library("TreeTools", quietly = TRUE) # for as.phylo
#' 
#' # jackTrees will usually be generated with Jackknife() or Resample(),
#' # but for simplicity:
#' jackTrees <- as.phylo(1:100, 8)
#' 
#' tree <- as.phylo(0, 8)
#' PresCont(tree, jackTrees)
#' 
#' tree$node.label <- PresCont(tree, jackTrees, plot = FALSE)
#' 
#' # Write the labelled tree to screen
#' ape::write.tree(tree)
#'
#' # Write labelled trees to a nexus file:
#' # write.nexus(tree, file = filename)
#' 
#' @inheritParams JackLabels
#' @inheritParams MostContradictedFreq
#' @template MRS
#' @family split support functions
#' @importFrom ape nodelabels
#' @export
PresCont <- function(tree, forest,
                     plot = TRUE,
                     add = FALSE,
                     adj = 0, col = NULL, frame = "none", pos = 2L,
                     ...) {
  gc <- (SplitFrequency(tree, forest) -
           MostContradictedFreq(tree, forest)) / length(forest)
  
  if (plot) {
    if (!add) plot(tree)
    if (is.null(col)) {
      col <- SupportColour(gc)
    }
    nodelabels(paste("\n\n", signif(gc, 2)),
               node = as.integer(names(gc)),
               adj = adj, col = col, pos = pos, frame = frame, ...)
    
    # Return:
    gc
  } else {
    ret <- character(tree[["Nnode"]])
    ret[as.integer(names(gc)) - NTip(tree)] <- gc
    
    # Return:
    ret
  }
}


#' Frequency of most common contradictory split
#' 
#' \insertCite{Goloboff2003;textual}{TreeSearch} propose comparing the frequency
#' of a split in a resampled population with the frequency of the most common
#' contradictory split.
#' @param forest a list of trees of class `phylo`, or a `multiPhylo` object;
#'  or a `Splits` object.
#' 
#' @returns `MostContradictedFreq()` returns, for each split in `tree`,
#' the number of times that its most common contradictory split occurs in
#' `forest`.
#' @importFrom TreeTools SplitConflicts
#' @export
MostContradictedFreq <- function(tree, forest) {
  referenceSplits <- as.Splits(tree)
  refLabels <- attr(referenceSplits, "tip.label")
  forest <- lapply(forest, KeepTip, refLabels)
  forestSplits <- as.Splits(forest, tipLabels = refLabels)
  
  setNames(
    vapply(seq_along(referenceSplits), function(i) {
      conflicts <- SplitConflicts(referenceSplits[[i]], forestSplits)
      badSplits <- Map(function(sp, keep) sp[[keep]], forestSplits, conflicts)
      strSplits <- unlist(lapply(badSplits, as.character), use.names = FALSE,
                          recursive = FALSE)
      max(tabulate(fmatch(strSplits, unique(strSplits))))
    }, integer(1)),
    names(referenceSplits))
}
