#' Group present or contradicted score
#' 
#' Implements the Groups Present / Contradicted (\acronym{GC}) measure
#' \insertCite{Goloboff2003}{TreeSearch}.
#' 
#' The GC score ranges from &#x2212;1 to 1, and is intended as an alternative
#' to bootstrap or jackknife support values.
#' 
#' The GC score counts the number of trees in `forest` that include a given
#' split, and subtracts the number of times that the *most frequent*
#' contradictory split occurs.  This value is then divided by the number of
#' trees in `forest`.
#' 
#' A score of 1 denotes that every tree in a forest (typically of bootstrap
#' or jackknife replicates) contains the split in question.
#' A score of &#x2212;1 denotes that a specific contradictory split occurs in
#' every tree in `forest`.
#' A score of zero indicates no support: i.e. that the split exhibits no more
#' support than its most common contradictory split.
#' 
#' The most frequent contradictory split is used to discriminate between
#' a scenario where a given split enjoys much more support than any other
#' alternative (even if many alternatives exist, each with low support),
#' and a scenario where the chosen split is scarcely any better supported than
#' a competing alternative.  The split is considered better supported than
#' the latter, where the runner-up may become preferred with a modest change to
#' the underlying dataset.
#' 
#' @template sprint
#' @returns `PresCont()`  returns a character vector that labels the nodes
#' in `tree` in sequence, suitable for passing to `nodelabels()` or
#' `tree$node.label`.
#' If `plot = TRUE`, it also plots `tree`, with splits labelled by their
#' groups present / contradicted score.
#' 
#' @inheritParams JackLabels
#' @param forest a list of trees of class `phylo`, or a `multiPhylo` object;
#'  or a `Splits` object.
# Remove when  https://github.com/r-lib/roxygen2/issues/1718 is closed:
# redundant to @inheritParams JackLabels
#' @param adj,col,frame,pos,\dots Parameters to pass to `nodelabels()`.
#' 
#' @seealso 
#' \code{\link[TreeTools]{SplitFrequency}()} and [`MostContradictedFreq()`] will
#' compute the number of trees that contain the split, and the frequency of the
#' most contradicted split, respectively.
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
#' # Compute the measure for further analysis
#' gpc <- (SplitFrequency(tree, jackTrees) -
#'   MostContradictedFreq(tree, jackTrees)) / length(jackTrees)
#' gpc
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
#' `MostContradictedFreq()` counts the occurrences of the single split that most
#' frequently contradicts each split in `tree`.
#' 
#' \insertCite{Goloboff2003;textual}{TreeSearch} propose comparing the frequency
#' of a split in a resampled population with the frequency of the most common
#' contradictory split.  This measure contributes to the "groups present /
#' contradicted" score.
#' 
#' @template sprint
#' @inheritParams PresCont
#' @returns `MostContradictedFreq()` returns, for each split in `tree`,
#' the number of times that its most common contradictory split occurs in
#' `forest`.
#' @seealso `PresCont()` calculates the "groups present / contradicted" score.
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
