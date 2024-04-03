#' Jackknife resampling
#' 
#' Resample trees using Jackknife resampling, i.e. removing a subset of
#' characters.
#' 
#' The function assumes 
#' that `InitializeData()` will return a morphy object; if this doesn't hold 
#' for you, post a [GitHub issue](https://github.com/ms609/TreeSearch/issues/new/)
#' or e-mail the maintainer.
#' 
#' @inheritParams Ratchet
#' @template EdgeSwapperParam
#' @param resampleFreq Double between 0 and 1 stating proportion of characters 
#' to resample.
#' @param jackIter Integer specifying number of jackknife iterations to conduct.
#' @return `Jackknife()` returns a list of trees recovered after jackknife
#' iterations.
#' @template MRS
#' @importFrom TreeTools RenumberEdges RenumberTips
#' @seealso 
#' - [`JackLabels()`]: Label nodes of a tree with jackknife supports.
#' @family split support functions
#' @family custom search functions
#' @export
Jackknife <- function (tree, dataset, resampleFreq = 2/3,
                       InitializeData = PhyDat2Morphy,
                       CleanUpData    = UnloadMorphy,
                       TreeScorer     = MorphyLength,
                       EdgeSwapper    = TBRSwap,
                       jackIter   = 5000L,
                       searchIter = 4000L, searchHits = 42L,
                       verbosity = 1L, ...) {
  # initialize tree and data
  if (dim(tree[["edge"]])[1] != 2 * tree[["Nnode"]]) {
    stop("tree must be bifurcating; try rooting with ape::root")
  }
  tree <- RenumberTips(tree, names(dataset))
  edgeList <- tree[["edge"]]
  edgeList <- RenumberEdges(edgeList[, 1], edgeList[, 2])
  
  morphyObj <- InitializeData(dataset)
  on.exit(morphyObj <- CleanUpData(morphyObj))
  
  startWeights <- MorphyWeights(morphyObj)["exact", ]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep.int(eachChar, startWeights)
  charsToKeep <- ceiling(resampleFreq * length(deindexedChars))
  if (charsToKeep < 1L) {
    stop("resampleFreq of ", resampleFreq, " is too low; can't keep 0 of ",
         length(deindexedChars), " characters.")
  } else if (charsToKeep >= length(deindexedChars)) {
    stop("resampleFreq of ", resampleFreq, " is too high; can't keep all ",
         length(deindexedChars), " characters.")
  }
  if (verbosity > 10L) { #nocov start
    message(" * Beginning search:")
  } #nocov end
  
  # Conduct jackIter replicates:
  jackEdges <- vapply(seq_len(jackIter), function (x) {
    if (verbosity > 0L) { #nocov start
      message(" * Jackknife iteration ", x, "/", jackIter)
    } #nocov end
    resampling <- tabulate(sample(deindexedChars, charsToKeep, replace = FALSE),
                           nbins = length(startWeights))
    errors <- vapply(eachChar, function (i) 
      mpl_set_charac_weight(i, resampling[i], morphyObj), integer(1))
    if (any(errors)) { #nocov start
      stop ("Error resampling morphy object: ", 
            mpl_translate_error(unique(errors[errors < 0L])))
    }
    if (mpl_apply_tipdata(morphyObj) -> error) {
      stop("Error applying tip data: ", mpl_translate_error(error))
    } #nocov end
    res <- EdgeListSearch(edgeList[1:2], morphyObj, EdgeSwapper = EdgeSwapper,
                          maxIter = searchIter, maxHits = searchHits,
                          verbosity = verbosity - 1L, ...)
    res[1:2]
  }, edgeList)
  
  jackTrees <- structure(apply(jackEdges, 2, function(edgeList) {
    ret <- tree
    ret[["edge"]] <- cbind(edgeList[[1]], edgeList[[2]])
    ret
  }), class = "multiPhylo")
}


#' Label nodes with jackknife support values
#' 
#' @template treeParam
#' @param jackTrees A list or `multiPhylo` object containing trees generated
#' by [`Jackknife()`].
#' @param add Logical specifying whether to add the labels to an existing
#' plot.
#' @param adj,col,frame,pos,\dots Parameters to pass to `nodelabels()`.
#' @param plot Logical specifying whether to plot results; if `FALSE`,
#' returns blank labels for nodes near the root that do not correspond to a
#' unique split.
#' 
#' @return A named vector specifying the proportion of jackknife trees 
#' consistent with each node in `tree`, as plotted.
#' If `plot = FALSE`, blank entries are included corresponding to nodes
#' that do not require labelling; the return value is in the value required
#' by `phylo$node.label`.
#' 
#' @examples
#' library("TreeTools", quietly = TRUE) # for as.phylo
#' 
#' # jackTrees will usually be generated with Jackknife(), but for simplicity:
#' jackTrees <- as.phylo(1:100, 8)
#' 
#' tree <- as.phylo(0, 8)
#' JackLabels(tree, jackTrees)
#' 
#' tree$node.label <- JackLabels(tree, jackTrees, plot = FALSE)
#' @template MRS
#' @importFrom ape nodelabels
#' @importFrom TreeTools SplitFrequency SupportColour
#' @seealso [`Jackknife()`]: Generate trees by jackknife resampling
#' @family split support functions
#' @export
JackLabels <- function (tree, jackTrees,
                        plot = TRUE,
                        add = FALSE,
                        adj = 0, col = NULL, frame = "none", pos = 2L,
                        ...) {
  jackSupport <- SplitFrequency(tree, jackTrees) / length(jackTrees)
  
  if (plot) {
    if (!add) plot(tree)
    if (is.null(col)) col <- SupportColour(jackSupport)
    ape::nodelabels(paste("\n\n", signif(jackSupport, 2)), 
                    node = as.integer(names(jackSupport)),
                    adj = adj, col = col, pos = pos, frame = frame, ...)
    
    # Return:
    jackSupport
  } else {
    ret <- character(tree[["Nnode"]])
    ret[as.integer(names(jackSupport)) - NTip(tree)] <- jackSupport
    
    # Return:
    ret
  }
}
