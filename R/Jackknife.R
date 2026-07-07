#' Jackknife resampling
#'
#' Resample trees using Jackknife resampling, i.e. removing a subset of
#' characters. For standard parsimony, [`Resample()`] is faster; use
#' `Jackknife()` when you need a custom `TreeScorer` or `EdgeSwapper`.
#'
#' @inheritParams Ratchet
#' @param resampleFreq Double between 0 and 1 stating proportion of characters
#' to resample.
#' @param Resampler Function that drops a random subset of characters from
#' `dataset` and returns a rearranged `edgeList`; see §Details for specification.
#' The default, [`JackknifeTree()`], requires `dataset` to be a
#' [`ParsimonyData`][PrepareData] object.
#' @param jackIter Integer specifying number of jackknife iterations to conduct.
#' @return `Jackknife()` returns a list of trees recovered after jackknife
#' iterations.
#' @details
#' Each argument from `Jackknife` is passed forward to `Resampler()`, except
#' `tree`, which is replaced with a two-element list of integer vectors;
#' the first entry lists the parent of each node in `tree`, i.e.
#' `tree$edge[, 1]`; the second entry lists the corresponding child.
#' `Resampler()` must return an equivalent list.
#' 
#' @template MRS
#' @seealso
#' - [`Resample()`]: Jackknife and bootstrap resampling using the C++ search
#'   engine.
#' - [`JackLabels()`]: Label nodes of a tree with jackknife supports.
#' @examples
#' set.seed(0)
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Asher2005"]]
#' trees <- MaximizeParsimony(dataset, inapp = "missing", verbosity = 0)
#' oneTree <- TreeTools::MakeTreeBinary(trees[[1]])
#' jackTrees <- Jackknife(oneTree, PrepareData(dataset),
#'                        jackIter = 10, searchIter = 120,
#'                        inapp = "missing", verbosity = 0)
#' JackLabels(oneTree, jackTrees)
#' @importFrom TreeTools RenumberEdges RenumberTips
#' @family split support functions
#' @family custom search functions
#' @export
Jackknife <- function(tree, dataset,
                      resampleFreq = 2 / 3,
                      InitializeData = NULL,
                      CleanUpData    = NULL,
                      TreeScorer     = EdgeListScore,
                      EdgeSwapper    = TBRSwap,
                      Resampler      = JackknifeTree,
                      jackIter = 5000L, searchIter = 4000L, searchHits = 42L,
                      verbosity = 1L, ...) {
  if (dim(tree[["edge"]])[1] != 2 * tree[["Nnode"]]) {
    stop("tree must be bifurcating; try RootTree() or MakeTreeBinary()")
  }

  tree <- RenumberTips(tree, .SearchTipLabels(dataset))
  edgeList <- tree[["edge"]]
  edgeList <- RenumberEdges(edgeList[, 1], edgeList[, 2])

  if (!is.null(InitializeData) || !is.null(CleanUpData)) {
    .DeprecatedSearchHooks("Jackknife")
    initializedData <- if (is.null(InitializeData)) dataset else
      InitializeData(dataset)
    if (!is.null(CleanUpData)) {
      on.exit(initializedData <- CleanUpData(initializedData))
    }
  } else {
    initializedData <- dataset
  }

  if (verbosity > 10L) { #nocov start
    message(" * Beginning search:")
  } #nocov end

  # Conduct jackIter replicates:
  jackEdges <- vapply(seq_len(jackIter), function (x) {
    if (verbosity > 0L) { #nocov start
      message(" * Jackknife iteration ", x, "/", jackIter)
    } #nocov end
    res <- Resampler(edgeList[1:2], initializedData, resampleFreq = resampleFreq,
                     TreeScorer = TreeScorer, EdgeSwapper = EdgeSwapper,
                     maxIter = searchIter, maxHits = searchHits,
                     verbosity = verbosity, ...)
    res[1:2]
  }, edgeList)

  jackTrees <- structure(apply(jackEdges, 2, function(edgeList) {
    ret <- tree
    ret[["edge"]] <- cbind(edgeList[[1]], edgeList[[2]])
    ret
  }), class = "multiPhylo")
}

#' @describeIn Jackknife Default `Resampler`: drops a random subset of
#' characters from a [`ParsimonyData`][PrepareData] object and conducts a
#' tree search on the reduced dataset.
#' @inheritParams EdgeListSearch
#' @param maxIter Numeric specifying maximum number of iterations to perform
#' in tree search.
#' @param maxHits Numeric specifying maximum number of hits to accomplish in
#' tree search.
#' @export
JackknifeTree <- function (edgeList, dataset, resampleFreq = 2 / 3,
                           TreeScorer = EdgeListScore, EdgeSwapper = NNISwap,
                           maxIter, maxHits, verbosity = 1L, ...) {
  startWeights <- dataset[["original_weight"]]
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

  resampling <- tabulate(sample(deindexedChars, charsToKeep, replace = FALSE),
                         nbins = length(startWeights))
  # R copy-on-modify: the caller's `dataset` is unchanged.
  dataset[["weight"]] <- as.integer(resampling)

  res <- EdgeListSearch(edgeList[1:2], dataset,
                        TreeScorer = TreeScorer, EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        verbosity = verbosity - 1L, ...)

  # Return:
  res[1:2]
}


#' Label nodes with jackknife support values
#' 
#' `JackLabels()` produces a list of node labels denoting split support from
#' a set of resampled trees, optionally printing them on a tree.
#' 
#' If an element of `jackTrees` contains multiple trees, then the iteration is
#' counted as supporting a split if all trees contain the split, and as
#' contradicting the split if no trees contain it.  If a split is only present
#' in a subset of trees, that iteration is considered not to be decisive, and
#' is ignored when calculating the support for that split.
#' 
#' @inheritParams TreeTools::Renumber
#' @param jackTrees A list or `multiPhylo` object containing trees generated
#' by [`Resample()`] or [`Jackknife()`].
#' @param add Logical specifying whether to add the labels to an existing
#' plot.
#' @param adj,col,frame,pos,\dots Parameters to pass to `nodelabels()`.
#' @param plot Logical specifying whether to plot results; if `FALSE`,
#' returns blank labels for nodes near the root that do not correspond to a
#' unique split.
#' @param showFraction Logical specifying whether to also annotate nodes
#' with the fraction of replicates that were decisive for the split.
#' @param format Character specifying return format.
#' `"character"` returns a character string suitable to add to the `node.labels`
#' attribute of a tree; 
#' "numeric" returns numeric values suitable for further analysis.
#' 
#' @return A named vector specifying the proportion of jackknife iterations 
#' consistent with each node in `tree`, as plotted.
#' If `format = "character"`, blank entries are included corresponding to nodes
#' that do not require labels, such that the return value is in the format
#' required by `phylo$node.label`.
#' If multiple trees are specified per iteration, the return value has an
#' attribute `decisive` listing, for each entry in the return value, how many
#' iterations were decisive for that split.
#' 
#' @examples
#' library("TreeTools", quietly = TRUE) # for as.phylo
#' 
#' # jackTrees will usually be generated with Jackknife() or Resample(),
#' # but for simplicity:
#' jackTrees <- as.phylo(1:100, 8)
#' 
#' tree <- as.phylo(0, 8)
#' JackLabels(tree, jackTrees)
#' 
#' tree$node.label <- JackLabels(tree, jackTrees, plot = FALSE)
#' 
#' # Write the labelled tree to screen
#' ape::write.tree(tree)
#'
#' # Write labelled trees to a nexus file:
#' # write.nexus(tree, file = filename)
#' @template MRS
#' @importFrom ape nodelabels
#' @importFrom TreeTools NSplits SplitFrequency SupportColour
#' @seealso
#' Generate trees by jackknife resampling using [`Resample()`] for standard
#' parsimony searches, or [`Jackknife()`] for custom search criteria.
#' @family split support functions
#' @export
JackLabels <- function(tree, jackTrees, plot = TRUE, add = FALSE, adj = 0,
                       col = NULL, frame = "none", pos = 2L,
                       showFraction = FALSE, format = "character", ...) {
  
  nJack <- length(jackTrees)
  multi <- vapply(jackTrees, inherits, TRUE, "multiPhylo")
  
  if (any(multi)) {
    jackTrees[!multi] <- lapply(jackTrees[!multi], c)
    supports <- vapply(jackTrees, function(trees) {
      freq <- SplitFrequency(tree, trees)
      ifelse(freq == 0, FALSE, ifelse(freq == length(trees), TRUE, NA))
      }, logical(NSplits(tree)))
    numerator <- rowSums(supports, na.rm = TRUE)
    denominator <- rowSums(!is.na(supports))
    jackSupport <- structure(numerator / denominator, decisive = denominator)
  } else {
    jackSupport <- SplitFrequency(tree, jackTrees) / nJack
  }
  
  fracText <- if(isTRUE(showFraction)) {
    if (!any(multi)) {
      numerator <- jackSupport * nJack
      denominator <- nJack
    }
    paste0("{", numerator, " / ", denominator, "}")
  } else {
    character(0)
  }
  
  if (plot) {
    if (!add) plot(tree)
    if (is.null(col)) {
      col <- SupportColour(jackSupport)
    }
    nodelabels(paste("\n\n", signif(jackSupport, 2),
                     gsub("{", "(", fixed = TRUE,
                          gsub("}", ")", fixed = TRUE, fracText))),
               node = as.integer(names(jackSupport)),
               adj = adj, col = col, pos = pos, frame = frame, ...)
  }

  numeric <- c("numeric", "number", "double")
  character <- c("character", "text")
  returnMode <- c(rep("numeric", length(numeric)),
                  rep("character", length(character)))[
    pmatch(tolower(format), c(numeric, character), duplicates.ok = TRUE)]
  
  # Return:
  switch(
    returnMode,
    "character" = {
      ret <- character(tree[["Nnode"]])
      idx <- as.integer(names(jackSupport)) - NTip(tree)
      
      ret[idx] <- if (isTRUE(showFraction)) {
        paste(jackSupport, fracText)
      } else {
        jackSupport
      }
      ret
    }, jackSupport
  )
}
