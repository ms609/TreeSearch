#' Rank taxa by their influence on phylogenetic results
#' 
#' `TaxonInfluence()` ranks taxa according to their influence on the most
#' parsimonious topology.
#' 
#' `TaxonInfluence()` follows the approach of 
#' \insertCite{Mariadassou2012sb;textual}{TreeSearch} in repeating tree search
#' whilst leaving each taxon in turn out of the analysis, and measuring
#' the distance of reconstructed trees from the optimal tree obtained when
#' all taxa are included in phylogenetic inference.
#' 
#' As \insertCite{Denton2018ee;textual}{TreeSearch} emphasize, the
#' Robinson&ndash;Foulds distance is unsuitable for this purpose; this function
#' allows the user to specify a preferred tree distance measure, defaulting
#' to the clustering information distance \insertCite{Smith2020}{TreeSearch}.
#' Because optimal parsimony trees are not equiprobable, taxon influence is
#' ranked based on the maximum and minimum tree-to-tree distances between
#' optimal trees.
#' 
#' @section Distance-weighted mean:
#' Sets of equally parsimonious trees are not statistical samples of tree space,
#' but are biased towards areas of uncertainty.
#' It is possible that a set of trees contains all possible resolutions of a
#' particular clade, and a single other topology in which that clade does not
#' exist &ndash; essentially two distinct solutions, one (_a_) which could be
#' summarised with a summary tree that contains a polytomy, and another (_b_) 
#' which could be summarized by a perfectly resolved tree.
#' Neither of these scenarios is preferable under the principles of parsimony;
#' but summary statistics (e.g. mean, median) will be strongly influenced by the
#' many trees in group _a_, thus underplaying the existence of solution _b_.
#' 
#' `TaxonInfluence()` uses an _ad hoc_ method to produce summary statistics
#' after weighting for trees' distance from other trees.  Trees that have few
#' close neighbours contribute more to the weighted mean, thus reducing the
#' influence of many trees that differ only in small details.
#' This distance-weighted mean is thus less prone to bias than a simple mean
#' &ndash; it is no more statistically valid, but (potentially) provides a more
#' representative summary of comparisons between sets of trees.
#' 
#' 
#' @inheritParams MaximizeParsimony
#' @param tree Optimal tree or summary tree (of class "phylo") or list of trees
#' (of class "list" or "multiPhylo") against which results should be evaluated.
#' If `NULL`, an optimal tree will be sought using parsimony search with 
#' the parameters provided in \code{\dots}.
#' @param Distance Function to calculate tree distance; default:
#' \code{\link[TreeDist:ClusteringInfoDistance]{ClusteringInfoDistance()}}.
#' @param calcWeighted Logical specifying whether to compute the
#' distance-weighted mean value.
#' @param savePath Character giving prefix of path to which reduced trees will
#' be saved (with \code{\link[ape:write.nexus]{write.nexus()}}).
#' File names will follow the pattern
#' `paste0(savePath, droppedTaxonName, ".nex")`; `savePath` should thus contain
#' a trailing `/` if writing to a directory, which will be created if it does
#' not exist.  Special characters will be removed from leaf labels when
#' creating the file path (using 
#' \code{\link[fs:path_sanitize]{path_sanitize()}}).
#' If `NULL`, computed trees will not be saved.
#' @param useCache Logical vector; if `TRUE`, previous tree search results will
#' be loaded from the location given by `savePath`, instead of running a fresh
#' search with the specified dataset and parameters.
#' 
#' @param verbosity,\dots Parameters for [`MaximizeParsimony()`].
#' Tree search will be conducted using `tree` as a starting tree.
#' 
#' @returns `TaxonInfluence()` returns a matrix listing the phylogenetic
#' influence of each taxon, measured in the units of the chosen tree distance
#' metric (default = bits).
#' Columns denote taxa; rows denote the maximum, distance-weighted mean,
#' and minimum distance between optimal tree sets.
#' 
#' @references \insertAllCited{}
#' 
#' @template MRS
#' @examples
#' #' # Load data for analysis in R
#' library("TreeTools")
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' 
#' # Small dataset for demonstration purposes
#' dataset <- congreveLamsdellMatrices[[42]][1:8, ]
#' bestTree <- MaximizeParsimony(dataset, verbosity = 0)[[1]]
#' 
#' # Calculate tip influence
#' influence <- TaxonInfluence(dataset, ratchIt = 0, startIt = 0, verbos = 0)
#' 
#' # Colour tip labels according to their influence
#' upperBound <- 2 * TreeDist::ClusteringEntropy(
#'   PectinateTree(NTip(dataset) - 1))
#' nBin <- 128
#' bin <- cut(
#'   influence["max", ],
#'   breaks = seq(0, upperBound, length.out = nBin),
#'   include.lowest = TRUE
#' )
#' palette <- hcl.colors(nBin, "inferno")
#' 
#' plot(bestTree, tip.color = palette[bin])
#' PlotTools::SpectrumLegend(
#'   "bottomleft",
#'   palette = palette,
#'   title = "Tip influence / bits",
#'   legend = signif(seq(upperBound, 0, length.out = 4), 3),
#'   bty = "n"
#' )
#' @family tree scoring
#' @importFrom ape read.nexus write.nexus
#' @importFrom cli cli_alert_info cli_h1
#' @importFrom fs path_sanitize
#' @importFrom stats weighted.mean
#' @importFrom TreeDist ClusteringInfoDistance
#' @encoding UTF-8
#' @export
TaxonInfluence <- function(
    dataset,
    tree = NULL,
    Distance = ClusteringInfoDistance,
    calcWeighted = TRUE,
    savePath = NULL,
    useCache = FALSE,
    verbosity = 3L,
    ...
  ) {
  if (!is.null(savePath)) {
    saveDir <- dirname(paste0(savePath, "taxon_name"))
    if (!dir.exists(saveDir)) {
      dir.create(saveDir, recursive = TRUE)
    }
  } else if (useCache) {
    stop("Specify cache path using `savePath` parameter")
  }
  
  if (is.null(tree)) {
    tree <- MaximizeParsimony(dataset, ...)
  }
  if (calcWeighted) {
    refWeights <- if (inherits(tree, "phylo") || length(tree) == 1) {
      1
    } else {
      rowSums(as.matrix(Distance(tree)))
    }
  }
  
  startTree <- MakeTreeBinary(if (inherits(tree, "phylo")) {
    tree
  } else {
    tree[[1]]
  })
  if (!inherits(startTree, "phylo")) {
    stop("`tree` must be an object / list of objects of class \"phylo\"")
  }
  
  # Return:
  vapply(names(dataset), function(leaf) {
    
    leafFile <- paste0(savePath, path_sanitize(leaf), ".nex")
    
    result <- if (useCache && file.exists(leafFile)) {
      if (verbosity > 1) {
        cli_alert_info(paste("Seeking cache: ", leafFile))
      }
      tryCatch(c(read.nexus(leafFile)),
               error = function(e) NULL)
    } else {
      NULL
    }
    
    if (is.null(result)) {
      if (verbosity > 0) {
        cli_h1(paste("Taxon influence search:", leaf))
      }
      result <- unique(MaximizeParsimony(
        dataset = dataset[setdiff(names(dataset), leaf)],
        tree = DropTip(startTree, leaf),
        verbosity = verbosity,
        ...
      ))
      if (!is.null(savePath)) {
        write.nexus(result, file = leafFile)
      }
    }
    d <- matrix(Distance(tree, result), length(result))
    dwMean <- if (calcWeighted) {
      resWeights <- if (length(result) > 1) {
        colSums(as.matrix(Distance(result)))
      } else {
        1
      }
      weighted.mean(d, outer(resWeights, refWeights))
    } else {
      NA_real_
    }
    c(min(d), dwMean, max(d))
  }, c(min = 0, dwMean = 0, max = 0))
}
