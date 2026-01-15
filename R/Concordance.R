#' Concordance factors
#'
#' Concordance measures the strength of support that characters in a dataset
#' present for each split (=edge/branch) in a tree
#' \insertCite{Minh2020;SmithConc}{TreeSearch}.
#'
# # Renumber before MaximizeParsimony, for `tree`
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]][, 1:20]
#' tree <- TreeSearch::referenceTree
#'
#' cc <- ClusteringConcordance(tree, dataset)
#' mcc <- MutualClusteringConcordance(tree, dataset)
#'
#' qc <- QuartetConcordance(tree, dataset)
#'
#' pc <- PhylogeneticConcordance(tree, dataset)
#' spc <- SharedPhylogeneticConcordance(tree, dataset)
#'
#' oPar <- par(mar = rep(0, 4), cex = 0.8) # Set plotting parameters
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(qc, 3), cex = 0.8)
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(cc, 3), cex = 0.8)
#' par(oPar) # Restore plotting parameters
#'
#' # Write concordance factors to file
#' labels <- paste0(cc, "/", qc, "/", pc) # "/" is a valid delimiter
#' # Identify the node that corresponds to each label
#' whichNode <- match(TreeTools::NTip(tree) + 1:tree$Nnode, names(qc))
#'
#' # The contents of tree$node.label will be written at each node
#' tree$node.label <- labels[whichNode]
#'
#' ape::write.tree(tree) # or write.nexus(tree, file = "mytree.nex")
#'
#' # Display correlation between concordance factors
#' pairs(cbind(cc, mcc, qc, pc, spc), asp = 1)
#' @template MRS
#' @family split support functions
#' @name SiteConcordance
NULL

#' @rdname SiteConcordance
#' @details
#' `ClusteringConcordance()` measures how well each split reflects the
#' grouping information present in each character
#' \insertCite{SmithConc}{TreeSearch}. It treats characters and
#' splits as clusterings of taxa, quantifying their agreement using normalized
#' mutual information. Values range from 1 (perfect agreement) to 0 (no more
#' agreement than expected by chance); negative values indicate less agreement
#' than expected. Summaries returned by `return = "edge"` or `"char"` report,
#' respectively, how well each split is supported by all characters, or how
#' well each character is reflected across all splits.
#'
#' @param return Character specifying whether to summarize support per
#' character (`"char"`) or per edge (`"edge"`). See below for details.
#' @param normalize Logical or numeric; if `TRUE` the mutual information will be
#' normalized such that zero corresponds to the expected mutual information of
#' a randomly drawn character with the same distribution of tokens.
#' If `FALSE`, zero will correspond to zero mutual information,
#' even if this is not achievable in practice.
#' The exact analytical solution, though fast, does not account for
#' non-independence between splits. This limitation is minor for larger
#' trees, and becomes negligible for trees with more than ~200 leaves.
#' For smaller trees, the expected value for random trees can be approximated
#' by resampling relabelled trees. Setting `normalize = n` will approximate the
#' expected value based on _n_ samples.
#'
#' For `return = "char"`, `"tree"`, values will be normalized such that 1
#' corresponds to the maximum possible value, and 0 to the expected value.
#' If `normalize = TRUE`, this will be the expected value for a random
#' character on the given tree. If `normalize` is numeric, the expected value
#' will be estimated by fitting the character to `n` uniformly random trees.
#'
#' @returns
#' `ClusteringConcordance(return = "all")` returns a 3D array where each
#' slice corresponds to a character (site), each column to a tree split, and
#' each row to a different information measure. The `normalized` row gives the
#' normalized mutual information between each split-character pair, scaled so
#' that 1.0 corresponds to `hBest` (the theoretical maximum mutual information,
#' being the minimum of `hSplit` and `hChar`) and 0.0 corresponds to `miRand`
#' (the expected mutual information under random association). `hSplit` gives
#' the entropy (information content) of each split's bipartition; `hChar` gives
#' the entropy of each character's state distribution; `hJoint` gives the joint
#' entropy of the split-character confusion matrix; `mi` gives the raw mutual
#' information; and `n` records the number of informative observations.
#' Negative normalized values indicate observed mutual information below random
#' expectation. `NA` is returned when `hBest = 0` (no information potential).
#'
#' `ClusteringConcordance(return = "edge")` returns a vector where each element
#' corresponds to a split (an edge of the tree) and gives the normalized mutual
#' information between that split and the character data, averaged across all
#' characters.
#' When `normalize = TRUE` (default), values are scaled relative to random
#' expectation; when `FALSE`, raw mutual information normalized by `hBest` is
#' returned.
#'
#' `ClusteringConcordance(return = "char")` returns a vector where each element
#' corresponds to a character (site) and gives the entropy-weighted average
#' normalized mutual information between that character and all tree splits.
#' Characters with higher information content receive proportionally more weight
#' from splits that can potentially convey more information about them.
#'
#' `ClusteringConcordance(return = "tree")` returns a single value representing
#' the overall concordance between the tree topology and the character data.
#' This averages the fit of the best-matching split for each character.
#' This is included for completeness, though it is not clear that this is a useful
#' or meaningful measure.
#
# I had previously considered calculating the entropy-weighted average of normalized
# mutual information across all split-character pairs, where each pair contributes
# proportionally to its potential information content.
# The problem here is that imperfect matches between compatible splits
# come to dominate, resulting in a small score that gets smaller as trees get
# larger, even with a perfect fit.
#'
#'
#' @seealso
#' - [Consistency()]
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ClusteringConcordance(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom abind abind
#' @importFrom stats setNames
#' @importFrom TreeDist ClusteringEntropy Entropy entropy_int
#' MutualClusteringInfo
#' @importFrom TreeTools as.Splits MatchStrings Subsplit TipLabels
#' @export
ClusteringConcordance <- function(
  tree,
  dataset,
  return = "edge",
  normalize = TRUE
) {
  # Check inputs
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (is.null(tree)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }

  keep <- MatchStrings(TipLabels(tree), names(dataset), warning)
  if (length(keep) == 0) {
    return(NULL)
  }
  dataset <- dataset[keep]

  # Prepare data
  splits <- as.logical(as.Splits(tree))

  at <- attributes(dataset)
  cont <- at[["contrast"]]
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1

  mat <- matrix(as.integer(unlist(dataset)), length(dataset), byrow = TRUE)
  mat[mat %in% which(ambiguous)] <- NA_integer_
  maxToken <- max(mat, na.rm = TRUE)
  tokens <- as.character(seq_len(maxToken))
  mat <- apply(mat, 2, function(x) {
    uniques <- tabulate(x, maxToken) == 1
    x[x %in% tokens[uniques]] <- NA_integer_
    x
  })

  # Calculate entropy
  h <- simplify2array(apply(mat, 2, function(char) {
    aChar <- !is.na(char)
    ch <- char[aChar]
    if (length(ch) == 0) {
      # All ambiguous
      n <- 0
      hChar <- 0
    } else {
      chMax <- max(1, ch)
      chTable <- tabulate(ch, chMax)
      n <- length(ch)
      hChar <- entropy_int(chTable)
    }

    hh <- apply(splits[, aChar, drop = FALSE], 1, function (spl) {
      spTable <- tabulate(spl + 1, 2)
      if (any(spTable < 2)) {
        c(hSplit = 0,
          hJoint = hChar,
          miRand = 0,
          n = n)
      } else {
        c(hSplit = entropy_int(spTable),
          hJoint = entropy_int(tabulate(ch + (spl * chMax), chMax + chMax)),
          miRand = .ExpectedMI(spTable, chTable),
          n = n)
      }
    })
    
    rbind(hChar = hChar, hh)
  }, simplify = FALSE))
  
  if (length(dim(h)) == 2) {
    # Matrix to 3D array
    dim(h) <- c(dim(h), 1)
  }

  h[abs(h) < sqrt(.Machine$double.eps)] <- 0
  hh <- h[, , at[["index"]], drop = FALSE]
  
  hBest <- `rownames<-`(pmin(hh["hChar", , , drop = FALSE],
                             hh["hSplit", , , drop = FALSE]), NULL)
  mi <- `rownames<-`(hh["hChar", , , drop = FALSE] + 
                       hh["hSplit", , , drop = FALSE] -
                       hh["hJoint", , , drop = FALSE], NULL)
  miRand <- `rownames<-`(hh["miRand", , , drop = FALSE], NULL)
  norm <- if (isFALSE(normalize)) {
    ifelse(hBest == 0, NA, mi / hBest)
  } else {
    ifelse(hBest == 0, NA, .Rezero(mi / hBest, miRand / hBest))
  }

  returnType <- pmatch(tolower(return), c("all", "edge", "character", "tree"),
                       nomatch = 1L)
  if (returnType %in% 3:4) { # character / tree
    charSplits <- apply(mat, 2, simplify = FALSE, function(x)
      as.Splits(x[!is.na(x)], tipLabels = keep[!is.na(x)]))
    charMax <- vapply(charSplits, ClusteringEntropy, double(1))[
      attr(dataset, "index")]
    charInfo <- MutualClusteringInfo(tree, charSplits)[at[["index"]]]
    if (is.numeric(normalize)) {
      rTrees <- replicate(normalize, RandomTree(tree), simplify = FALSE)
      randInfo <- MutualClusteringInfo(rTrees, charSplits)[, attr(dataset, "index")]
      randMean <- colMeans(randInfo)
      var <- rowSums((t(randInfo) - randMean) ^ 2) / (normalize - 1)
      mcse <- sqrt(var / normalize)
      randTreeInfo <- rowSums(randInfo)
      randTreeMean <- mean(randTreeInfo)
      treeVar <- var(randTreeInfo)
      mcseTree <- sqrt(treeVar / normalize)
    }
  }

  # Return:
  switch(returnType,
         # all
         abind(
           along = 1,
           normalized = norm,
           hh,
           hBest = hBest,
           mi = mi,
           miRand = miRand
         ), {
           # edge
           best <- rowSums(hBest[1, , , drop = FALSE], dims = 2)
           ifelse(!is.na(best) & best == 0,
                  NA_real_,
                  if (isTRUE(normalize)) {
                    .Rezero(
                      rowSums(mi[1, , , drop = FALSE], dims = 2) / best,
                      rowSums(miRand[1, , , drop = FALSE], dims = 2) / best
                    )
                  } else {
                    rowSums(mi[1, , , drop = FALSE], dims = 2) / best
                  })[1, ]
         }, {
           # char
         
           # one <- hh["hChar", 1, , drop = TRUE] # All rows equal
           one <- charMax
           zero <- if (isFALSE(normalize)) {
             0
           } else if (isTRUE(normalize)) {
             apply(hh["miRand", , ], 2, max)
           } else {
             randMean
           }
           ret <- (charInfo - zero) / (one - zero)
           if (is.numeric(normalize)) {
             mcseInfo <- ((one - charInfo) / (one - zero) ^ 2) * mcse
             mcseInfo[mcseInfo < sqrt(.Machine$double.eps)] <- 0
             structure(ret, hMax = charMax, mcse = mcseInfo)
           } else {
             structure(charInfo / charMax, hMax = charMax)
           }
         }, {
           # tree: Entropy-weighted mean across best character-split pairs
           warning("I'm not aware of a situation in which this is a useful measure.")
           # Random normalization doesn't work if we pick the best matching
           # split; one will be randomly better than another, even in the
           # absence of signal.
           norm[is.na(norm)] <- 0
           bestMatch <- apply(norm, 3, which.max)
           idx <- cbind(1, bestMatch, seq_along(bestMatch))
           return(weighted.mean(norm[idx], hBest[idx]))
          
           if (isFALSE(normalize)) {
           } else {
             one <- sum(hh["hChar", 1, ])
             zero <- if (isTRUE(normalize)) {
               sum(apply(hh["miRand", , ], 2, max))
             } else {
               randTreeMean
             }
             ret <- (sum(charInfo) - zero) / (one - zero)
             if (is.numeric(normalize)) {
               mcseInfo <- ((one - sum(charInfo)) / (one - zero) ^ 2) * mcseTree
               mcseInfo[mcseInfo < sqrt(.Machine$double.eps)] <- 0
               structure(ret, mcse = mcseInfo)
             } else {
               ret
             }
           }
         }
  )
}

#' Generate colour to depict the amount and quality of observations
#' @param amount Numeric vector of values between 0 and 1, denoting the relative
#' amount of information
#' @param quality Numeric vector of values between -1 and 1, denoting the
#' quality of observations, where 0 is neutral.
#' @return `QACol()` returns an RGB hex code for a colour, where lighter colours
#' correspond to entries with a higher `amount`; unsaturated colours denote
#' a neutral `quality`; and red/cyan colours denote low/high `quality`.
#' @examples
#' amount <- runif(80, 0, 1)
#' quality <- runif(80, -1, 1)
#' plot(amount, quality, col = QACol(amount, quality), pch = 15)
#' abline(h = 0)
#' @template MRS
#' @importFrom colorspace hex polarLUV
#' @export
QACol <- function(amount, quality) {
  h <- 80 + (quality * 140)
  l <- amount * 88 # < 100: white can take no hue
  c <- abs(quality) * .MaxChroma(h, l)
  # Saturation higher than 1 risks overflowing the colour space
  # Small overflows are caught via `fixup = TRUE`; large overflows will produce
  # bright red errors
  saturation <- 0.999 # Safe if max_chroma(floor = FALSE) slightly overestimates
  saturation <- 1.16
  
  hex(polarLUV(
    H = as.numeric(h),
    C = as.numeric(c) * saturation,
    L = as.numeric(l)
  ), fixup = TRUE)
}

#' @importFrom colorspace max_chroma
.MaxChroma <- function(h, l) {
  ret <- `length<-`(double(0), length(h))
  applicable <- !is.na(h) & !is.na(l)
  ret[applicable] <- max_chroma(h[applicable], l[applicable])
  ret
}

#' @rdname QACol
#' @return `QCol()` returns an RGB hex code for a colour, where darker,
#' unsaturated colours denote a neutral `quality`;
#' and red/cyan colours denote low/high `quality`. `amount` is ignored.
#' @export
QCol <- function(amount, quality) {
  h <- 80 + (quality * 140)
  l <- abs(quality) * 88 # < 100: white can take no hue
  c <- abs(quality) * .MaxChroma(h, l)
  # Saturation higher than 1 risks overflowing the colour space
  # Small overflows are caught via `fixup = TRUE`; large overflows will produce
  # bright red errors
  saturation <- 0.999 # Safe if max_chroma(floor = FALSE) slightly overestimates
  
  hex(polarLUV(
    H = as.numeric(h),
    C = as.numeric(c) * saturation,
    L = as.numeric(l)
  ), fixup = TRUE)
}

#' @rdname QACol
#' @param where Location of legend, passed to `par(fig = where)`
#' @param n Integer vector giving number of cells to plot in swatch for
#' `quality` and `amount`.
#' @inheritParams ConcordanceTable
#' @export
QALegend <- function(where = c(0.1, 0.3, 0.1, 0.3), n = 5, Col = QACol) {
  oPar <- par(fig = where, new = TRUE, mar = rep(0, 4), xpd = NA)
  on.exit(par(oPar))
  n <- rep(n, length.out = 2)
  nA <- n[[2]]
  nQ <- n[[1]]
  amount <- seq(0, 1, length.out = nA)
  quality <- seq(-1, 1, length.out = nQ)
  mat <- outer(amount, quality,
               Vectorize(function (a, q) Col(a, q)))
  image(x = amount, y = quality, z = matrix(1:prod(n), nA, nQ),
        col = mat, axes = FALSE, xlab = "", ylab = "")
  mtext("Amount \U2192", side = 1, line = 1)
  mtext("Quality \U2192", side = 2, line = 1)
}

#' Plot concordance table
#'
#' `ConcordanceTable()` plots a concordance table
#' \insertCite{SmithConc}{TreeSearch}.
#'
#' @inheritParams ClusteringConcordance
#' @param Col Function that takes vectors `amount` and `quality` and returns
#' a vector of colours. [QCol] colours by data quality (concordance);
#' [QACol] by quality and amount of information.
#' @param largeClade Integer; if greater than 1, vertical lines will be drawn
#' at edges whose descendants are both contain more than `largeClade` leaves.
#' @param xlab Character giving a label for the x axis.
#' @param ylab Character giving a label for the y axis.
#' @param \dots Arguments to `abline`, to control the appearance of vertical
#' lines marking important edges.
#' @returns `ConcordanceTable()` invisibly returns an named list containing:
#' - `"info"`: The amount of information in each character-edge pair, in bits;
#' - `"relInfo"`: The information, normalized to the most information-rich pair;
#' - `"quality"`: The normalized mutual information of the pair;
#' - `"col"`: The colours used to plot the table.
#'
#' @references \insertAllCited{}
#' @examples
#' # Load data and tree
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]][, 1:20]
#' tree <- referenceTree
#'
#' # Plot tree and identify nodes
#' library("TreeTools", quietly = TRUE)
#' plot(tree)
#' nodeIndex <- as.integer(rownames(as.Splits(tree)))
#' nodelabels(seq_along(nodeIndex), nodeIndex, adj = c(2, 1),
#'            frame = "none", bg = NULL)
#' QALegend(where = c(0.1, 0.4, 0.1, 0.3))
#'
#' # View information shared by characters and edges
#' ConcordanceTable(tree, dataset, largeClade = 3, col = 2, lwd = 3)
#' axis(1)
#' axis(2)
#'
#' # Visualize dataset
#' image(t(`mode<-`(PhyDatToMatrix(dataset), "numeric")), axes = FALSE,
#'       xlab = "Leaf", ylab = "Character")
#' @importFrom graphics abline image mtext
#' @importFrom TreeTools CladeSizes NTip
#' @export
ConcordanceTable <- function(tree, dataset, Col = QACol, largeClade = 0,
                             xlab = "Edge", ylab = "Character", ...) {
  cc <- ClusteringConcordance(tree, dataset, return = "all")
  nodes <- seq_len(dim(cc)[[2]])
  info <- cc["hBest", , ] * cc["n", , ]
  amount <- info / max(info, na.rm = TRUE)
  amount[is.na(amount)] <- 0
  quality <- cc["normalized", , ]
  # Plot points with incalculable quality as black, not transparent.
  amount[is.na(quality)] <- 0
  quality[is.na(quality)] <- 0

  col <- matrix(Col(amount, quality), dim(amount)[[1]], dim(amount)[[2]])
  image(nodes, seq_len(dim(cc)[[3]]),
        matrix(1:prod(dim(amount)), dim(amount)[[1]]),
        frame.plot = FALSE, axes = FALSE,
        col = col, xlab = xlab, ylab = ylab)
  
  if (largeClade > 1) {
    cladeSize <- CladeSizes(tree)
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    bigNode <- vapply(as.integer(colnames(cc)), function (node) {
      all(cladeSize[child[parent == parent[child == node]]] >= largeClade)
    }, logical(1))
    abline(v = nodes[bigNode] - 0.5, ...)
  }
  invisible(list(info = info, amount = amount, quality = quality, col = col))
}

#' @rdname SiteConcordance
#' @details
#' `MutualClusteringConcordance()` provides a character‑wise summary that
#' emphasises each character’s best‑matching split(s). It treats each character
#' as a simple tree and computes the mutual clustering information between this
#' character‑tree and the supplied phylogeny. High values identify characters
#' whose signal is well represented anywhere in the tree, even if concentrated
#' on a single edge.
#'
#' @return `MutualClusteringConcordance()` returns the mutual clustering
#' concordance of each character in `dataset` with `tree`.
#' The attribute `weighted.mean` gives the mean value, weighted by the
#' information content of each character.
#' @importFrom TreeTools MatchStrings
#' @importFrom TreeDist ClusteringEntropy MutualClusteringInfo
#' @export
MutualClusteringConcordance <- function(tree, dataset) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }

  dataset <- dataset[MatchStrings(TipLabels(tree), names(dataset))]
  splits <- as.multiPhylo(as.Splits(tree))
  characters <- as.multiPhylo(dataset)
  
  support <- rowSums(vapply(characters, function (char) {
    trimmed <- KeepTip(splits, TipLabels(char))
    cbind(mi = MutualClusteringInfo(char, trimmed),
          possible = ClusteringEntropy(trimmed))
  }, matrix(NA_real_, length(splits), 2)), dims = 2)
  
  ret <- support[, 1] / support[, 2]
  # Return:
  structure(ret, weighted.mean = weighted.mean(ret, support[, 2]))
}

#' @rdname SiteConcordance
#' @details
#' `QuartetConcordance()` is the proportion of quartets (sets of four leaves)
#' that are decisive for a split which are also concordant with it
#' (the site concordance factor \insertCite{Minh2020}{TreeSearch}).
#' For example, a quartet with the characters `0 0 0 1` is not decisive, as
#' all relationships between those leaves are equally parsimonious.
#' But a quartet with characters `0 0 1 1` is decisive, and is concordant
#' with any tree that groups the first two leaves together to the exclusion
#' of the second.
#'
#' By default, the reported value weights each site by the number of quartets
#' it is decisive for.  This value can be interpreted as the proportion of
#' all decisive quartets that are concordant with a split.
#' If `weight = FALSE`, the reported value is the mean of the concordance
#' value for each site.
#' Consider a split associated with two sites:
#' one that is concordant with 25% of 96 decisive quartets, and
#' a second that is concordant with 75% of 4 decisive quartets.
#' If `weight = TRUE`, the split concordance will be 24 + 3 / 96 + 4 = 27%.
#' If `weight = FALSE`, the split concordance will be mean(75%, 25%) = 50%.
#'
#' `QuartetConcordance()` is computed exactly, using all quartets, where as
#' other implementations (e.g. IQ-TREE) follow
#' \insertCite{@Minh2020;textual}{TreeSearch} in using a random subsample
#'  of quartets for a faster, if potentially less accurate, computation.
#' Ambiguous and inapplicable tokens are treated as containing no grouping
#' information (i.e. `(02)` or `-` are each treated as `?`).
#' @return
#' `QuartetConcordance(return = "edge")` returns a numeric vector giving the
#' concordance index at each split across all sites; names specify the number of
#' each corresponding split in `tree`.
#'
#' `QuartetConcordance(return = "char")` returns a numeric vector giving the
#' concordance index calculated at each site, averaged across all splits.
#'
#' @param weight Logical specifying whether to weight sites according to the
#' number of quartets they are decisive for.
#' @importFrom ape keep.tip
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom utils combn
#' @importFrom TreeTools as.Splits PhyDatToMatrix TipLabels
#' @export
QuartetConcordance <- function(
  tree,
  dataset = NULL,
  weight = TRUE,
  return = "edge"
) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }
  tipLabels <- intersect(TipLabels(tree), names(dataset))
  if (!length(tipLabels)) {
    warning("No overlap between tree labels and dataset.")
    return(NULL)
  }
  dataset <- dataset[tipLabels, drop = FALSE]
  splits <- as.Splits(tree, dataset)
  logiSplits <- vapply(seq_along(splits), function (i) as.logical(splits[[i]]),
                       logical(NTip(dataset)))
  
  characters <- PhyDatToMatrix(dataset, ambigNA = TRUE)
  charLevels <- attr(dataset, "allLevels")
  isAmbig <- rowSums(attr(dataset, "contrast")) > 1
  isInapp <- charLevels == "-"
  nonGroupingLevels <- charLevels[isAmbig | isInapp]
  characters[characters %in% nonGroupingLevels] <- NA

  charInt <- `mode<-`(characters, "integer")
  raw_counts <- quartet_concordance(logiSplits, charInt)

  num <- raw_counts$concordant
  den <- raw_counts$decisive
  options <- c("character", "site", "default")
  return <- options[[pmatch(tolower(trimws(return)), options,
                            nomatch = length(options))]]
  

  if (return == "default") {
    if (isTRUE(weight)) {
      # Sum numerator and denominator across sites (columns), then divide
      # This matches weighted.mean(num/den, den) == sum(num) / sum(den)
      split_sums_num <- rowSums(num)
      split_sums_den <- rowSums(den)
      ret <- ifelse(
        split_sums_den == 0,
        NA_real_,
        split_sums_num / split_sums_den
      )
    } else {
      # Mean of ratios per site
      # Avoid division by zero (0/0 -> NaN -> NA handled by na.rm)
      ratios <- num / den
      # Replace NaN/Inf with NA for rowMeans calculation
      ratios[!is.finite(ratios)] <- NA
      ret <- rowMeans(ratios, na.rm = TRUE)
    }

    setNames(ret, names(splits))
  } else {
    # return = "char"
    p <- num / den
    if (isTRUE(weight)) {
      vapply(
        seq_len(dim(num)[[2]]),
        function(i) {
          weighted.mean(num[, i] / den[, i], den[, i])
        },
        double(1)
      )
    } else {
      vapply(
        seq_len(dim(num)[[2]]),
        function(i) {
          mean(num[den[, i] > 0, i] / den[den[, i] > 0, i])
        },
        double(1)
      )
    }
  }
}

#' @importFrom fastmap fastmap
.ExpectedMICache <- fastmap()

# @param a must be a vector of length <= 2
# @param b may be longer
#' @importFrom base64enc base64encode
.ExpectedMI <- function(a, b) {
  if (length(a) < 2 || length(b) < 2) {
    0
  } else {
    key <- base64enc::base64encode(mi_key(a, b))
    if (.ExpectedMICache$has(key)) {
      .ExpectedMICache$get(key)
    } else {
      ret <- expected_mi(a, b)

      # Cache:
      .ExpectedMICache$set(key, ret)
      # Return:
      ret
    }
  }
}

#' Re-zero a value by normalization
#' @param value value ranging from zero to one
#' @param zero new value to set as zero
#' @keywords internal
.Rezero <- function(value, zero) {
  (value - zero) / (1 - zero)
}

#' @rdname SiteConcordance
#' @importFrom TreeTools as.multiPhylo CladisticInfo CompatibleSplits
#' MatchStrings
#' @export
PhylogeneticConcordance <- function(tree, dataset) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  dataset <- dataset[MatchStrings(TipLabels(tree), names(dataset))]
  splits <- as.Splits(tree)
  if (is.null(names(splits))) {
    names(splits) <- paste0("sp", seq_along(splits))
  }
  characters <- as.multiPhylo(dataset)
  
  blankRet <- matrix(0, length(splits), 2,
                     dimnames = list(names(splits),
                                     c("concordant", "possible")))
  
  support <- rowSums(vapply(characters, function (char) {
    ret <- blankRet
    if (NTip(char) > 3L) {
      thinned <- Subsplit(splits, TipLabels(char))
      compatible <- CompatibleSplits(thinned, char)
      if (length(compatible)) {
        ci <- CladisticInfo(thinned)
        ret[names(thinned), "concordant"] <- ci * apply(compatible, 1, all)
        ret[names(thinned), "possible"] <- ci
      }
    }
    # Return:
    ret
  }, blankRet), dims = 2)
  
  # Return:
  support[, 1] / support[, 2]
}

#' @rdname SiteConcordance
#' @details
#' `SharedPhylogeneticConcordance()` treats each character as a simple tree.
#' Each token in the character corresponds to a node whose pendant edges are the
#' taxa with that token.
#' The Shared Phylogenetic Concordance for each character in `dataset` is then
#' the Shared Phylogenetic Information \insertCite{Smith2020}{TreeSearch} of
#' this tree and `tree`.
#' @return `SharedPhylogeneticConcordance()` returns the shared phylogenetic
#' concordance of each character in `dataset` with `tree`.
#' The attribute `weighted.mean` gives the mean value, weighted by the
#' information content of each character.
#' @importFrom TreeTools as.multiPhylo MatchStrings
#' @importFrom TreeDist ClusteringInfo SharedPhylogeneticInfo
#' @export
SharedPhylogeneticConcordance <- function(tree, dataset) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  dataset <- dataset[MatchStrings(TipLabels(tree), names(dataset))]
  splits <- as.multiPhylo(as.Splits(tree))
  characters <- as.multiPhylo(dataset)
  
  support <- rowSums(vapply(characters, function (char) {
    trimmed <- lapply(splits, keep.tip, TipLabels(char))
    cbind(mi = SharedPhylogeneticInfo(char, trimmed),
          possible = ClusteringInfo(trimmed))
  }, matrix(NA_real_, length(splits), 2)), dims = 2)
  
  ret <- support[, 1] / support[, 2]
  # Return:
  structure(ret, weighted.mean = weighted.mean(ret, support[, 2]))
}

#' Evaluate the concordance of information between a tree and a dataset
#'
#' Details the amount of information in a phylogenetic dataset that is
#' consistent with a specified phylogenetic tree, and the signal:noise
#' ratio of the character matrix implied if the tree is true.
#'
#' Presently restricted to datasets whose characters contain a maximum of
#' two parsimony-informative states.
#'
#' @return `ConcordantInformation()` returns a named vector with elements:
#'
#' - `informationContent`: cladistic information content of `dataset`
#' - `signal`, `noise`: amount of cladistic information that represents
#' phylogenetic signal and noise, according to `tree`
#' - `signalToNoise`: the implied signal:noise ratio of `dataset`
#' - `treeInformation`: the cladistic information content of a bifurcating tree
#' on `dataset`; this is the minimum amount of information necessary to resolve
#' a bifurcating tree, assuming no duplicate information or noise
#' - `matrixToTree`: the ratio of the cladistic information content of the
#' matrix to the cladistic information content of the tree, a measure of the
#' redundancy of the matrix
#' - `ignored`: information content of characters whose signal and noise could
#' not be calculated (too many states) and so are not included in the totals
#' above.
#'
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ConcordantInformation(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom TreeTools Log2UnrootedMult Log2Unrooted MatchStrings
#' @export
ConcordantInformation <- function(tree, dataset) {
  dataset <- dataset[MatchStrings(TipLabels(tree), names(dataset))]
  originalInfo <- sum(apply(PhyDatToMatrix(dataset), 2, CharacterInformation))
  dataset <- PrepareDataProfile(dataset)

  extraSteps <- CharacterLength(tree, dataset, compress = TRUE) -
    MinimumLength(dataset, compress = TRUE)
  chars <- matrix(unlist(dataset), attr(dataset, "nr"))
  ambiguousToken <- which(attr(dataset, "allLevels") == "?")
  asSplits <- apply(chars, 1, function(x) {
    ret <- table(x)
    if (length(ambiguousToken) != 0) {
      ret[names(ret) != ambiguousToken]
    } else {
      ret
    }
  })
  if (is.matrix(asSplits)) {
    asSplits <- lapply(seq_len(dim(asSplits)[2]), function(i) asSplits[, i])
  }
  ic <- vapply(asSplits, function (split) 
    Log2Unrooted(sum(split)) - Log2UnrootedMult(split),
    double(1))
  
  infoLosses <- apply(chars, 1, StepInformation, 
                      ambiguousToken = ambiguousToken) # , drop = FALSE
  if (is.matrix(infoLosses)) {
    infoLosses <- lapply(seq_len(dim(infoLosses)[2]),
                         function (i) infoLosses[, i])
  }
  
  signal <- vapply(seq_along(extraSteps), function (i) {
    infoLosses[[i]][extraSteps[i] + 1L]
  }, double(1))
  noise <- ic - signal
  noise[noise < sqrt(.Machine[["double.eps"]])] <- 0

  index <- attr(dataset, "index")
  if (any(is.na(signal))) {
    na <- is.na(signal)
    icA <- ic
    icA[na] <- 0
    totalInfo <- sum(ic[index])
    kept <- sum(icA[index])
    discarded <- totalInfo - kept
    warning("Could not calculate signal for characters ",
            paste0(match(which(na), index), collapse = ", "),
            "; discarded ", signif(discarded), " bits from totals.")
    totalNoise <- sum(noise[index], na.rm = TRUE)
    totalSignal <- sum(signal[index], na.rm = TRUE)
    signalNoise <- totalSignal / totalNoise

    infoNeeded <- Log2Unrooted(length(dataset))
    infoOverkill <- totalInfo / infoNeeded
    
    message("`dataset` contains ",
            signif(totalInfo), " bits (after discarding ",
            signif(discarded), "), of which ",
            signif(totalSignal), " signal, ",
            signif(totalNoise), " noise, ",
            signif(infoNeeded), " needed.  ",
            "S:N = ", signif(signalNoise), "\n")
    
  } else {
    totalInfo <- sum(ic[index])
    totalNoise <- sum(noise[index])
    totalSignal <- sum(signal[index])
    signalNoise <- totalSignal / totalNoise
    discarded = 0

    infoNeeded <- Log2Unrooted(length(dataset))
    infoOverkill <- totalInfo / infoNeeded
    discarded <- originalInfo - totalInfo
    if (discarded < sqrt(.Machine[["double.eps"]])) discarded <- 0
    
    message("dataset contains ",
            signif(totalInfo), " bits",
            if (totalInfo != originalInfo) {
              paste0(" (after discarding ", signif(originalInfo - totalInfo),
                     " bits)")
            }, ", of which ", 
            signif(totalSignal), " signal, ",
            signif(totalNoise), " noise, ",
            signif(infoNeeded), " needed.  ",
            "S:N = ", signif(signalNoise), "\n")
  }

  # Return:
  c(
    informationContent = totalInfo,
    signal = totalSignal,
    noise = totalNoise,
    signalToNoise = signalNoise,

    treeInformation = infoNeeded,
    matrixToTree = infoOverkill,
    ignored = discarded
  )
}

#' @rdname ConcordantInformation
#' @export
Evaluate <- function(tree, dataset) {
  .Deprecated("ConcordantInformation()")
  ConcordantInformation(tree, dataset)
}

#' @rdname ConcordantInformation
#' @export
ConcordantInfo <- ConcordantInformation
