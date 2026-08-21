#' Concordance factors
#'
#' Concordance measures the strength of support that characters in a dataset
#' present for each split (=edge/branch) in a tree
#' \insertCite{Minh2020,SmithConc}{TreeSearch}.
#'
# # Renumber before MaximizeParsimony, for `tree`
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \dontshow{
#'   knitr::opts_chunk$set(fig.width = 6, fig.height = 6)
#'   oPar <- par(mar = rep(0, 4), cex = 0.8)
#' }data("congreveLamsdellMatrices", package = "TreeSearch")
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
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(qc, 3), cex = 0.8)
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(cc, 3), cex = 0.8)
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
#' \dontshow{par(oPar)
#' }#' # Display correlation between concordance factors
#' pairs(cbind(cc, mcc, qc, pc, spc), asp = 1)
#' @template MRS
#' @family split support functions
#' @name SiteConcordance
NULL

#' @rdname SiteConcordance
#' @details
#' `ClusteringConcordance()` measures how well each tree split reflects the
#' grouping structure implied by each character
#' \insertCite{SmithConc}{TreeSearch}.
#' Characters and splits are treated as clusterings of taxa, and their agreement
#' is quantified using mutual information (MI). All reported values are scaled 
#' so that 1 corresponds to the maximum possible mutual information for each
#' split–character pair (`hBest`).
#'
#'
#' @param return Character specifying the summary to return. Options are:
#'   - `"edge"`: average concordance of each tree split across all characters;
#'   - `"char"`: concordance of each character averaged over all splits, weighted
#'     by each character's information content;
#'   - `"tree"`: an overall tree‑level concordance score;
#'   - `"all"`: a full array of MI components and normalized values for every
#'     split–character pair.
#'
#' @param chanceCorrect Sets the zero point of the scale.
#' If `FALSE`, zero corresponds to zero MI.
#' If `TRUE`, zero is the value expected when each character's tokens are
#' reassigned at random across the leaves, holding its state frequencies and the
#' split sizes fixed: `QuartetConcordance()` computes this expectation exactly,
#' whereas `ClusteringConcordance()` approximates it, accurately for large trees
#' (~200+ taxa) but neglecting correlation between splits.
#' If a positive integer `n`, the expectation is sampled over `n` random reassignments
#' of each character's tokens (`QuartetConcordance()`), or against `n` uniformly random
#' trees (`ClusteringConcordance()`).
#' _Hint: Clamp chance-corrected values to \eqn{[-1, 1]} before plotting with [QCol()] /
#' [QACol()]_.
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
#' When `chanceCorrect = TRUE` (default), values are scaled relative to random
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
#' - [ConcordanceTable()]: Visualize concordance.
#' - [Consistency()]: Measure consistency of characters with dataset.
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ClusteringConcordance(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom abind abind
#' @importFrom stats setNames
#' @importFrom TreeDist ClusteringEntropy Entropy entropy_int
#' @importFrom TreeDist MutualClusteringInfo
#' @importFrom TreeTools as.Splits MatchStrings Subsplit TipLabels
#' @export
ClusteringConcordance <- function(
  tree,
  dataset,
  return = "edge",
  chanceCorrect = TRUE
) {
  # Check inputs
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (is.null(tree)) {
    warning("Cannot calculate concordance without `tree`.")
    return(NULL)
  }

  keep <- MatchStrings(TipLabels(tree), names(dataset), warning)
  if (length(keep) == 0) {
    return(NULL)
  }
  dataset <- dataset[keep]

  # Prepare data
  # `tree` may carry tips absent from `dataset`; Subsplit restricts `splits` to the
  # shared taxa without renumbering nodes.
  splits <- as.logical(Subsplit(as.Splits(tree), keep))
  # Recover each surviving split's original node number
  if (is.null(rownames(splits))) {
    fullRestricted <- as.logical(as.Splits(tree))[, TipLabels(tree) %in% keep,
                                                  drop = FALSE]
    rownames(splits) <- vapply(seq_len(nrow(splits)), function(i) {
      row <- splits[i, ]
      hit <- apply(fullRestricted, 1, function(r) all(r == row) || all(r != row))
      rownames(fullRestricted)[which(hit)[1]]
    }, character(1))
  }

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
  norm <- if (isFALSE(chanceCorrect)) {
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
    if (is.numeric(chanceCorrect)) {
      rTrees <- replicate(chanceCorrect, RandomTree(tree), simplify = FALSE)
      # Score each random tree against `charSplits` separately: characters with
      # ambiguous tokens yield splits over different tip subsets, and the
      # vectorised `MutualClusteringInfo(<list of trees>, <list of splits>)`
      # path cannot reconcile a single label set across them ("Old and new
      # labels must match"). Looping one tree at a time mirrors the working
      # `charInfo` call above.
      randInfo <- t(vapply(
        rTrees,
        function(rt) MutualClusteringInfo(rt, charSplits),
        double(length(charSplits))
      ))[, attr(dataset, "index"), drop = FALSE]
      randMean <- colMeans(randInfo)
      var <- rowSums((t(randInfo) - randMean) ^ 2) / (chanceCorrect - 1)
      mcse <- sqrt(var / chanceCorrect)
      randTreeInfo <- rowSums(randInfo)
      randTreeMean <- mean(randTreeInfo)
      treeVar <- var(randTreeInfo)
      mcseTree <- sqrt(treeVar / chanceCorrect)
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
                  if (isTRUE(chanceCorrect)) {
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
           zero <- if (isFALSE(chanceCorrect)) {
             0
           } else if (isTRUE(chanceCorrect)) {
             apply(.ConcSlice(hh, "miRand"), 2, max)
           } else {
             randMean
           }
           ret <- (charInfo - zero) / (one - zero)
           if (is.numeric(chanceCorrect)) {
             mcseInfo <- ((one - charInfo) / (one - zero) ^ 2) * mcse
             mcseInfo[mcseInfo < sqrt(.Machine$double.eps)] <- 0
             structure(ret, hMax = charMax, mcse = mcseInfo)
           } else {
             # The characterwise return is deliberately NOT random-expectation
             # normalized for logical `chanceCorrect`: `charInfo` is
             # MutualClusteringInfo() against the whole tree, whereas the
             # analytic `zero` baseline above is per-single-split expected MI, so
             # subtracting it would mix incompatible quantities (and the
             # entropy-weighted variant was abandoned -- see the note below the
             # @return docs). Only the Monte-Carlo path (numeric
             # `chanceCorrect`) offers a same-scale empirical baseline, so we
             # return charInfo scaled
             # by its maximum (hBest-like), as shipped since the original
             # implementation (#205).
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
          
           if (isFALSE(chanceCorrect)) {
           } else {
             one <- sum(hh["hChar", 1, ])
             zero <- if (isTRUE(chanceCorrect)) {
               sum(apply(hh["miRand", , ], 2, max))
             } else {
               randTreeMean
             }
             ret <- (sum(charInfo) - zero) / (one - zero)
             if (is.numeric(chanceCorrect)) {
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
#' @param quality Numeric vector of values between -1 and 1, denoting the
#' quality of observations, where 0 is neutral.
#' @param amount Numeric vector of values between 0 and 1, denoting the relative
#' amount of information.
#' @return `QACol()` returns an RGB hex code for a colour, where lighter colours
#' correspond to entries with a higher `amount`; unsaturated colours denote
#' a neutral `quality`; and red/cyan colours denote low/high `quality`.
#' @examples
#' amount <- runif(80, 0, 1)
#' quality <- runif(80, -1, 1)
#' plot(quality, amount, col = QACol(quality, amount), pch = 15)
#' abline(v = 0)
#' @template MRS
#' @importFrom colorspace hex polarLUV
#' @family utility functions
#' @export
QACol <- function(quality, amount) {
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
QCol <- function(quality, amount) {
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
#' @param \dots Additional parameters to [mtext()].
#' @inheritParams ConcordanceTable
#' @export
QALegend <- function(where = c(0.1, 0.3, 0.1, 0.3), n = 5, Col = QACol,
                     xlab = "Amount \U2192", ylab = "Quality \U2192", ...) {
  oPar <- par(fig = where, new = TRUE, mar = rep(0, 4), xpd = NA)
  on.exit(par(oPar))
  n <- rep(n, length.out = 2)
  nA <- n[[2]]
  nQ <- n[[1]]
  amount <- seq(0, 1, length.out = nA)
  quality <- seq(-1, 1, length.out = nQ)
  mat <- outer(amount, quality,
               Vectorize(function (a, q) Col(q, a)))
  image(x = amount, y = quality, z = matrix(1:prod(n), nA, nQ),
        col = mat, axes = FALSE, xlab = "", ylab = "")
  mtext(xlab, side = 1, line = 1, ...)
  mtext(ylab, side = 2, line = 1, ...)
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
#' @param plot Logical specifying whether to draw the plot.
#' @param marginSize Integer scalar or vector controlling summary margin strips.
#' If a scalar (length 1) and greater than zero, both a left strip and a bottom
#' strip are added, each `marginSize` grid cells wide/tall.
#' If a vector (length > 1), each entry controls one side following the usual
#' `par(mar)` order — `c(bottom, left, top, right)` — where a positive value
#' enables that strip with the given width/height and `NA` or `0` suppresses it.
#' The left and right strips are coloured by the characterwise concordance
#' (weighted mean across edges); the bottom and top strips by the edgewise
#' concordance (weighted mean across characters).
#' One blank cell separates each strip from the main grid.
#' @param paintSize Integer scalar or vector.  Adds a painted strip OUTSIDE any
#'   `marginSize` strip, using hue from [TreeTools::PaintTree()] (edges) and the
#'   [PaintCharacters()] algorithm (characters).  A scalar `> 0` adds a right
#'   strip (characters) and a top strip (edges), each `paintSize` cells wide/tall.
#'   A length-4 vector follows `c(bottom, left, top, right)` like `marginSize`;
#'   `NA` or `0` suppresses that side.  One blank cell separates each paint strip
#'   from the adjacent margin strip (or main grid if no margin exists on that side).
#' @param palette Palette specification passed to [TreeTools::PaintTree()].
#'   Either a character string (`"default"`, `"protanopia"`, `"tritanopia"`) or
#'   a function `function(h, s)`.  Ignored when `paintSize` is zero on all sides.
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
#' paint <- PaintTree(tree)
#' plot(tree, edge.col = paint$edgeCol, tip.col = paint$tipCol, edge.width = 2)
#' nodeIndex <- as.integer(rownames(as.Splits(tree)))
#' nodelabels(seq_along(nodeIndex), nodeIndex, adj = c(2, 1),
#'            frame = "none", bg = NULL)
#' QALegend(where = c(0.1, 0.4, 0.1, 0.3))
#'
#' # View information shared by characters and edges
#' ConcordanceTable(tree, dataset, largeClade = 3, col = 2, lwd = 3,
#'                  marginSize = c(0, 0, 1, 2), paintSize = c(1, 2, 0, 0))
#' axis(1)
#' axis(2)
#'
#' # Visualize dataset
#' image(t(`mode<-`(PhyDatToMatrix(dataset), "numeric")), axes = FALSE,
#'       xlab = "Leaf", ylab = "Character")
#' @importFrom graphics abline image mtext
#' @importFrom grDevices col2rgb convertColor rgb
#' @importFrom TreeTools CladeSizes NTip PaintTree
#' @family split support functions
#' @seealso
#' - [SiteConcordance()]: compute underlying concordance values.
#' @export
ConcordanceTable <- function(tree, dataset, Col = QACol, largeClade = 0,
                             xlab = "Edge", ylab = "Character",
                             chanceCorrect = TRUE, plot = TRUE,
                             marginSize = 0L, paintSize = 0L,
                             palette = "default", ...) {
  cc <- ClusteringConcordance(tree, dataset, return = "all",
                              chanceCorrect = chanceCorrect)
  nodes <- seq_len(dim(cc)[[2]])
  info <- .ConcSlice(cc, "hBest") * .ConcSlice(cc, "n")
  amount <- info / max(info, na.rm = TRUE)
  amount[is.na(amount)] <- 0
  quality <- .ConcSlice(cc, "normalized")
  # Plot points with incalculable quality as black, not transparent.
  amount[is.na(quality)] <- 0
  quality[is.na(quality)] <- 0

  col <- matrix(Col(quality, amount), dim(amount)[[1]], dim(amount)[[2]])

  # Parse marginSize: scalar → bottom + left; vector → c(bottom, left, top, right)
  ms <- as.integer(marginSize)
  if (length(ms) == 1L) {
    ms_bottom <- if (!is.na(ms) && ms > 0L) ms else 0L
    ms_left   <- ms_bottom
    ms_top    <- 0L
    ms_right  <- 0L
  } else {
    ms_bottom <- if (!is.na(ms[1L]) && ms[1L] > 0L) ms[1L] else 0L
    ms_left   <- if (length(ms) >= 2L && !is.na(ms[2L]) && ms[2L] > 0L) ms[2L] else 0L
    ms_top    <- if (length(ms) >= 3L && !is.na(ms[3L]) && ms[3L] > 0L) ms[3L] else 0L
    ms_right  <- if (length(ms) >= 4L && !is.na(ms[4L]) && ms[4L] > 0L) ms[4L] else 0L
  }

  # Parse paintSize: scalar → top + right; vector → c(bottom, left, top, right)
  ps <- as.integer(paintSize)
  if (length(ps) == 1L) {
    ps_top    <- if (!is.na(ps) && ps > 0L) ps else 0L
    ps_right  <- ps_top
    ps_bottom <- 0L
    ps_left   <- 0L
  } else {
    ps_bottom <- if (!is.na(ps[1L]) && ps[1L] > 0L) ps[1L] else 0L
    ps_left   <- if (length(ps) >= 2L && !is.na(ps[2L]) && ps[2L] > 0L) ps[2L] else 0L
    ps_top    <- if (length(ps) >= 3L && !is.na(ps[3L]) && ps[3L] > 0L) ps[3L] else 0L
    ps_right  <- if (length(ps) >= 4L && !is.na(ps[4L]) && ps[4L] > 0L) ps[4L] else 0L
  }

  # Paint is outermost; its width is prepended/appended to the margin offset.
  ps_x_offset <- if (ps_left   > 0L) ps_left   + 1L else 0L
  ps_y_offset <- if (ps_bottom > 0L) ps_bottom + 1L else 0L
  ps_x_suffix <- if (ps_right  > 0L) ps_right  + 1L else 0L
  ps_y_suffix <- if (ps_top    > 0L) ps_top    + 1L else 0L

  x_offset <- ps_x_offset + if (ms_left   > 0L) ms_left   + 1L else 0L
  y_offset <- ps_y_offset + if (ms_bottom > 0L) ms_bottom + 1L else 0L
  x_suffix <- (if (ms_right  > 0L) ms_right  + 1L else 0L) + ps_x_suffix
  y_suffix <- (if (ms_top    > 0L) ms_top    + 1L else 0L) + ps_y_suffix

  if (ms_left > 0L || ms_bottom > 0L || ms_top > 0L || ms_right > 0L ||
      ps_left > 0L || ps_bottom > 0L || ps_top > 0L || ps_right > 0L) {
    n_edges <- dim(cc)[[2]]
    n_chars <- dim(cc)[[3]]

    # Marginal concordance: hBest-weighted average of normalized MI
    hBest_w <- .ConcSlice(cc, "hBest")
    hBest_w[is.na(hBest_w)] <- 0
    # `quality` already has NAs zeroed above

    # Extended layout (x = left→right, y = bottom→top):
    #   x: [paint_left] [blank] [margin_left] [blank] [grid] [blank] [margin_right] [blank] [paint_right]
    #   y: [paint_bottom] [blank] [margin_bottom] [blank] [grid] [blank] [margin_top] [blank] [paint_top]
    nx <- x_offset + n_edges + x_suffix
    ny <- y_offset + n_chars + y_suffix
    ext_col <- matrix("#FFFFFF", nx, ny)

    xi <- (x_offset + 1L):(x_offset + n_edges)  # x indices of main grid
    yi <- (y_offset + 1L):(y_offset + n_chars)   # y indices of main grid
    ext_col[xi, yi] <- col

    if (ms_left > 0L || ms_right > 0L) {
      denom_c <- colSums(hBest_w)
      char_conc <- pmax(-1, pmin(1,
        ifelse(denom_c == 0, 0, colSums(quality * hBest_w) / denom_c)))
      charInfo <- cc["hChar", 1, ] * cc["n", 1, ]
      char_cols <- Col(char_conc, charInfo / max(charInfo))
      if (ms_left > 0L) {
        for (i in seq_len(ms_left)) ext_col[ps_x_offset + i, yi] <- char_cols
      }
      if (ms_right > 0L) {
        for (i in seq_len(ms_right)) ext_col[x_offset + n_edges + 1L + i, yi] <- char_cols
      }
    }
    if (ms_bottom > 0L || ms_top > 0L) {
      denom_e <- rowSums(hBest_w)
      edge_conc <- pmax(-1, pmin(1,
        ifelse(denom_e == 0, 0, rowSums(quality * hBest_w) / denom_e)))
      edge_cols <- Col(edge_conc, rowMeans(.ConcSlice(cc, "hSplit")))
      if (ms_bottom > 0L) {
        for (j in seq_len(ms_bottom)) ext_col[xi, ps_y_offset + j] <- edge_cols
      }
      if (ms_top > 0L) {
        for (j in seq_len(ms_top)) ext_col[xi, y_offset + n_chars + 1L + j] <- edge_cols
      }
    }

    if (ps_left > 0L || ps_right > 0L || ps_top > 0L || ps_bottom > 0L) {
      paint <- PaintTree(tree, palette)
      ctNodes <- as.integer(rownames(info))
      edgeIdx <- match(ctNodes, tree[["edge"]][, 2L])
      edge_paint_cols <- paint$edgeCol[edgeIdx]

      if (ps_left > 0L || ps_right > 0L) {
        # Per-character colours: Lab-weighted mean of edge paint colours,
        # reusing `amount` (= relInfo) and `quality` already NA-zeroed above.
        labMat <- matrix(
          convertColor(t(col2rgb(edge_paint_cols)) / 255, from = "sRGB", to = "Lab"),
          ncol = 3L
        )
        wMat_p   <- pmax(quality, 0) * amount
        wSum_p   <- colSums(wMat_p)
        noInfo_p <- wSum_p == 0
        labAvg_p <- t(t(labMat) %*% wMat_p) / ifelse(noInfo_p, 1, wSum_p)
        rgbAvg_p <- matrix(
          pmax(0, pmin(1, convertColor(labAvg_p, from = "Lab", to = "sRGB"))),
          ncol = 3L
        )
        char_paint_cols <- rgb(rgbAvg_p[, 1L], rgbAvg_p[, 2L], rgbAvg_p[, 3L])
        char_paint_cols[noInfo_p] <- "#888888"

        if (ps_left  > 0L)
          for (i in seq_len(ps_left))  ext_col[i, yi] <- char_paint_cols
        if (ps_right > 0L)
          for (i in seq_len(ps_right)) ext_col[nx - ps_right + i, yi] <- char_paint_cols
      }
      if (ps_bottom > 0L)
        for (j in seq_len(ps_bottom)) ext_col[xi, j] <- edge_paint_cols
      if (ps_top    > 0L)
        for (j in seq_len(ps_top))    ext_col[xi, ny - ps_top + j] <- edge_paint_cols
    }

    image(seq_len(nx), seq_len(ny),
          matrix(seq_len(nx * ny), nx, ny),
          col = as.vector(ext_col),
          frame.plot = FALSE, axes = FALSE,
          xlab = xlab, ylab = ylab)
  } else {
    image(nodes, seq_len(dim(cc)[[3]]),
          matrix(1:prod(dim(amount)), dim(amount)[[1]]),
          frame.plot = FALSE, axes = FALSE,
          col = col, xlab = xlab, ylab = ylab)
  }

  if (largeClade > 1) {
    cladeSize <- CladeSizes(tree)
    edge <- tree[["edge"]]
    parent <- edge[, 1]
    child <- edge[, 2]
    bigNode <- vapply(as.integer(colnames(cc)), function (node) {
      all(cladeSize[child[parent == parent[child == node]]] >= largeClade)
    }, logical(1))
    abline(v = nodes[bigNode] + x_offset - 0.5, ...)
  }
  invisible(list(info = info, relInfo = amount, quality = quality, col = col))
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
#' `NaN` is returned for a character that no split in `tree` is informative
#' about (zero possible information, a 0 / 0 division).
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
#' `QuartetConcordance()` measures the agreement between each character and each
#' split in the currency of quartets (sets of four leaves).
#' A quartet with the characters `0 0 0 1` is not decisive, as
#' all relationships between those leaves are equally parsimonious.
#' But a quartet with characters `0 0 1 1` is decisive, and is concordant
#' with any tree that groups the first two leaves together to the exclusion
#' of the second.
#' In contrast to the site concordance factor
#' \insertCite{Minh2020}{TreeSearch}, `QuartetConcordance()` considers all
#' quartets that are decisive for a branch.
#' Doing so circumvents the criticisms of \insertCite{Goloboff2024;textual}{TreeSearch}.
#'
#' The quartets that a split resolves are not logically independent: because
#' \eqn{(ab, cd)} and \eqn{(ab, ce)} jointly entail \eqn{(ab, de)}
#' \insertCite{Nelson1992}{TreeSearch}, a split of sizes \eqn{(k, t - k)}
#' resolves \eqn{\binom{k}{2}\binom{t - k}{2}} quartets, of which only
#' \eqn{(k - 1)(t - k - 1)} are non-redundant.
#' `unit` selects which of the two is counted.  By default, agreement is
#' measured against the non-redundant content, so that a character scores full
#' marks only where it displays the split, rather than merely agreeing with a
#' large combinatorial volume of its quartets.
#'
#' By default, the reported value weights each site by the quartet content it
#' shares with the split.
#' If `weight = FALSE`, the reported value is the mean of the concordance
#' value for each site.
#' Consider a split associated with two sites (counting in `unit = "quartet"`):
#' one that is concordant with 25&percnt; of 96 decisive quartets, and
#' a second that is concordant with 75&percnt; of 4 decisive quartets.
#' If `weight = TRUE`, the split concordance will be 24 + 3 / 96 + 4 = 27&percnt;.
#' If `weight = FALSE`, the split concordance will be mean(75&percnt;, 25&percnt;) = 50&percnt;.
#'
#' `QuartetConcordance()` is computed exactly, using all quartets, 
#' rather than a random subsample \insertCite{@cf. @Minh2020}{TreeSearch}.
#' Ambiguous and inapplicable tokens are treated as containing no grouping
#' information (i.e. `(02)` or `-` are each treated as `?`).
#'
#' `return` is matched (case-insensitively, and partially) against `"edge"`
#' and `"char"`; any other value raises an error.
#'
#' @return
#' `QuartetConcordance(return = "edge")` returns a numeric vector giving the
#' concordance index at each split across all sites; names specify the number of
#' each corresponding split in `tree`.
#'
#' `QuartetConcordance(return = "char")` returns a numeric vector giving the
#' concordance index calculated at each site, averaged across all splits.
#' Unlike the `"edge"` result, this vector is unnamed: characters retain no
#' persistent identifier once collapsed to patterns.
#'
#' With `unit = "nrqs"`, per-character values are much smaller than 1 even
#' where a character perfectly matches a split in a tree. Values rank characters
#' on a given tree, but are not comparable between trees.
#'
#' @param weight Logical specifying whether to weight sites according to the
#' quartet content that they share with each split.
#' @param unit Character specifying the currency in which quartets are counted:
#'   - `"nrqs"` (default): non-redundant quartet statements.
#'     At a given split, only a character whose own split is identical scores
#'     full marks; nested (compatible) characters receive partial support, and 
#'     incompatible characters score lower still.
#'   - `"quartet"`: each resolved quartet counts once,
#'     analogous to the site concordance factor \insertCite{Minh2020}{TreeSearch}.
#' @references \insertAllCited{}
#' @importFrom ape keep.tip
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom utils combn
#' @importFrom TreeTools as.Splits PhyDatToMatrix TipLabels
#' @export
QuartetConcordance <- function(
  tree,
  dataset = NULL,
  weight = TRUE,
  return = "edge",
  unit = c("nrqs", "quartet"),
  chanceCorrect = TRUE
) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }
  unit <- match.arg(unit)
  if (!isFALSE(chanceCorrect)) {
    if (!isTRUE(chanceCorrect)) {
      if (!is.numeric(chanceCorrect) || length(chanceCorrect) != 1L ||
          is.na(chanceCorrect) || chanceCorrect < 1) {
        stop("`chanceCorrect` must be FALSE, TRUE, or a positive integer.")
      }
    }
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
  
  contrast <- attr(dataset, "contrast")
  charLevels <- attr(dataset, "allLevels")
  
  appCols <- colnames(contrast) != "-"
  isInapp <- charLevels == "-"
  # A level combining an applicable state with "-" (e.g. `{0,-}`) sets
  # exactly one non-"-" column, so it would otherwise pass the `rowSums`
  # check below as if it were a pure single-state (grouping) level.
  combinesInapp <- rowSums(contrast[, !appCols, drop = FALSE] > 0) > 0 & !isInapp
  isAmbig <- rowSums(contrast[, appCols, drop = FALSE]) > 1 | combinesInapp
  isGrouping <- !isAmbig & !isInapp
  
  # For each grouping level, which column of the contrast matrix does it uniquely set?
  groupingCols <- apply(contrast[isGrouping, , drop = FALSE] > 0, 1, which)
  
  levelToInt <- rep(NA_integer_, length(charLevels))
  levelToInt[isGrouping] <- as.integer(groupingCols)
  
  characters <- PhyDatToMatrix(dataset)
  charInt <- array(
    levelToInt[match(characters, charLevels)],
    dim      = dim(characters),
    dimnames = dimnames(characters)
  )
  
  options <- c("edge", "char")
  matched <- pmatch(tolower(trimws(return)), options, nomatch = NA_integer_)
  if (is.na(matched)) {
    stop("`return` must (partially) match one of ",
         paste(sQuote(options), collapse = ", "))
  }
  return <- options[[matched]]

  if (unit == "nrqs") {
    # Return:
    return(.NrqsConcordance(logiSplits, charInt, weight, return, splits,
                            chanceCorrect))
  }

  raw_counts <- quartet_concordance(logiSplits, charInt)

  num <- raw_counts$concordant
  den <- raw_counts$decisive

  # Chance correction (option): re-zero the observed concordant/decisive ratio
  # against the ratio expected under the same fixed-marginal null used for the
  # NRQS currency.  Only `conc` and `dec` vary under the null (the split sizes
  # and state counts are fixed), so we need E[conc] and E[dec] per (split, char),
  # computed exactly from the hypergeometric pmf (no floors here, unlike NRQS)
  # or by Monte-Carlo tip-shuffle, then re-zero the pooled ratio.
  doNorm <- !isFALSE(chanceCorrect)
  if (doNorm) {
    base <- if (isTRUE(chanceCorrect)) {
      .QuartetExpect(charInt, logiSplits)
    } else {
      .QuartetMC(charInt, logiSplits, chanceCorrect)
    }
    eNum <- base[["concordant"]]
    eDen <- base[["decisive"]]
  }

  if (return == "edge") {
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
      if (doNorm) {
        bDen <- rowSums(eDen)
        base <- ifelse(bDen == 0, NA_real_, rowSums(eNum) / bDen)
        ret <- .RezeroGuarded(ret, base)
      }
    } else {
      # Mean of ratios per site
      # Avoid division by zero (0/0 -> NaN -> NA handled by na.rm)
      ratios <- num / den
      # Replace NaN/Inf with NA for rowMeans calculation
      ratios[!is.finite(ratios)] <- NA
      if (doNorm) {
        bRatios <- eNum / eDen
        bRatios[!is.finite(bRatios)] <- NA
        # Average observed and baseline over the same (split, char) cells.
        naMask <- is.na(ratios) | is.na(bRatios)
        ratios[naMask] <- NA
        bRatios[naMask] <- NA
        ret <- .RezeroGuarded(rowMeans(ratios, na.rm = TRUE),
                              rowMeans(bRatios, na.rm = TRUE))
      } else {
        ret <- rowMeans(ratios, na.rm = TRUE)
      }
    }

    setNames(ret, names(splits))
  } else {
    # return = "char"
    if (isTRUE(weight)) {
      if (doNorm) {
        obs <- vapply(seq_len(ncol(num)), function(i) {
          d <- sum(den[, i]); if (d == 0) NA_real_ else sum(num[, i]) / d
        }, double(1))
        b <- vapply(seq_len(ncol(eNum)), function(i) {
          d <- sum(eDen[, i]); if (d == 0) NA_real_ else sum(eNum[, i]) / d
        }, double(1))
        .RezeroGuarded(obs, b)
      } else {
        vapply(
          seq_len(dim(num)[[2]]),
          function(i) {
            weighted.mean(num[, i] / den[, i], den[, i])
          },
          double(1)
        )
      }
    } else {
      if (doNorm) {
        sObs <- ifelse(den > 0, num / den, NA_real_)
        sBase <- ifelse(eDen > 0, eNum / eDen, NA_real_)
        naMask <- is.na(sObs) | is.na(sBase)
        sObs[naMask] <- NA_real_
        sBase[naMask] <- NA_real_
        obs <- vapply(seq_len(ncol(sObs)), function(i) {
          m <- mean(sObs[, i], na.rm = TRUE); if (is.nan(m)) NA_real_ else m
        }, double(1))
        b <- vapply(seq_len(ncol(sBase)), function(i) {
          m <- mean(sBase[, i], na.rm = TRUE); if (is.nan(m)) NA_real_ else m
        }, double(1))
        .RezeroGuarded(obs, b)
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
}

# Expected raw concordant + decisive quartet counts under the fixed-marginal
# null.  The exact expectation (per state-pair, an exact sum over the
# trivariate-hypergeometric pmf) is computed in C++ (`quartet_expect` in
# src/concordance_expect.cpp) for speed; `.QuartetExpect` is a thin R wrapper.
# The `.QuartetMC` Monte-Carlo baseline reshuffles tokens and re-scores through
# the C++ `quartet_concordance` kernel.  (The C++ path is validated bit-for-bit
# against the earlier R implementation; see dev/benchmarks/frac-quart.)
.QuartetExpect <- function(charInt, logiSplits) {
  quartet_expect(logiSplits, charInt)
}

.QuartetMC <- function(charInt, logiSplits, nRelabel) {
  nSplit <- ncol(logiSplits)
  nChar <- ncol(charInt)
  accConc <- accDec <- matrix(0, nSplit, nChar)
  for (s in seq_len(nRelabel)) {
    shuffled <- charInt
    for (ci in seq_len(nChar)) {
      col <- charInt[, ci]
      scored <- !is.na(col)
      shuffled[scored, ci] <- sample(col[scored])
    }
    kc <- quartet_concordance(logiSplits, shuffled)
    accConc <- accConc + kc[["concordant"]]
    accDec <- accDec + kc[["decisive"]]
  }
  list(concordant = accConc / nRelabel, decisive = accDec / nRelabel)
}

# Nelson-Ladiges fractional ("nrqs") currency for QuartetConcordance().
#
# The quartets a bipartition of sizes (k, t - k) resolves are the 4-cycles of
# the complete bipartite graph K_{k, t-k}; the Nelson-Ladiges entailment
# (ab,cd) + (ab,ce) -> (ab,de) is GF(2) cycle addition, so only the cyclomatic
# number (k - 1)(t - k - 1) of them are non-redundant (the NRQS).
#
# A multi-state character is the disjoint union of its state-pairs: a quartet is
# decisive only when two taxa share one state and two share another, so every
# decisive quartet lives in the K_{n_i, n_j} block between two states.  Those
# blocks are edge-disjoint and the entailment never crosses them (it forces the
# shared "c, d, e" taxa into a single state), so NRQS add over state pairs:
#   W_char = sum_{i<j} (n_i - 1)(n_j - 1),   A = sum_{i<j} A_ij
# with a per state-pair weight of 4 / (n_i n_j).  (Verified by GF(2) rank of the
# decisive / concordant quartet sets.)  Each state-pair is therefore an
# independent binary sub-problem; a binary character is the single-pair case.
#
# Scoring uses coverage ("option b"): per pair, concordant NRQS over the
# *reported* unit's own content (Wk for edges, Wc for characters), so that
# A_ij <= min(Wc, Wk) keeps each ratio in [0, 1] and only an identical
# character-split *pair* scores 1.
#
# That per-pair bound does NOT carry over to the value reported for
# `return = "char"`.  The character score is the M-weighted mean of those
# per-pair ratios across every split in the tree, and a character agrees
# exactly with at most one split while merely being *compatible* with (silent
# about) the rest.  Coverage scores silence as failure to cover, so even a
# character identical to one of the tree's own splits averages well below 1:
# 0.49 for a 22|26 split of a 48-leaf tree, 0.25 for a 2|46 split.  The
# attainable maximum therefore depends on the tree, and character scores are
# comparable to one another but not to an absolute ceiling of 1.
#
# There is no cheap repair.  Substituting the character's concordant NRQS with
# respect to the whole tree would require counting a union of per-split
# concordant sets, and those sets overlap heavily: summing A_ij across the 45
# splits of the reference tree used in testing overshoots Wc by ~11x, so the
# whole-tree count is not additive over splits and has no closed form here.
# Taking max_e A_ij / Wc instead does attain 1, but best-match selection is
# upward-biased in the absence of signal (see the `tree` return, which carries
# the same warning) and E[max] is not available from the analytic per-split
# expectations in `nrqs_expect()`, so `chanceCorrect = TRUE` could not support
# it.  Callers who need a per-character score bounded at 1 should use
# `unit = "quartet"`, whose char path reaches 1 for an identical character.
#
# Pairs are pooled by the shared information
# M = min(Wc, Wk) (the hBest analogue), which is symmetric between character and
# split so both `return`s pool by the same amount.  Ambiguous / inapplicable /
# absent tokens drop out per character (treated as "?", as in the quartet path),
# giving each character its own effective taxon count.
.NrqsConcordance <- function(logiSplits, charInt, weight, return, splits,
                             chanceCorrect) {
  nSplit <- ncol(logiSplits)
  nChar <- ncol(charInt)
  pos <- function(z) {
    z[z < 0] <- 0
    z
  }

  # Observed contributions to the pools, summed over each character's
  # state-pairs; the chance baseline (if requested) uses the same accumulators
  # for the expected pools, so the re-zero divides like against like.
  numEdge <- numChar <- denM <- matrix(0, nSplit, nChar)
  wcTot <- numeric(nChar)                  # character NRQS content, split-free
  doNorm <- !isFALSE(chanceCorrect)
  if (doNorm) {
    baseNumEdge <- baseNumChar <- baseDenM <- matrix(0, nSplit, nChar)
    if (isTRUE(chanceCorrect)) {
      # Exact hypergeometric expectation for all characters at once (C++;
      # `nrqs_expect` in src/concordance_expect.cpp).  Uninformative characters
      # yield zero pools, matching the observed loop's `wcTot > 0` guard.
      be <- nrqs_expect(logiSplits, charInt)
      baseNumEdge <- be[["numEdge"]]
      baseNumChar <- be[["numChar"]]
      baseDenM <- be[["denM"]]
    }
  }

  for (ci in seq_len(nChar)) {
    col <- charInt[, ci]
    obs <- .CharNrqsContrib(col, logiSplits, pos)
    numEdge[, ci] <- obs[["numEdge"]]
    numChar[, ci] <- obs[["numChar"]]
    denM[, ci] <- obs[["denM"]]
    wcTot[ci] <- obs[["wcTot"]]
    if (doNorm && !isTRUE(chanceCorrect) && obs[["wcTot"]] > 0) {
      # Monte-Carlo shuffle
      base <- .CharNrqsMC(col, logiSplits, pos, chanceCorrect)
      baseNumEdge[, ci] <- base[["numEdge"]]
      baseNumChar[, ci] <- base[["numChar"]]
      baseDenM[, ci] <- base[["denM"]]
    }
  }

  informative <- wcTot > 0

  if (return == "edge") {
    # edge: one value per split
    if (isTRUE(weight)) {
      denom <- rowSums(denM)
      ret <- ifelse(denom == 0, NA_real_, rowSums(numEdge) / denom)
      if (doNorm) {
        bDen <- rowSums(baseDenM)
        base <- ifelse(bDen == 0, NA_real_, rowSums(baseNumEdge) / bDen)
        ret <- .RezeroGuarded(ret, base)
      }
    } else {
      # Mean per-site quality over informative characters; uninformative
      # characters carry no NRQS and are dropped, as in the quartet path.
      sEdge <- ifelse(denM > 0, numEdge / denM, NA_real_)
      if (doNorm) {
        sBase <- ifelse(baseDenM > 0, baseNumEdge / baseDenM, NA_real_)
        # Average observed and baseline over the SAME (split, char) cells: a
        # degenerate multi-state pair can be dropped from one but not the other
        # (observed wk = 0 yet E[m] > 0), which would otherwise re-zero mismatched
        # populations.
        naMask <- is.na(sEdge) | is.na(sBase)
        sEdge[naMask] <- NA_real_
        sBase[naMask] <- NA_real_
        ret <- .RezeroGuarded(.RowMeanInformative(sEdge, informative, nSplit),
                              .RowMeanInformative(sBase, informative, nSplit))
      } else {
        ret <- .RowMeanInformative(sEdge, informative, nSplit)
      }
    }
    setNames(ret, names(splits))
  } else {
    # char: one value per character
    if (isTRUE(weight)) {
      denom <- colSums(denM)
      ret <- ifelse(denom == 0, NA_real_, colSums(numChar) / denom)
      if (doNorm) {
        bDen <- colSums(baseDenM)
        base <- ifelse(bDen == 0, NA_real_, colSums(baseNumChar) / bDen)
        ret <- .RezeroGuarded(ret, base)
      }
      ret
    } else {
      sChar <- ifelse(denM > 0, numChar / denM, NA_real_)
      if (doNorm) {
        sBase <- ifelse(baseDenM > 0, baseNumChar / baseDenM, NA_real_)
        # Match cells, as in the edge return above.
        naMask <- is.na(sChar) | is.na(sBase)
        sChar[naMask] <- NA_real_
        sBase[naMask] <- NA_real_
        ret <- .RezeroGuarded(.ColMeanInformative(sChar, informative, nChar),
                              .ColMeanInformative(sBase, informative, nChar))
      } else {
        ret <- .ColMeanInformative(sChar, informative, nChar)
      }
      ret
    }
  }
}

# Observed per-character NRQS contributions, summed over the character's
# state-pairs, returned as per-split vectors.  Shared by the observed pass and
# the Monte-Carlo baseline so the two use byte-identical arithmetic.
.CharNrqsContrib <- function(col, logiSplits, pos) {
  nSplit <- ncol(logiSplits)
  numEdge <- numChar <- denM <- numeric(nSplit)
  wcTot <- 0
  scored <- !is.na(col)
  states <- sort(unique(col[scored]))
  nStates <- length(states)
  if (nStates >= 2L) {
    for (a in seq_len(nStates - 1L)) {
      for (b in seq(a + 1L, nStates)) {
        inI <- scored & col == states[a]
        inJ <- scored & col == states[b]
        aI <- colSums(logiSplits & inI)      # state i, side A
        bI <- colSums(!logiSplits & inI)     # state i, side B
        aJ <- colSums(logiSplits & inJ)      # state j, side A
        bJ <- colSums(!logiSplits & inJ)     # state j, side B
        nI <- aI + bI                        # n_i (constant across splits)
        nJ <- aJ + bJ                        # n_j
        mA <- aI + aJ                        # taxa of this pair on side A
        tP <- nI + nJ                        # taxa scored in this pair
        # Concordant NRQS; the (x - 1)_+ floors stop self-agreement exceeding 1.
        aij <- pos(aI - 1) * pos(bJ - 1) + pos(bI - 1) * pos(aJ - 1)
        wc <- pos(nI - 1) * pos(nJ - 1)      # pair's character content
        wk <- pos(mA - 1) * pos(tP - mA - 1) # pair's split content
        m <- pmin(wc, wk)                    # shared information
        denM <- denM + m
        numEdge <- numEdge + ifelse(wk > 0, m * aij / wk, 0)
        numChar <- numChar + ifelse(wc > 0, m * aij / wc, 0)
        wcTot <- wcTot + wc[1]               # wc is constant across splits
      }
    }
  }
  list(numEdge = numEdge, numChar = numChar, denM = denM, wcTot = wcTot)
}

# The exact NRQS expectation E[m], E[m*A/wk], E[m*A/wc] is computed in C++
# (`nrqs_expect` in src/concordance_expect.cpp) for all characters at once, and
# consumed directly by `.NrqsConcordance`; there is no per-character R exact
# helper.  The Monte-Carlo baseline below is the opt-in
# `chanceCorrect = <int>` path.

# Monte-Carlo baseline: average `.CharNrqsContrib` over `nRelabel` random
# reassignments of the character's tokens across its scored leaves (the split,
# the scored set and the state counts are held fixed).  Averaging the pools
# (rather than per-split ratios) mirrors the exact ratio-of-expectations path.
.CharNrqsMC <- function(col, logiSplits, pos, nRelabel) {
  nSplit <- ncol(logiSplits)
  scored <- !is.na(col)
  tokens <- col[scored]
  accEdge <- accChar <- accDen <- numeric(nSplit)
  for (i in seq_len(nRelabel)) {
    shuffled <- col
    shuffled[scored] <- sample(tokens)
    cc <- .CharNrqsContrib(shuffled, logiSplits, pos)
    accEdge <- accEdge + cc[["numEdge"]]
    accChar <- accChar + cc[["numChar"]]
    accDen <- accDen + cc[["denM"]]
  }
  list(numEdge = accEdge / nRelabel,
       numChar = accChar / nRelabel,
       denM = accDen / nRelabel)
}

# Mean of per-split quality across informative characters, matching the
# quartet path's treatment of uninformative characters (dropped).
.RowMeanInformative <- function(mat, informative, nSplit) {
  ret <- if (any(informative)) {
    rowMeans(mat[, informative, drop = FALSE], na.rm = TRUE)
  } else {
    rep(NA_real_, nSplit)
  }
  ret[is.nan(ret)] <- NA_real_
  ret
}

.ColMeanInformative <- function(mat, informative, nChar) {
  vapply(seq_len(nChar), function(ci) {
    if (informative[ci]) {
      m <- mean(mat[, ci], na.rm = TRUE)
      if (is.nan(m)) NA_real_ else m
    } else {
      NA_real_
    }
  }, double(1))
}

# Re-zero to a chance baseline, guarding the degenerate `zero -> 1` case
# (returns NA rather than dividing by ~0).  Values are returned UNCLAMPED, as
# in `ClusteringConcordance`; clamping to [-1, 1] is deferred to plotting.
.RezeroGuarded <- function(value, zero) {
  denom <- 1 - zero
  ifelse(is.na(value) | is.na(zero) | denom < sqrt(.Machine[["double.eps"]]),
         NA_real_, (value - zero) / denom)
}

.ExpectedMICache <- new.env(hash = TRUE, parent = emptyenv())
# Bound on cache size: a session computing concordance for many differently
# -sized trees/characters would otherwise grow this cache without limit.
# Simplest possible bounded policy: wipe the whole cache once full, rather
# than tracking per-entry recency for an LRU/LFU scheme.  The entry count is
# tracked separately (`.ExpectedMICacheSize`) rather than read via `length()`
# on every miss: `length()` on a hashed environment is O(n), which would make
# filling the cache to its bound an O(limit^2) operation.
.ExpectedMICacheLimit <- 100000L
.ExpectedMICacheSize <- new.env(parent = emptyenv())
.ExpectedMICacheSize$n <- 0L

# @param a must be a vector of length <= 2
# @param b may be longer
.ExpectedMI <- function(a, b) {
  if (length(a) < 2 || length(b) < 2) {
    0
  } else {
    key <- mi_key(a, b)
    if (!is.null(.ExpectedMICache[[key]])) {
      .ExpectedMICache[[key]]
    } else {
      ret <- expected_mi(a, b)

      # Cache:
      if (.ExpectedMICacheSize$n >= .ExpectedMICacheLimit) {
        rm(list = ls(.ExpectedMICache, all.names = TRUE),
           envir = .ExpectedMICache)
        .ExpectedMICacheSize$n <- 0L
      }
      .ExpectedMICache[[key]] <- ret
      .ExpectedMICacheSize$n <- .ExpectedMICacheSize$n + 1L
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

# Extract row `name` of a (measure x split x character) concordance array as
# a (split x character) matrix, preserving dimnames.  Plain `[` silently
# drops the split or character axis whenever it has extent one, corrupting
# shape and rownames downstream (e.g. a single-split tree, or a
# single-pattern dataset).
.ConcSlice <- function(x, name) {
  array(x[name, , , drop = FALSE], dim(x)[2:3], dimnames(x)[2:3])
}

#' @rdname SiteConcordance
#' 
#' @details
#' `PhylogeneticConcordance()` treats each character in `dataset` as a
#' phylogenetic hypothesis and measures the extent to which it supports the
#' splits of `tree`. Each character is first interpreted as a tree (or set of
#' trees) in which taxa sharing the same token form a clade. Only splits for
#' which the character contains at least four relevant taxa can contribute
#' information.
#'
#' For each split, the function identifies which characters could potentially
#' support that split (i.e. those for which the induced subtrees contain
#' informative structure), and among these, which characters are actually
#' compatible with the split. The concordance value for each split is the
#' proportion of informative characters that support it.
#' A value of 1 indicates that all characters informative for that subset of
#' taxa support the split; a value of 0 indicates that none do. Characters that
#' contain only ambiguous or uninformative states for the relevant taxa do not
#' affect the result.
#' 
#' @return `PhylogeneticConcordance()` returns a numeric vector giving the
#' phylogenetic information of each split in `tree`, named according to the
#' split's internal numbering.
#' `NaN` is returned for a split that no character is informative about
#' (zero possible information, a 0 / 0 division).
#'
#' @importFrom TreeTools as.multiPhylo CladisticInfo CompatibleSplits
#' @importFrom TreeTools MatchStrings
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
#' `NaN` is returned for a character that no split in `tree` is informative
#' about (zero possible information, a 0 / 0 division).
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
#' not be calculated, and so are not included in the totals above.
#'
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ConcordantInformation(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom TreeTools CharacterInformation Log2UnrootedMult Log2Unrooted
#' @importFrom TreeTools MatchStrings
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
            paste0(which(index %in% which(na)), collapse = ", "),
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
