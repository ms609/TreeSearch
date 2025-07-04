#' Calculate site concordance factor
#' 
#' The site concordance factor \insertCite{Minh2020}{TreeSearch} is a measure 
#' of the strength of support that the dataset presents for a given split in a
#' tree.
#' 
#' `QuartetConcordance()` is the proportion of quartets (sets of four leaves)
#' that are decisive for a split which are also concordant with it.
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
#' 
# `ClusteringConcordance()` and `PhylogeneticConcordance()` respectively report
# the proportion of clustering information and phylogenetic information 
# \insertCite{@as defined in @Vinh2010, @SmithDist}{TreeDist} within a dataset
# that is reflected in each split.
# These give smaller values because a split may be compatible with a character
# without being identical to it.
#TODO More thought / explanation needed.
#' 
#TODO Finally, `ProfileConcordance()` (to follow)
#' 
#' **NOTE:** These functions are under development. They are incompletely
#' tested, and may change without notice.
#' Complete documentation and discussion will follow in due course.
#' 
# # Renumber before MaximizeParsimony, for `tree`
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @param weight Logical specifying whether to weight sites according to the
#' number of quartets they are decisive for.
#' 
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @examples 
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]][, 1:20]
#' tree <- referenceTree
#' qc <- QuartetConcordance(tree, dataset)
#' cc <- ClusteringConcordance(tree, dataset)
#' pc <- PhylogeneticConcordance(tree, dataset)
#' spc <- SharedPhylogeneticConcordance(tree, dataset)
#' mcc <- MutualClusteringConcordance(tree, dataset)
#' 
#' oPar <- par(mar = rep(0, 4), cex = 0.8) # Set plotting parameters
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(qc, 3), cex = 0.8)
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(cc, 3), cex = 0.8)
#' par(oPar) # Restore plotting parameters
#' 
#' # Write concordance factors to file
#' labels <- paste0(qc, "/", cc, "/", pc) # "/" is a valid delimiter
#' # Identify the node that corresponds to each label
#' whichNode <- match(TreeTools::NTip(tree) + 1:tree$Nnode, names(qc))
#' 
#' # The contents of tree$node.label will be written at each node
#' tree$node.label <- labels[whichNode]
#' 
#' ape::write.tree(tree) # or write.nexus(tree, file = "mytree.nex")
#' 
#' # Display correlation between concordance factors
#' pairs(cbind(qc, cc, pc, spc, mcc), asp = 1)
#' @template MRS
#' @importFrom ape keep.tip
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom utils combn
#' @importFrom TreeTools as.Splits PhyDatToMatrix TipLabels
#' @name SiteConcordance
#' @family split support functions
#' @export
QuartetConcordance <- function (tree, dataset = NULL, weight = TRUE) {
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
  
  cli_progress_bar(name = "Quartet concordance", total = dim(logiSplits)[[2]])
  setNames(apply(logiSplits, 2, function (split) {
    cli_progress_update(1, .envir = parent.frame(2))
    quarts <- apply(characters, 2, function (char) {
      tab <- table(split, char)
      nCol <- dim(tab)[[2]]
      if (nCol > 1L) {
        # Consider the case
        # split  0   1   2
        # FALSE  2   3   0
        #  TRUE  0   4   2
        # 
        # Concordant quartets with bin i = 1 (character 0) and bin j = 2
        # (character 1) include any combination of the two taxa from 
        # [FALSE, 0], and the two taxa chosen from the four in [TRUE, 1],
        # i.e. choose(2, 2) * choose(4, 2)
        concordant <- sum(vapply(seq_len(nCol), function (i) {
          inBinI <- tab[1, i]
          iChoices <- choose(inBinI, 2)
          sum(vapply(seq_len(nCol)[-i], function (j) {
            inBinJ <- tab[2, j]
            iChoices * choose(inBinJ, 2)
            }, 1))
        }, 1))
        
        # To be discordant, we must select a pair of taxa from TT and from FF;
        # and the character states must group each T with an F
        # e.g. T0 T1  F0 F1
        # T0 T1 F0 F2 would not be discordant - just uninformative
        discordant <- sum(apply(combn(nCol, 2), 2, function (ij) prod(tab[, ij])))
        
        # Only quartets that include two T and two F can be decisive
        # Quartets must also include two pairs of characters
        decisive <- concordant + discordant
        
        # Return the numerator and denominatory of equation 2 in
        # Minh et al. 2020
        c(concordant, decisive)
      } else {
        c(0L, 0L)
      }
    })
    if (isTRUE(weight)) {
      quartSums <- rowSums(quarts)
      ifelse(is.nan(quartSums[[2]]), NA_real_, quartSums[[1]] / quartSums[[2]])
    } else {
      mean(ifelse(is.nan(quarts[2, ]), NA_real_, quarts[1, ] / quarts[2, ]),
           na.rm = TRUE)
    }
  }), names(splits))
}

#' @importFrom TreeDist Entropy
.Entropy <- function(...) {
  Entropy(c(...) / sum(...))
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
#' @param return Character specifying what to return.
#' - `"mean"` returns the mean concordance index at each split across all sites.
#' - `"all"` returns all values calculated during the working for each site at
#'  each split.
#' @param normalize Logical; if `TRUE` the mutual information will be
#' normalized such that zero corresponds to the expected mutual information of
#' a randomly drawn character with the same distribution of tokens.
#' If `FALSE`, zero will correspond to zero mutual information, 
#' even if this is not possible to accomplish in practice.
#' @returns 
#' `ClusteringConcordance(return = "all")` returns a 3D array where each
#' slice corresponds to a site; each row to a split; and each row to a
#' measure of information: `normalized` gives the mutual information (`mi`)
#' normalized such that a value of one corresponds to `hBest`, 
#' which is the lower of `hSplit`, the clustering information
#' (entropy) of the split, and `hChar`, the clustering information of the
#' site / character; and zero corresponds to `miRand`, the expected mutual
#' information of a randomly drawn character with the same distribution of
#' tokens. Negative values denote that the observed tokens contain less mutual
#' information than a random draw.
#' `NA` is returned where $hBest = 0$.
#' `hJoint` gives the joint entropy – the entropy of the
#' confusion matrix of the split and character considered together.
#' 
#' `ClusteringConcordance(return = "mean")` returns a matrix or vector listing
#' for each site the proportion of clustering information across all sites that
#' held in common with the split.
#' 
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ClusteringConcordance(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom abind abind
#' @importFrom pbapply pbapply
#' @importFrom stats setNames
#' @importFrom TreeTools Subsplit
#' @export
ClusteringConcordance <- function (tree, dataset, return = "mean",
                                   normalize = TRUE) {
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
  splits <- as.logical(as.Splits(tree))
  
  at <- attributes(dataset)
  cont <- at[["contrast"]]
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1
  
  mat <- matrix(unlist(dataset), length(dataset), byrow = TRUE)
  mat[mat %in% which(ambiguous)] <- NA_real_
  maxToken <- max(mat, na.rm = TRUE)
  tokens <- as.character(seq_len(maxToken))
  mat <- apply(mat, 2, function (x) {
    uniques <- tabulate(x, maxToken) == 1
    x[x %in% tokens[uniques]] <- NA_real_
    x
  })
  
  h <- simplify2array(apply(mat, 2, function(char) {
    aChar <- !is.na(char)
    ch <- char[aChar]
    chMax <- max(1, ch)
    chTable <- tabulate(ch, chMax)
    n <- length(ch)
    hChar <- Entropy(chTable / n)
    hh <- apply(splits[, aChar, drop = FALSE], 1, function (spl) {
      spTable <- tabulate(spl + 1, 2)
      if (any(spTable < 2)) {
        c(hSplit = 0,
          hJoint = hChar,
          miRand = 0)
      } else {
        c(hSplit = Entropy(spTable / n),
          hJoint = Entropy(tabulate(ch + (spl * chMax), chMax + chMax) / n),
          miRand = .ExpectedMI(spTable, chTable))
      }
    })
    
    rbind(hChar = hChar, hh)
  }, simplify = FALSE))
  
  if (length(dim(h)) == 2) {
    dim(h) <- c(dim(h), 1)
  }
  
  hh <- h[, , at[["index"]], drop = FALSE]
  hBest <- `rownames<-`(pmin(hh["hChar", , , drop = FALSE],
                             hh["hSplit", , , drop = FALSE]), NULL)
  mi <- `rownames<-`(hh["hChar", , , drop = FALSE] + 
                       hh["hSplit", , , drop = FALSE] -
                       hh["hJoint", , , drop = FALSE], NULL)
  miRand <- `rownames<-`(hh["miRand", , , drop = FALSE], NULL)
  
  # Return:
  switch(pmatch(tolower(return), c("all", "mean"), nomatch = 1L),
         # all
         abind(
           along = 1,
           normalized = ifelse(hBest == 0, NA,
                               .Rezero(mi / hBest, miRand / hBest)),
           hh,
           hBest = hBest,
           mi = mi,
           miRand = miRand
         ), {
           # mean
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
         }
  )
}


if (packageVersion("TreeTools") < "1.14.0.9000") {
  MatchStrings <- function(x, table, Fail = stop, max.distance = 0.5, ...) {
    matches <- match(x, table)
    missing <- is.na(matches)
    if (any(missing)) {
      nearMiss <- unlist(lapply(x[missing], agrep, table,
                                max.distance = max.distance, ...),
                         use.names = FALSE, recursive = FALSE)
      message <- paste0("Could not find '", paste(x[missing], collapse = "', '"), 
                        "' in ", deparse(substitute(table)), ".  ",
                        if (length(nearMiss)) {
                          paste0("Did you mean '", 
                                 paste(table[nearMiss], collapse = "', '"), "'?")
                        })
      Fail(message)
    }
    table[matches[!missing]]
  }
} else {
  MatchStrings <- TreeTools::MatchStrings
}

#' Generate colour to depict the amount and quality of observations
#' @param amount Numeric vector of values between 0 and 1, denoting the relative
#' amount of information
#' @param quality Numeric vector of values between -1 and 1, denoting the 
#' quality of observations, where 0 is neutral.
#' @param lMax Maximum lightness of the colour, between 0 and 100.
#' @return `QACol()` returns a colour in HCL space, where darker colours
#' correspond to entries with a higher `amount`; unsaturated colours denote
#' a neutral `quality`; and redder/bluer colours denote low or high `quality`.
#' @examples
#' amount <- runif(80, 0, 1)
#' quality <- runif(80, -1, 1)
#' plot(amount, quality, col = QACol(amount, quality), pch = 15)
#' abline(h = 0)
#' @template MRS
#' @export
QACol <- function(amount, quality, lMax = 70) {
  hcl(
    h = 80 + (quality * 140),
    c = abs(quality) * 100,
    l = (100 - lMax) + ((1 - amount) * lMax)
  )
}

#' @rdname QACol
#' @param where Location of legend, passed to `par(fig = where)`
#' @param n Integer vector giving number of cells to plot in swatch for 
#' `quality` and `amount`.
#' @inheritParams ConcordanceTable
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
  image(x = amount, y = quality,
        z = matrix(1:prod(n), nA, nQ),
        col = mat, axes = FALSE)
  mtext("Amount", side = 1, line = 1)
  mtext("Quality", side = 2, line = 1)
}

#' Plot concordance table
#' 
#' @inheritParams ClusteringConcordance
#' @param Col Function that takes vectors `amount` and `quality` and returns
#' a vector of colours.
#' @param largeClade Integer; if greater than 1, vertical lines will be drawn
#' at edges whose descendants are both contain more than `largeClade` leaves.
#' @param xlab Character giving a label for the x axis.
#' @param ylab Character giving a label for the y axis.
#' @param \dots Arguments to `abline`, to control the appearance of vertical
#' lines marking important edges.
#' @examples
#' # Load data and tree
#' data("congreveLamsdellMatrices", package = "TreeSearch")
#' dataset <- congreveLamsdellMatrices[[1]][, 1:20]
#' tree <- referenceTree
#' 
#' # Plot tree and identify nodes
#' plot(tree)
#' nodeID <- seq_len(tree$Nnode - 1)
#' nodelabels(nodeID, NTip(tree) + 1 + nodeID, adj = c(2, 1),
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
#' @importFrom TreeTools CladeSizes NTip
#' @export
ConcordanceTable <- function(tree, dataset, Col = QACol, largeClade = 0,
                             xlab = "Edge", ylab = "Character", ...) {
  cc <- ClusteringConcordance(tree, dataset, return = "all")
  nodes <- seq_len(dim(cc)[[2]] - 1) # Omit root
  amount <- cc["hBest", -1, ] / max(cc["hBest", -1, ], na.rm = TRUE)
  amount[is.na(amount)] <- 0
  quality <- cc["normalized", -1, ]
  
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
    bigNode <- vapply(nodes + NTip(tree) + 1, function (node) {
      all(cladeSize[child[parent == parent[child == node]]] >= largeClade)
    }, logical(1))
    abline(v = nodes[bigNode] - 0.5, ...)
  }
}

#' @rdname SiteConcordance
#' @importFrom TreeTools as.multiPhylo CladisticInfo CompatibleSplits
#' @export
PhylogeneticConcordance <- function (tree, dataset) {
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
# Mutual clustering information of each split with the split implied by each character
#' @importFrom TreeDist ClusteringEntropy MutualClusteringInfo
#' @export
MutualClusteringConcordance <- function (tree, dataset) {
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
  
  # Return:
  support[, 1] / support[, 2]
}

#' @rdname SiteConcordance
#' @importFrom TreeTools as.multiPhylo
#' @importFrom TreeDist ClusteringInfo SharedPhylogeneticInfo
#' @export
SharedPhylogeneticConcordance <- function (tree, dataset) {
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
  
  # Return:
  support[, 1] / support[, 2]
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
#' @importFrom TreeTools Log2UnrootedMult Log2Unrooted
#' @export
ConcordantInformation <- function (tree, dataset) {
  dataset <- dataset[MatchStrings(TipLabels(tree), names(dataset))]
  originalInfo <- sum(apply(PhyDatToMatrix(dataset), 2, CharacterInformation))
  dataset <- PrepareDataProfile(dataset)
  
  extraSteps <- CharacterLength(tree, dataset, compress = TRUE) -
    MinimumLength(dataset, compress = TRUE)
  chars <- matrix(unlist(dataset), attr(dataset, "nr"))
  ambiguousToken <- which(attr(dataset, "allLevels") == "?")
  asSplits <- apply(chars, 1, function (x) {
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
  c(informationContent = totalInfo,
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
Evaluate <- function (tree, dataset) {
  .Deprecated("ConcordantInformation()")
  ConcordantInformation(tree, dataset)
}

#' @rdname ConcordantInformation
#' @export
ConcordantInfo <- ConcordantInformation
