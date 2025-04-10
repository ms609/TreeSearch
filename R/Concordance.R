#' Calculate site concordance factor
#' 
#' The site concordance factor is a measure 
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
#' The `methods` argument selects from different methods for quantifying the
#' concordance associated with a given edge of a tree.
#' 
#' If `method = "split"` (the default), the concordance is calculated by
#' counting all quartets that are decisive for a split, i.e. for the split
#' $A|B$, any quartet that contains two leaves from set $A$ and two from set $B$.
#' This value can be interpreted as the proportion of all decisive quartets
#' that are concordant with a split.
#' If `method = "sitemean"`, the reported value is the mean of the concordance
#' value for each site.
#' Consider a split associated with two sites:
#' one that is concordant with 25% of 96 decisive quartets, and
#' a second that is concordant with 75% of 4 decisive quartets.
#' `method = "split"` returns a value of 24 + 3 / 96 + 4 = 27%.
#' `method = "sitemean"` returns mean(75%, 25%) = 50%.
#' 
#' `method = "minh"` uses the approach of
#' \insertCite{Minh2020;textual}{TreeSearch} to compute the site concordance
#' factor.
#' Briefly, this interprets each edge as defining *four* clades, and counts
#' the status of quartets that contain exactly one leaf from each of these
#' clades.  This is a subset of the quartets considered by other methods.
#' `method = "iqtree"` uses the \insertCite{Minh2020;textual}{TreeSearch}
#' approach as [implemented in IQ-TREE](
#' https://github.com/iqtree/iqtree2/issues/415).
#' 
#' 
#' `QuartetConcordance()` is computed exactly, using all quartets.
#' Other implementations (e.g. IQ-TREE) follow
#' \insertCite{@Minh2020;textual}{TreeSearch} in using a random subsample
#' of quartets for a faster computation, potentially at the expense of accuracy.
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
#' For an overview of the use and interpretation of concordance factors,
#' see \insertCite{Lanfear2024}{TreeSearch}.
#' 
#' **NOTE:** These functions are under development. They are incompletely
#' tested, and may change without notice.
#' Complete documentation and discussion will follow in due course.
#' 
# # Renumber before MaximizeParsimony, for `tree`
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @param method Character vector specifying which concordance measures to
#' calculate.  See details section for available options.
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
#' @importFrom TreeTools as.Splits DescendantTips PhyDatToMatrix TipLabels
#' @name SiteConcordance
#' @family split support functions
#' @export
# TODO support `method = c("split", "sitemean", "minh")`
QuartetConcordance <- function (tree, dataset = NULL, method = "split",
                                n = 250) {
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
  characters <- PhyDatToMatrix(dataset, ambigNA = TRUE)
  
  if (method %in% c("minh", "minhq", "iqtree")) {
    quarters <- .TreeQuarters(tree)
    
    cli_progress_bar(name = "Quartet concordance", total = dim(quarters)[[2]])
    on.exit(cli_progress_done(.envir = parent.frame(2)))
    if (method == "minh") {
      sCF <- apply(quarters, 2, function(abcd) {
        resample <- function(x) x[sample.int(length(x), 1)]
        leaves <- replicate(n, {
          vapply(0:3, function(q) resample(which(abcd == q)), integer(1))
        })
        cfQ <- apply(leaves, 2, function(q) {
          qCh <- `mode<-`(characters[q, , drop = FALSE], "numeric")
          concordant <- sum(qCh[1, ] == qCh[2, ] &
            qCh[2, ] != qCh[3, ] &
            qCh[3, ] == qCh[4, ])
          decisive <- dim(qCh)[[2]]
          decisive <- sum(colSums(apply(qCh + 1, 2, tabulate, 4) >= 2) > 1)
          decisive <- sum(colSums(apply(qCh + 1, 2, tabulate, 4) >= 2) > 0)
          concordant / decisive
        })
        mean(cfQ, na.rm = TRUE)
      })
    } else {
      concDeci <- apply(quarters, 2, function (q) {
        #cli_progress_update(1, .envir = parent.frame(2))
        quarts <- apply(characters, 2, function (char) {
          tab <- table(q, char)
          nCol <- dim(tab)[[2]]
          if (nCol > 1L) {
            # RULES:
            # We must pick one leaf per row
            # Only parsimony informative are decisive
            # To be concordant, leaves 1 and 2 must come from the same column
            
            # Example table:
            #   q  0   1   2
            #   A  2   3   0
            #   B  2   3   2
            #   C  0   1   2
            #   D  2   1   2
            # 
            
            concordant <- sum(apply(combn(seq_len(nCol), 2), 2, function(ij) {
              i <- ij[[1]]
              j <- ij[[2]]
              # Taxa from A and B from column i, C and D from j
              prod(tab[1:2, i], tab[3:4, j]) +
                # Taxa from A and B from column j, C and D from i
                prod(tab[1:2, j], tab[3:4, i])
            }))
            
            if (method == "minhq") {
              # To be parsimony informative, we must pick two leaves from each of
              #   two columns
              
              
              decisive <- sum(apply(combn(seq_len(nCol), 2), 2, function(ij) {
                # With three columns, we'll encounter i = 1, 1, 2; j = 2, 3, 3
                i <- ij[[1]]
                j <- ij[[2]]
                sum(apply(combn(4, 2), 2, function(kl) {
                  # We can pick taxa AB, AC, AD, BC, BD, CD from column i,
                  # and the other pair from col j
                  pair1 <- kl
                  pair2 <- (1:4)[-kl]
                  prod(tab[pair1, i], tab[pair2, j])
                }))
              }))
            } else {
              # IQ-TREE seems to use "variant" in place of "parsimony informative"
              
              nPatterns <- prod(rowSums(tab))
              allSame <- sum(apply(tab, 2, prod))
              
              decisive <- nPatterns# - allSame
              
            }
            c(concordant, decisive)
          } else {
            c(0L, 0L)
          }
        })
        rowSums(quarts)
      })
      concDeci[1, ] / concDeci[2, ]
    }
  } else {
    splits <- as.Splits(tree, dataset)
    logiSplits <- vapply(seq_along(splits), function (i) as.logical(splits[[i]]),
                         logical(NTip(dataset)))
    cli_progress_bar(name = "Quartet concordance", total = dim(logiSplits)[[2]])
    on.exit(cli_progress_done(.envir = parent.frame(2)))
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
      if (method == "split") {
        quartSums <- rowSums(quarts)
        ifelse(is.nan(quartSums[[2]]), NA_real_, quartSums[[1]] / quartSums[[2]])
      } else if (method == "sitemean") {
        mean(ifelse(is.nan(quarts[2, ]), NA_real_, quarts[1, ] / quarts[2, ]),
             na.rm = TRUE)
      } else {
        stop("Unrecognized method: ", method)
      }
    }), names(splits))
  }
}

#' @importFrom TreeDist Entropy
.Entropy <- function (...) {
  Entropy(c(...) / sum(...))
}

#' @rdname SiteConcordance
#' @importFrom TreeTools Subsplit
#' @importFrom stats setNames
#' @export
ClusteringConcordance <- function (tree, dataset) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  dataset <- dataset[TipLabels(tree)]
  splits <- as.logical(as.Splits(tree))
  
  at <- attributes(dataset)
  cont <- at[["contrast"]]
  if ("-" %in% colnames(cont)) {
    cont[cont[, "-"] > 0, ] <- 1
  }
  ambiguous <- rowSums(cont) != 1
  
  mat <- matrix(unlist(dataset), length(dataset), byrow = TRUE)
  mat[mat %in% which(ambiguous)] <- NA
  mat <- apply(mat, 2, function (x) {
    uniques <- table(x) == 1
    x[x %in% names(uniques[uniques])] <- NA
    x
  })
  
  h <- apply(mat, 2, function (char) {
    aChar <- !is.na(char)
    ch <- char[aChar]
    hChar <- .Entropy(table(ch))
    h <- apply(splits[, aChar], 1, function (spl) {
      c(hSpl = .Entropy(table(spl)), hJoint =  .Entropy(table(ch, spl)))
    })
    
    cbind(hSum = hChar + h["hSpl", ], joint = h["hJoint", ])
  })
  
  splitI <- seq_len(dim(splits)[1])
  both <- rowSums(h[splitI, at[["index"]]])
  joint <- rowSums(h[-splitI, at[["index"]]])
  mi <- both - joint
  
  # Return:
  setNames(mi / joint, rownames(splits))
}

#' @rdname SiteConcordance
#' @importFrom TreeTools as.multiPhylo CladisticInfo CompatibleSplits
#' @export
PhylogeneticConcordance <- function (tree, dataset) {
  if (is.null(dataset)) {
    warning("Cannot calculate concordance without `dataset`.")
    return(NULL)
  }
  dataset <- dataset[TipLabels(tree)]
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

.TreeQuarters <- function(tree) {
  if (attr(tree, "order") != "preorder") {
    stop("Tree must be in preorder; try `tree <- Preorder(tree)`")
  }
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]
  nTip <- NTip(tree)
  nodes <- min(parent):max(parent)
  parentEdge <- match(nodes, child)
  descended <- DescendantTips(parent, child)
  childEdges <- lapply(nodes, function(x) which(parent == x))
  daughters <- lapply(childEdges, function(e) child[e])
  
  rootNode <- nTip + which(is.na(parentEdge))
  rootChildren <- child[parent == rootNode]
  nonRootChildren <- setdiff(nodes, c(rootNode, rootChildren))
  
  # Return:
  cbind(
    if (all(rootChildren > nTip)) {
      ret <- integer(nTip)
      groups <- childEdges[[rootChildren - nTip]]
      for (i in 0:3) {
        ret[descended[groups[[i + 1]], ]] <- i
      }
      `dimnames<-`(cbind(ret), list(TipLabels(tree), rootChildren[[1]]))
    },
    `dimnames<-`(vapply(nonRootChildren, function(n) {
      ret <- integer(nTip)
      if (parent[[parentEdge[[n - nTip]]]] == rootNode) {
        stop("Ought something to be here?")
      }
      ret[apply(
        descended[childEdges[[parent[parentEdge[n - nTip]] - nTip]], ],
        2,
        any)] <- 1L
      ce <- childEdges[[n - nTip]]
      # Overprint with own children
      ret[descended[ce[[1]], ]] <- 2L
      ret[descended[ce[[2]], ]] <- 3L
      ret
    }, integer(nTip)), list(TipLabels(tree), nonRootChildren))
  )
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
  dataset <- dataset[TipLabels(tree)]
  splits <- as.multiPhylo(as.Splits(tree))
  characters <- as.multiPhylo(dataset)
  
  support <- rowSums(vapply(characters, function (char) {
    trimmed <- lapply(splits, keep.tip, TipLabels(char))
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
  dataset <- dataset[TipLabels(tree)]
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
  dataset <- dataset[TipLabels(tree)]
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
