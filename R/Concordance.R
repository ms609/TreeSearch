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
#TODO
#' `QuartetConcordance()` in principle supports ambiguous character states,
#' but this has not yet been tested.
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
#' tree <- TreeSearch::referenceTree
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
QuartetConcordance <- function(tree, dataset = NULL, weight = TRUE,
                               return = "mean") {
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
  charLevels[rowSums(attr(dataset, "contrast")) > 1] <- NA
  charLevels <- setdiff(charLevels, "-")
  
  charInt <- matrix(match(characters, charLevels),
                    nrow = nrow(characters),
                    dimnames = dimnames(characters))
  raw_counts <- quartet_concordance(logiSplits, charInt)
  
  num <- raw_counts$concordant
  den <- raw_counts$decisive
  
  options <- c("char", "default")
  return <- options[[pmatch(tolower(trimws(return)), options,
                            nomatch = length(options))]]
  
  if (return == "default") {
    if (isTRUE(weight)) {
      # Sum numerator and denominator across sites (columns), then divide
      # This matches weighted.mean(num/den, den) == sum(num) / sum(den)
      split_sums_num <- rowSums(num)
      split_sums_den <- rowSums(den)
      ret <- ifelse(split_sums_den == 0, NA_real_, split_sums_num / split_sums_den)
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
    p <- num / den
    if (isTRUE(weight)) {
      vapply(seq_len(dim(num)[[2]]), function(i) {
        weighted.mean(num[, i] / den[, i], den[, i])
      }, double(1))
    } else {
      vapply(seq_len(dim(num)[[2]]), function(i) {
        mean(num[den[, i] > 0, i] / den[den[, i] > 0, i])
      }, double(1))
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
  both <- rowSums(h[splitI, at[["index"]], drop = FALSE])
  joint <- rowSums(h[-splitI, at[["index"]], drop = FALSE])
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
