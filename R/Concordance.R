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
#' NOTE: These functions are under development, and may be incompletely tested
#' or change without notice.
#' Complete documentation and discussion will follow soon.
#' 
#' @template treeParam
#' @template datasetParam
#' 
#' 
#' 
#' @references 
#' \insertAllCited{}
#' 
#' @examples 
#' data('congreveLamsdellMatrices', package = 'TreeSearch')
#' dataset <- congreveLamsdellMatrices[[1]][, 1:20]
#' tree <- referenceTree
#' qc <- QuartetConcordance(tree, dataset)
#' cc <- ClusteringConcordance(tree, dataset)
#' pc <- PhylogeneticConcordance(tree, dataset)
#' spc <- SharedPhylogeneticConcordance(tree, dataset)
#' mcc <- MutualClusteringConcordance(tree, dataset)
#' 
#' oPar <- par(mar = rep(0, 4), cex = 0.8)
#' plot(tree)
#' TreeTools::LabelSplits(tree, signif(qc, 3))
#' TreeTools::LabelSplits(tree, signif(cc, 3))
#' TreeTools::LabelSplits(tree, signif(pc, 3))
#' par(oPar)
#' 
#' pairs(cbind(qc, cc, pc, spc, mcc))
#' @template MRS
#' @importFrom ape keep.tip
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom utils combn
#' @importFrom TreeTools as.Splits PhyDatToMatrix TipLabels
#' @name SiteConcordance
#' @family split support functions
#' @export
QuartetConcordance <- function (tree, dataset) {
  dataset <- dataset[tree$tip.label]
  splits <- as.Splits(tree, dataset)
  logiSplits <- vapply(seq_along(splits), function (i) as.logical(splits[[i]]),
                       logical(NTip(dataset)))
  characters <- PhyDatToMatrix(dataset)
  
  cli_progress_bar(name = 'Quartet concordance', total = dim(logiSplits)[2])
  setNames(apply(logiSplits, 2, function (split) {
    cli_progress_update(1, .envir = parent.frame(2))
    quarts <- rowSums(apply(characters, 2, function (char) {
      tab <- table(split, char)
      nCol <- dim(tab)[2]
      concordant <- sum(vapply(seq_len(nCol), function (i) {
        inBinI <- tab[1, i]
        iChoices <- choose(inBinI, 2)
        sum(vapply(seq_len(nCol)[-i], function (j) {
          inBinJ <- tab[2, j]
          iChoices * choose(inBinJ, 2)
          }, 1))
      }, 1))
      discordant <- sum(apply(combn(nCol, 2), 2, function (ij) prod(tab[, ij])))
      decisive <- concordant + discordant
      c(concordant, decisive)
    }))
    quarts[1] / quarts[2]
  }), names(splits))
}

#' @importFrom TreeDist Entropy
.Entropy <- function (...) {
  Entropy (c(...) / sum(...))
}

#' @rdname SiteConcordance
#' @importFrom TreeTools Subsplit
#' @importFrom stats setNames
#' @export
ClusteringConcordance <- function (tree, dataset) {
  dataset <- dataset[tree$tip.label]
  splits <- as.logical(as.Splits(tree))
  
  at <- attributes(dataset)
  cont <- at$contrast
  if ('-' %in% colnames(cont)) {
    cont[cont[, '-'] > 0, ] <- 1
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
    
    cbind(hSum = hChar + h['hSpl', ], joint = h['hJoint', ])
  })
  
  splitI <- seq_len(dim(splits)[1])
  both <- rowSums(h[splitI, at$index])
  joint <- rowSums(h[-splitI, at$index])
  mi <- both - joint
  
  # Return:
  setNames(mi / joint, rownames(splits))
}

#' @rdname SiteConcordance
#' @importFrom TreeTools as.multiPhylo CladisticInfo CompatibleSplits
#' @export
PhylogeneticConcordance <- function (tree, dataset) {
  dataset <- dataset[tree$tip.label]
  splits <- as.Splits(tree)
  characters <- as.multiPhylo(dataset)
  
  blankRet <- matrix(0, length(splits), 2,
                     dimnames = list(names(splits),
                                     c('concordant', 'possible')))
  
  support <- rowSums(vapply(characters, function (char) {
    ret <- blankRet
    if (NTip(char) > 3L) {
      thinned <- Subsplit(splits, TipLabels(char))
      compatible <- CompatibleSplits(thinned, char)
      if (length(compatible)) {
        ci <- CladisticInfo(thinned)
        ret[names(thinned), 'concordant'] <- ci * apply(compatible, 1, all)
        ret[names(thinned), 'possible'] <- ci
      }
    }
    # Return:
    ret
  }, blankRet), dims = 2)
  
  # Return:
  support[, 1] / support[, 2]
}

#' @rdname SiteConcordance
#' @importFrom TreeDist ClusteringEntropy MutualClusteringInfo
#' @export
MutualClusteringConcordance <- function (tree, dataset) {
  dataset <- dataset[tree$tip.label]
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
  dataset <- dataset[tree$tip.label]
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
#' @template treeParam
#' @template datasetParam
#' @examples
#' data(congreveLamsdellMatrices)
#' myMatrix <- congreveLamsdellMatrices[[10]]
#' ConcordantInformation(TreeTools::NJTree(myMatrix), myMatrix)
#' @template MRS
#' @importFrom TreeTools Log2UnrootedMult Log2Unrooted
#' @export
ConcordantInformation <- function (tree, dataset) {
  dataset <- dataset[tree$tip.label]
  originalInfo <- sum(apply(PhyDatToMatrix(dataset), 2, CharacterInformation))
  dataset <- PrepareDataProfile(dataset)
  
  extraSteps <- CharacterLength(tree, dataset, compress = TRUE) -
    MinimumLength(dataset, compress = TRUE)
  chars <- matrix(unlist(dataset), attr(dataset, 'nr'))
  ambiguousToken <- which(attr(dataset, 'allLevels') == "?")
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
  noise[noise < sqrt(.Machine$double.eps)] <- 0
  
  
  index <- attr(dataset, 'index')
  if (any(is.na(signal))) {
    na <- is.na(signal)
    icA <- ic
    icA[na] <- 0
    totalInfo <- sum(ic[index])
    kept <- sum(icA[index])
    discarded <- totalInfo - kept
    warning("Could not calculate signal for characters ",
            paste0(match(which(na), index), collapse = ', '),
            '; discarded ', signif(discarded), " bits from totals.")
    totalNoise <- sum(noise[index], na.rm = TRUE)
    totalSignal <- sum(signal[index], na.rm = TRUE)
    signalNoise <- totalSignal / totalNoise
    
    infoNeeded <- Log2Unrooted(length(dataset))
    infoOverkill <- totalInfo / infoNeeded
    
    message('`dataset` contains ',
            signif(totalInfo), ' bits (after discarding ',
            signif(discarded), '), of which ',
            signif(totalSignal), ' signal, ',
            signif(totalNoise), ' noise, ',
            signif(infoNeeded), ' needed.  ',
            'S:N = ', signif(signalNoise), "\n")
    
  } else {
    totalInfo <- sum(ic[index])
    totalNoise <- sum(noise[index])
    totalSignal <- sum(signal[index])
    signalNoise <- totalSignal / totalNoise
    discarded = 0
    
    infoNeeded <- Log2Unrooted(length(dataset))
    infoOverkill <- totalInfo / infoNeeded
    discarded <- originalInfo - totalInfo
    if (discarded < sqrt(.Machine$double.eps)) discarded <- 0
    
    message('dataset contains ',
            signif(totalInfo), ' bits',
            if (totalInfo != originalInfo) {
              paste0(' (after discarding ', signif(originalInfo - totalInfo),
                     ' bits)')
            }, ', of which ', 
            signif(totalSignal), ' signal, ',
            signif(totalNoise), ' noise, ',
            signif(infoNeeded), ' needed.  ',
            'S:N = ', signif(signalNoise), "\n")
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
  .Deprecated('ConcordantInformation()')
  ConcordantInformation(tree, dataset)
}

#' @rdname ConcordantInformation
#' @export
ConcordantInfo <- ConcordantInformation