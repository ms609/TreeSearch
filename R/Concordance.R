#' Calculate site concordance factor
#' 
#' The site concordance factor (Minh, Hahn & Lanfear 2020) is a measure of the
#' strength of support that the dataset presents for a given split in a tree.
#' 
#' `QuartetConcordance()` is the proportion of quartets (sets of four leaves) that 
#' are decisive for a split which are also concordant with it.
#' For example, a quartet with the characters `0 0 0 1` is not decisive, as 
#' all relationships between those leaves are equally parsimonious.
#' But a quartet with characters `0 0 1 1` is decisive, and is concordant
#' with any tree that groups the first two leaves together to the exclusion
#' of the second.
#' 
# `ClusteringConcordance()` and `PhylogeneticConcordance()` respectively report
# the proportion of clustering information and phylogenetic information 
# (as defined in Vinh 2010, Smith 2020) within a dataset that is reflected in each split.
# These give smaller values because a split may be compatible with a character
# without being identical to it.
#TODO More thought / explanation needed.
#' 
#TODO Finally, `ProfileConcordance()` (to follow)
#' 
#' 
#' @template treeParam
#' @template datasetParam
#' 
#' 
#' 
#' @references 
#' 
#' \insertRef{Minh2020}{TreeSearch}
#' 
#' \insertRef{SmithDist}{TreeDist}
#' 
#' \insertRef{Vinh2010}{TreeDist}
#' @examples 
#' data('Lobo', package = 'TreeTools')
#' dataset <- Lobo.phy
#' tree <- TreeTools::NJTree(dataset)
#' qc <- QuartetConcordance(tree, dataset)
#' cc <- ClusteringConcordance(tree, dataset)
#' pc <- PhylogeneticConcordance(tree, dataset)
#' spc <- SharedPhylogeneticConcordance(tree, dataset)
#' mcc <- MutualClusteringConcordance(tree, dataset)
#' 
#' oPar <- par(mar = rep(0, 4), cex = 0.8)
#' plot(tree)
#' LabelSplits(tree, signif(qc, 3))
#' LabelSplits(tree, signif(cc, 3))
#' LabelSplits(tree, signif(pc, 3))
#' par(oPar)
#' 
#' pairs(cbind(qc, cc, pc, spc, mcc))
#' @template MRS
#' @importFrom ape keep.tip
#' @importFrom TreeTools as.Splits PhyDatToMatrix TipLabels
#' @importFrom Quartet SingleTreeQuartetAgreement
#' @name SiteConcordance
#' @family split support functions
#' @export
QuartetConcordance <- function (tree, dataset) {
  splits <- as.multiPhylo(as.Splits(tree))
  characters <- as.multiPhylo(dataset)
  
  status <- rowSums(vapply(characters, function (char) {
    trimmed <- lapply(splits, keep.tip, TipLabels(char))
    status <- SingleTreeQuartetAgreement(trimmed, char)
    s <- status[, 's']
    cbind(concordant = s, decisive = s + status[, 'd'])
  }, matrix(NA_real_, length(splits), 2)), dims = 2)
  
  # Return:
  status[, 1] / status[, 2]
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
  splits <- as.logical(as.Splits(tree))
  
  at <- attributes(dataset)
  cont <- at$contrast
  if ('-' %in% colnames(cont)) cont[cont[, '-'] > 0, ] <- 1
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
