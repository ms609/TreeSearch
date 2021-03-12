#' Calculate site concordance factor
#' 
#' The site concordance factor (Binh, Hahn & Lanfear 2020) is a measure of the
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
#' `ClusteringConcordance()` and `PhylogeneticConcordance()` respectively report
#' the proportion of clustering information and phylogenetic information 
#' (as defined in Smith 2020) within a dataset that is reflected in each split.
#' These give smaller values because a split may be compatible with a character
#' without being identical to it.
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
#' \insertRef{Binh2020}{TreeSearch}
#' 
#' \insertRef{Smith2020}{TreeSearch}
#' @examples 
#' data('Lobo')
#' dataset <- Lobo.phy
#' tree <- NJTree(dataset)
#' qc <- QuartetConcordance(tree, dataset)
#' cc <- ClusteringConcordance(tree, dataset)
#' pc <- PhylogeneticConcordance(tree, dataset)
#' 
#' oPar <- par(mar = rep(0, 4), cex = 0.8)
#' plot(tree)
#' LabelSplits(tree, signif(qc, 3))
#' LabelSplits(tree, signif(cc, 3))
#' LabelSplits(tree, signif(pc, 3))
#' par(oPar)
#' 
#' pairs(cbind(qc, cc, pc))
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

#' @rdname SiteConcordance
#' @importFrom TreeDist ClusteringEntropy MutualClusteringInfo
#' @export
ClusteringConcordance <- function (tree, dataset) {
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
#' @importFrom TreeDist ClusteringInfo SharedPhylogeneticInfo
#' @export
PhylogeneticConcordance <- function (tree, dataset) {
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
