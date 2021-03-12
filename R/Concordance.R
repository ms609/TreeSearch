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

#' Convert object to multiPhylo class
#' @param x Object to be converted
#' @return `as.multiPhylo` returns an object of class `multiPhylo`
#' @examples 
#' library("TreeTools") # TODO move this function to TreeTools
#' as.multiPhylo(BalancedTree(8))
#' as.multiPhylo(list(BalancedTree(8), PectinateTree(8)))
#' data('Lobo')
#' as.multiPhylo(Lobo.phy)
#' @export
as.multiPhylo <- function (x) UseMethod('as.multiPhylo')

#' @rdname as.multiPhylo
#' @export
as.multiPhylo.phylo <- function (x) c(x)


#' @rdname as.multiPhylo
#' @export
as.multiPhylo.list <- function (x) structure(x, class = 'multiPhylo')

#' @rdname as.multiPhylo
#' @return `as.multiPhylo.phyDat()` returns a list of trees, each corresponding
#' to the partitions implied by each non-ambiguous character in `x`.
#' @importFrom ape read.tree
#' @export
as.multiPhylo.phyDat <- function (x) {
  at <- attributes(x)
  cont <- at$contrast
  if ('-' %in% colnames(cont)) cont[cont[, '-'] > 0, ] <- 1
  ambiguous <- rowSums(cont) != 1
  labels <- names(x)
  
  mat <- matrix(unlist(x), length(x), byrow = TRUE)
  mat[ambiguous] <- NA
  structure(apply(mat, 2, function (split) {
    a <- !is.na(split)
    aSplit <- split[a]
    tokens <- unique(aSplit)
    aLabels <- labels[a]
    if (length(tokens) == 1L) {
      read.tree(text = paste0('(', paste(aLabels, collapse = ', '), ');'))
    } else {
      read.tree(text = paste0('((', 
             paste(vapply(unique(aSplit), function (token) {
               paste(aLabels[aSplit == token], collapse = ', ')
              }, character(1)), collapse = '), ('),
            '));'))
    }
  }), class = 'multiPhylo')
}

#' @rdname as.multiPhylo
#' @importFrom ape read.tree
#' @importFrom TreeTools TipLabels
#' @export
as.multiPhylo.Splits <- function (x) {
  labels <- TipLabels(x)
  apply(as.logical(x), 1, function (a) {
      read.tree(text = paste0(
        '((', 
        paste0(labels[a], collapse = ','),
        '),(',
        paste0(labels[!a], collapse = ','),
        '));'))
    })
}