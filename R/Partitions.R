#' Calculate bipartitions using phangorn's bipCPP function
#' 
#' @param orig Tree edge matrix (`someTree$edge`)
#' @param nTips number of tips
#' 
#' @return A list specifying, for each partition, the included taxa.
#' 
#' @export
#' @keywords internal
phangorn_bipCPP <- function(orig, nTips) {
  .Call(`C__TreeSearch_phangorn_bipCPP`, orig, nTips)
}

#' Tree2Splits
#' 
#' Converts a phylogenetic tree to an array of bipartition splits.
#' 
#' @param tr A tree of class \code{\link[ape:read.tree]{phylo}}, with tips 
#' bearing integer labels (i.e. `tr$tip.label == 1:N`).
#' @return Returns a two-dimensional array.  Columns correspond to unique
#'  bipartitions, named with the number of a node that denotes the partition.
#'  Rows correspond to tips `1:N`.
#'
#' @author Martin R. Smith
#' 
#' @examples Tree2Splits(ape::rtree(6, tip.label=1:6, br=NULL))
#'
#' @importFrom ape reorder.phylo
#' @export
Tree2Splits <- function (tr) {
  tr <- reorder.phylo(tr, 'postorder')
  tip_label <- tr$tip.label
  n_tip <- as.integer(length(tip_label))
  root <- length(tip_label) + 1
  bipartitions <- phangorn_bipCPP(tr$edge, n_tip)
  ret <- vapply(bipartitions[-seq_len(root)], 
                function (x) seq_len(n_tip) %in% x, 
                logical(n_tip))[seq_len(n_tip), , drop=FALSE]
  rownames(ret) <- tip_label
  colnames(ret) <- seq_len(ncol(ret)) + root
  
  ret <- UniqueSplits(ret)
  # Return:
  DropSingleSplits(ret)
}
#' @rdname Tree2Splits
#' @export
#' @keywords internal
Tree2Bipartitions <- Tree2Splits
