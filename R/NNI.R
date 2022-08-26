#' Nearest neighbour interchange (NNI)
#' 
#' `NNI()`performs a single iteration of the nearest-neighbour interchange
#' algorithm; `RootedNNI()` retains the position of the root.
#' These functions are based on equivalents in the \pkg{phangorn} package.
#' `cNNI()` is an equivalent function coded in C, that runs much faster.
#' 
#' Branch lengths are not supported.
#' 
#' 
#' 
#' All nodes in a tree must be bifurcating; [ape::collapse.singles()] and
#' [ape::multi2di()] may help.
#' 
#' @param tree A tree of class `phylo`. 
#' For `cNNI()`, this must be a binary tree rooted on a single leaf, whose root
#' node is the lowest numbered internal node.
#' @template treeParam
#' @param edgeToBreak In (`Rooted`)`NNI()`, an optional integer specifying the
#' index of an edge to bisect/prune, generated randomly if not specified.
#' If \code{-1}, a complete list of all trees one step from the input tree
#' will be returned.
#' In `cNNI()`, an integer from zero to `nEdge(tree) - nTip(tree) - 2`,
#' specifying which internal edge to break.
#' 
#' @return Returns a tree with class \code{phylo} (if \code{returnAll = FALSE}) or 
#'         a set of trees, with class \code{multiPhylo} (if \code{returnAll = TRUE}).
#'
#' @references
#' The algorithm is summarized in
#'  \insertRef{Felsenstein2004}{TreeSearch}
#' 
#' 
#' @examples
#' tree <- TreeTools::BalancedTree(8)
#' # A random rearrangement
#' NNI(tree)
#' cNNI(tree)
#' 
#' # All trees one NNI rearrangement away
#' NNI(tree, edgeToBreak = -1)
#' 
#' # Manual random sampling
#' cNNI(tree, sample.int(14 - 8 - 1, 1), sample.int(2, 1))
#' 
#' # A specified rearrangement
#' cNNI(tree, 0, 0)
#' 
#' # If a tree may not be binary, collapse nodes with
#' tree <- TreeTools::MakeTreeBinary(tree)
#' 
#' # If a tree may be improperly rooted, use
#' tree <- TreeTools::RootTree(tree, 1)
#' 
#' # If a tree may exhibit unusual node ordering, this can be addressed with
#' tree <- TreeTools::Preorder(tree)
#' @template MRS
#'
#' @family tree rearrangement functions
#' @export
NNI <- function (tree, edgeToBreak = NULL) {
  edge <- tree$edge
  parent <- edge[, 1]
  StopUnlessBifurcating(parent)
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    child <- edge[, 2]
    nTips <- (length(parent) / 2L) + 1L
    samplable <- child > nTips
    # newEdges <- vapply(which(samplable), DoubleNNI, parent=parent, child=child, list(matrix(0L, nEdge, 2), matrix(0L, nEdge, 2)))
    newEdges <- unlist(lapply(which(samplable), DoubleNNI,
                              parent = parent, child = child), recursive = FALSE) # Quicker than vapply, surprisingly
    newTrees <- structure(lapply(newEdges, function (edges) {tree$edge <- edges; tree}), # Quicker than vapply, surprisingly
                          class = 'multiPhylo')
    # Return:
    newTrees
  } else {
    newEdge <- NNISwap(parent, edge[, 2], edgeToBreak = edgeToBreak)
    tree$edge <- cbind(newEdge[[1]], newEdge[[2]])
    
    # Return:
    tree
  }
}


#' @param whichSwitch Integer from zero to one, specifying which way to re-build
#' the broken internal edge.
#' 
#' @return `cNNI()` returns a tree of class `phylo`, rooted on the same leaf,
#' on which the specified rearrangement has been conducted.
#' @rdname NNI
#' @importFrom TreeTools NTip
#' @export
cNNI <- function (tree, edgeToBreak = NULL, whichSwitch = NULL) {
  edge <- tree$edge
  if (is.null(edgeToBreak)) edgeToBreak <- sample.int(dim(edge)[1] - NTip(tree) - 1L, 1L)
  if (is.null(whichSwitch)) whichSwitch <- sample.int(2L, 1L)
  tree$edge <- nni(edge, edgeToBreak, whichSwitch)
  
  # Return:
  tree
}
  
#' @describeIn NNI faster version that takes and returns parent and child parameters
#' @template treeParent
#' @template treeChild
#' @param nTips (optional) Number of tips.
#' @return `NNISwap()` returns a list containing two elements, corresponding in
#' turn to the  rearranged parent and child parameters.
#' @importFrom TreeTools SampleOne
#' @export
NNISwap <- function (parent, child, nTips = (length(parent) / 2L) + 1L,
                     edgeToBreak = NULL) {
  rootNode  <- nTips + 1L
  samplable <- child > nTips
  if (!any(samplable)) stop("Not enough edges to allow NNI rearrangement")
  
  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak must be an internal edge")
  }

  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash = FALSE)]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  RenumberEdges(parent, child)
}

## TODO use RenumberList
#' Double NNI
#' 
#' Returns the edge parameter of the two trees consistent with the speficied \acronym{NNI} rearrangement
#'
#' @template treeParent
#' @template treeChild
#' @template edgeToBreakParam
#'
#' @return the \code{tree$edge} parameter of the two trees consistent with the specified rearrangement
#'
#' @keywords internal
#' @importFrom TreeTools RenumberTree
#' @author Martin R. Smith
#' 
DoubleNNI <- function (parent, child, edgeToBreak) {
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2.3 <- which(parent == end2)
  ind2   <- ind2.3[1]
  ind3   <- ind2.3[2]

  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  child2 <- child
  childSwap <- child[newInd]
  child2[oldInd] <- childSwap
  
  newInd <- c(ind3, ind1)
  oldInd <- c(ind1, ind3)
  childSwap <- child[newInd]
  child[oldInd] <- childSwap
  
  nEdge <- length(parent)
  
  # Return:
  list(RenumberTree(parent, child), RenumberTree(parent, child2))
}

#' Rooted NNI 
#' @describeIn NNI Perform \acronym{NNI} rearrangement, retaining position of root
#' @export
RootedNNI <- function (tree, edgeToBreak=NULL) {
  edge <- tree$edge
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    parent <- edge[, 1]
    child  <- edge[, 2]
    nTips <- (length(parent) / 2L) + 1L
    rootNode <- nTips + 1L
    samplable <- parent != rootNode & child > nTips
    newEdges <- unlist(lapply(which(samplable), DoubleNNI, 
                              parent = parent, child = child), 
                       recursive = FALSE) # Quicker than vapply, surprisingly
    newTrees <- lapply(newEdges, function (edges) {tree$edge <- edges; tree}) # Quicker than vapply, surprisingly
    
    # Return:
    newTrees
  } else {
    newEdge <- RootedNNISwap(edge[, 1], edge[, 2], edgeToBreak=edgeToBreak)
    tree$edge <- cbind(newEdge[[1]], newEdge[[2]])
    
    # Return:
    tree
  }
}

#' @describeIn NNI faster version that takes and returns parent and child parameters
#' @return a list containing two elements, corresponding in turn to the rearranged parent and child parameters
#' @export
RootedNNISwap <- function (parent, child, nTips = (length(parent) / 2L) + 1L,
                           edgeToBreak = NULL) {
  rootNode <- nTips + 1L
  
  samplable <- parent != rootNode & child > nTips

  if (is.null(edgeToBreak)) { 
    edgeToBreak <- SampleOne(which(samplable))
  } else if (!samplable[edgeToBreak]) {
    stop("edgeToBreak cannot include a tip or the root node")
  }
  
  if (is.na(edgeToBreak)) stop("Cannot find a valid rearrangement")
  
  end1   <- parent[edgeToBreak]
  end2   <- child[edgeToBreak]
  ind1   <- which(parent == end1)
  ind1   <- ind1[ind1 != edgeToBreak][1]
  ind2   <- which(parent == end2)[sample.int(2L, 1L, useHash=FALSE)]
  
  newInd <- c(ind2, ind1)
  oldInd <- c(ind1, ind2)
  
  child_swap <- child[newInd]
  child[oldInd] <- child_swap
  RenumberEdges(parent, child)
}
