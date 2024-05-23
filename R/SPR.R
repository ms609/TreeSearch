#' @describeIn TBRWarning for SPR rearrangements
#' @keywords internal
#' @export
SPRWarning <- function (parent, child, error) {
  warning ("No SPR operation performed.\n  > ", error)
  
  # Return:
  list(parent, child)
}


#' Non-duplicate root
#'
#' Identify, for each edge, whether it denotes a different partition from
#' the root edge.
#' The first edge of the input tree must be a root edge; this can be
#' accomplished using `Preorder()`.
#'
#' This function is a copy of a deprecated ancestor in TreeTools; see 
#' [#32](https://github.com/ms609/TreeTools/issues/32).
#'
#' @template treeParent
#' @template treeChild
#' @template treeNEdgeOptional
#'
#' @return `.NonDuplicateRoot()` returns a logical vector of length `nEdge`,
#' specifying `TRUE` unless an edge identifies the same partition as
#' the root edge.
#'
#' @examples
#' tree <- TreeTools::Preorder(TreeTools::BalancedTree(8))
#' edge <- tree$edge
#' parent <- edge[, 1]
#' child <- edge[, 2]
#'
#' which(!.NonDuplicateRoot(parent, child))
#' @keywords internal
#' @importFrom TreeTools DescendantEdges
#' @template MRS
#' @export
.NonDuplicateRoot <- function(parent, child, nEdge = length(parent)) {
  notDuplicateRoot <- !logical(nEdge)
  rightSide <- DescendantEdges(edge = 1, parent, child, nEdge = nEdge)
  nEdgeRight <- sum(rightSide)
  if (nEdgeRight == 1) {
    notDuplicateRoot[2] <- FALSE
  } else if (nEdgeRight == 3) {
    notDuplicateRoot[4] <- FALSE
  } else {
    notDuplicateRoot[1] <- FALSE
  }
  notDuplicateRoot
}

#' Subtree pruning and rearrangement (SPR)
#'
#' Perform one \acronym{SPR} rearrangement on a tree
#' 
#' Equivalent to `kSPR()` in the \pkg{phangorn} package, but faster.
#' Note that rearrangements that only change the position of the root WILL be returned by 
#' \code{SPR}.  If the position of the root is irrelevant (as in Fitch parsimony, for example)
#' then this function will occasionally return a functionally equivalent topology.  
#' \code{RootIrrelevantSPR} will search tree space more efficiently in these cases.
#' Branch lengths are not (yet) supported.
#'
#' All nodes in a tree must be bifurcating; [ape::collapse.singles] and
#' [ape::multi2di] may help.
#'
#' @template treeParam
#' @param edgeToBreak the index of an edge to bisect, generated randomly if not specified.
#' @param mergeEdge the index of an edge on which to merge the broken edge.
#' @return This function returns a tree in \code{phyDat} format that has undergone one \acronym{SPR} iteration.
#' 
#' @references The \acronym{SPR} algorithm is summarized in
#'  \insertRef{Felsenstein2004}{TreeSearch}
#' 
#' @author Martin R. Smith
#' 
#' @seealso 
#' - [`RootedSPR()`]: useful when the position of the root node should be retained.
#' @family tree rearrangement functions
#' 
#' @examples{
#' tree <- ape::rtree(20, br=FALSE)
#' SPR(tree)
#' }
#' @importFrom ape root
#' @importFrom TreeTools Preorder
#' @export
SPR <- function(tree, edgeToBreak = NULL, mergeEdge = NULL) {
  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") {
    tree <- Preorder(tree)
  }
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  StopUnlessBifurcating(parent)
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    child <- edge[, 2]
    nEdge <- length(parent)
    stop("Negative edgeToBreak not yet supported; on TODO list for next release")
    notDuplicateRoot <- .NonDuplicateRoot(parent, child, nEdge)
    # Return:
    unique(unlist(lapply(which(notDuplicateRoot), AllSPR,
      parent=parent, child=child, nEdge=nEdge, notDuplicateRoot=notDuplicateRoot),
      recursive=FALSE)) # TODO the fact that we need to use `unique` indicates that 
                         #      we're being inefficient here.
  } else {
    newEdge <- SPRSwap(parent, edge[, 2], edgeToBreak = edgeToBreak,
                       mergeEdge = mergeEdge)
    tree[["edge"]] <- cbind(newEdge[[1]], newEdge[[2]])
    # Return:
    tree
  }
}

#' @rdname SPR
#' @return `TBRMoves()` returns a list of all trees one SPR move away from
#'  `tree`, with edges and nodes in preorder, rooted on the first-labelled tip.
#' @export
SPRMoves <- function (tree, edgeToBreak = integer(0)) UseMethod("SPRMoves")

#' @rdname SPR
#' @importFrom TreeTools Preorder RootTree
#' @export
SPRMoves.phylo <- function (tree, edgeToBreak = integer(0)) {
  tree <- Preorder(RootTree(tree, tree[["tip.label"]][1]))
  edges <- unique(.all_spr(tree[["edge"]], edgeToBreak))
  structure(lapply(edges, function (edg) {
    tree[["edge"]] <- edg
    tree
  }), class = "multiPhylo", tip.label = tree[["tip.label"]])
}

# error checking for all_spr
.all_spr <- function (edge, break_order) {
  nEdge <- dim(edge)[1]
  if (nEdge < 5) {
    stop("No SPR rearrangements possible on a tree with < 5 edges");
  }
  nInternal <- floor(nEdge / 2)
  nTip <- nInternal + 1L
  rootNode <- nTip + 1L
  
  if (edge[1] != rootNode) {
    stop("edge[1,] must connect root to leaf. Try Preorder(root(tree)).");
  }
  if (edge[2] != rootNode) {
    stop("edge[2,] must connect root to leaf. Try Preorder(root(tree)).");
  }
  # Return:
  all_spr(edge, break_order)
}

#' @rdname SPR
#' @export
SPRMoves.matrix <- function (tree, edgeToBreak = integer(0)) {
  tree <- Preorder(RootTree(tree, 1))
  unique(.all_spr(tree, edgeToBreak))
}

## TODO Do edges need to be pre-ordered before coming here?
#' @describeIn SPR faster version that takes and returns parent and child parameters
#' @template treeParent
#' @template treeChild
#' @template treeNEdgeOptional
#' @param nNode (optional) Number of nodes.
#' @return a list containing two elements, corresponding in turn to the
#'  rearranged parent and child parameters
#' @importFrom TreeTools DescendantEdges
#' @export
SPRSwap <- function (parent, child, nEdge = length(parent), nNode = nEdge / 2L,
                     edgeToBreak = NULL, mergeEdge = NULL) {
  
  if (nEdge < 5) {
    # TODO we need to re-root this tree...
    return(list(parent, child))
  }
  
  notDuplicateRoot <- .NonDuplicateRoot(parent, child, nEdge)
  
  if (is.null(edgeToBreak)) {
    # Pick an edge at random
    edgeToBreak <- SampleOne(which(notDuplicateRoot), len=nEdge - 1L)
  } else if (edgeToBreak > nEdge) {
    return(SPRWarning(parent, child, "edgeToBreak > nEdge"))
  } else if (edgeToBreak < 1) {
    return(SPRWarning(parent, child, "edgeToBreak < 1"))
  }
  # Breaking a single edge
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
    
  edgesCutAdrift <- DescendantEdges(edge = edgeToBreak, parent = parent,
                                    child = child, nEdge = nEdge)
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters
  if (breakingRootEdge <- !any(brokenEdgeParent)) { 
    # Edge to break is the Root Node.
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
  
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (mergeEdge > nEdge) {
      return(SPRWarning(parent, child, "mergeEdge value > number of edges"))
    } else if (length(mergeEdge) !=  1) {
        return(SPRWarning(
          parent, child, 
          paste0("mergeEdge value ", paste(mergeEdge, collapse="|"),
                 " invalid; must be NULL or a vector of length 1\n")
        ))
    } else if (nearBrokenEdge[mergeEdge]) {
      return(SPRWarning(parent, child,
                        "Selected mergeEdge will not change tree topology."))
    } else if (DescendantEdges(edge = edgeToBreak, parent = parent,
                               child = child, nEdge = nEdge)[[mergeEdge]]) {
      stop("mergeEdge is within pruned subtree")
    }
  } else {
    mergeEdge <- which(!nearBrokenEdge & !edgesOnAdriftSegment & notDuplicateRoot)
    nCandidates <- length(mergeEdge)
    #####Assert(nCandidates > 0)
    if (nCandidates > 1) {
      mergeEdge <- SampleOne(mergeEdge, len = nCandidates)
    }
  }
  
  if (breakingRootEdge) {
    parent[brokenRootDaughters] <- brokenEdge.parentNode
    spareNode <- child[brokenEdgeSister]
    child [brokenEdgeSister] <- child[mergeEdge]
    parent[brokenEdge | brokenEdgeSister] <- spareNode
    child[mergeEdge] <- spareNode
  } else {
    parent[brokenEdgeSister] <- parent[brokenEdgeParent]
    parent[brokenEdgeParent] <- parent[[mergeEdge]]
    parent[[mergeEdge]] <- brokenEdge.parentNode
  }
  
  #####Assert(identical(unique(table(parent)), 2L))
  #####Assert(identical(unique(table(child)),  1L))
  # Return:
  RenumberEdges(parent, child)
}


#' `cSPR()` expects a tree rooted on a single tip. 
#' @template treeParam
#' @param whichMove Integer specifying which SPR move index to perform.
#' @examples 
#' tree <- TreeTools::BalancedTree(8)
#' 
#' # Tree must be rooted on leaf
#' tree <- TreeTools::RootTree(tree, 1)
#'
#' # Random rearrangement
#' cSPR(tree)
#'
#' # Specific rearrangement
#' cSPR(tree, 9)
#' @template MRS
#' @importFrom TreeTools NTip
#' @export
cSPR <- function (tree, whichMove = NULL) {
  edge <- tree[["edge"]]
  if (is.null(whichMove)) whichMove <- sample.int(2147483647L, 1L)
  tree[["edge"]] <- spr(edge, whichMove)
  
  # Return:
  tree
}


#' All SPR trees
#'
#' @template treeParent
#' @template treeChild
#' @template treeNEdge
#' @param notDuplicateRoot logical vector of length `nEdge`, specifying for each
#' edge whether it is the second edge leading to the root (in which case
#' its breaking will be equivalent to breaking the other root edge... 
#' except insofar as it moves the position of the root.)
#' @template edgeToBreakParam
#' 
#' @return `AllSPR()` returns a list of edge matrices for all trees one SPR 
#' rearrangement from the starting tree
#'
#' @author Martin R. Smith
#' 
AllSPR <- function (parent, child, nEdge, notDuplicateRoot, edgeToBreak) {

  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
    
  edgesCutAdrift <- DescendantEdges(edge = edgeToBreak, parent = parent,
                                    child = child, nEdge = nEdge)
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | brokenEdgeSister | brokenEdgeParent |
    brokenEdgeDaughters
  if (breakingRootEdge <- !any(brokenEdgeParent)) { 
    # Edge to break is the Root Node.
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
  
  mergeEdges <- which(!nearBrokenEdge & !edgesOnAdriftSegment & 
                        notDuplicateRoot)
  nCandidates <- length(mergeEdges)
  
  if (breakingRootEdge) {
    newEdges <- lapply(mergeEdges, function (mergeEdge) {
      newParent <- parent
      newChild <- child
      newParent[brokenRootDaughters] <- brokenEdge.parentNode
      newChild [brokenEdgeSister] <- child[mergeEdge]
      newParent[brokenEdge | brokenEdgeSister] <- child[brokenEdgeSister]
      newChild[mergeEdge] <- child[brokenEdgeSister]
      # Return:
      RenumberTree(newParent, newChild)
    }) # lapply faster than vapply
  } else {
    newEdges <- lapply(mergeEdges, function (mergeEdge) {
      newParent <- parent
      newParent[brokenEdgeSister] <- parent[brokenEdgeParent]
      newParent[brokenEdgeParent] <- newParent[mergeEdge]
      newParent[mergeEdge] <- brokenEdge.parentNode
      # Return:
      RenumberTree(newParent, child)
    }) # lapply faster than vapply
  }
  # Return:
  lapply(newEdges, function (newEdge) {tree[["edge"]] <- newEdge; tree})
}

#' Rooted SPR 
#' @describeIn SPR Perform \acronym{SPR} rearrangement, retaining position of root
#' @importFrom TreeTools Preorder
#' @export
RootedSPR <- function(tree, edgeToBreak = NULL, mergeEdge = NULL) {
  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") {
    tree <- Preorder(tree)
  }
  edge <- tree[["edge"]]
  newEdge <- RootedSPRSwap(edge[, 1], edge[, 2], edgeToBreak = edgeToBreak,
                           mergeEdge = mergeEdge)
  tree[["edge"]] <- cbind(newEdge[[1]], newEdge[[2]])
  return (tree)
}

## TODO Do edges need to be pre-ordered before coming here?
#' @describeIn SPR faster version that takes and returns parent and child parameters
#' @return a list containing two elements, corresponding in turn to the rearranged parent and child parameters
#' @export
RootedSPRSwap <- function (parent, child, nEdge = length(parent), nNode = nEdge / 2L,
                     edgeToBreak=NULL, mergeEdge=NULL) {
  
  if (nEdge < 5) return (SPRWarning(parent, child, "Too few tips to rearrange."))
  
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  breakable <- !logical(nEdge) & !rootEdges
  
  
  if (!is.null(edgeToBreak) && edgeToBreak == -1) {
    notDuplicateRoot <- .NonDuplicateRoot(parent, child, nEdge)
    return(unique(unlist(lapply(which(breakable), AllSPR,
      parent=parent, child=child, nEdge=nEdge, notDuplicateRoot=notDuplicateRoot),
      recursive=FALSE))) # TODO the fact that we need to use `unique` indicates that 
                         #      we're being inefficient here.
  }
  
  rightSide <- DescendantEdges(edge = 1, parent, child, nEdge = nEdge)
  leftSide  <- !rightSide
  nEdgeRight <- which(rootEdges)[2] - 1
  nEdgeLeft <- nEdge - nEdgeRight
  if (nEdgeRight < 4) {
    if (nEdgeLeft < 4) {
      return(SPRWarning(parent, child, "No rearrangement possible with this root position."))
    }

    breakable <- breakable & !rightSide
    rightHalfOfLeftSide <- DescendantEdges(edge = nEdgeRight + 2L, parent,
                                           child, nEdge = nEdge)
     leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide & !rootEdges
      if (sum(rightHalfOfLeftSide) == 1) breakable[nEdgeRight + 3] <- FALSE
      if (sum( leftHalfOfLeftSide) == 1) breakable[nEdgeRight + 2] <- FALSE
  } else {
    if (nEdgeLeft < 4) {
      breakable <- breakable & rightSide
    } else {
      rightHalfOfLeftSide <- DescendantEdges(edge = nEdgeRight + 2L, parent,
                                             child, nEdge = nEdge)
       leftHalfOfLeftSide <- leftSide & !rightHalfOfLeftSide & !rootEdges
      if (sum(rightHalfOfLeftSide) == 1) breakable[nEdgeRight + 3] <- FALSE
      if (sum( leftHalfOfLeftSide) == 1) breakable[nEdgeRight + 2] <- FALSE
    }
    rightHalfOfRightSide <- DescendantEdges(edge = 2L, parent = parent, 
                                            child = child, nEdge = nEdge)
     leftHalfOfRightSide <- rightSide & !rightHalfOfRightSide & !rootEdges
    if (sum(rightHalfOfRightSide) == 1) {
      breakable[3] <- FALSE
    }
    if (sum(leftHalfOfRightSide) == 1) {
      breakable[2] <- FALSE
    }
  }  
  
  if (is.null(edgeToBreak)) {
    # Pick an edge at random
    edgeToBreak <- SampleOne(which(breakable))
  } else {
    if (!breakable[edgeToBreak]) {
      return(SPRWarning(
        parent, child,
        paste("Nowhere to regraft if pruning on edge", edgeToBreak)))
    } else if (edgeToBreak > nEdge) {
      return(SPRWarning(parent, child, "edgeToBreak > nEdge"))
    } else if (edgeToBreak < 1) {
      return(SPRWarning(parent, child, "edgeToBreak < 1"))
    }
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesCutAdrift <- DescendantEdges(edge = edgeToBreak, parent= parent,
                                    child = child, nEdge = nEdge)
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | 
    brokenEdgeSister | 
    brokenEdgeParent | 
    brokenEdgeDaughters
    
  if (!is.null(mergeEdge)) { # Quick sanity checks
    if (mergeEdge > nEdge) {
      return(SPRWarning(parent, child, "mergeEdge value > number of edges"))
    } else if (length(mergeEdge) !=  1) {
        return(SPRWarning(
          parent, child, 
          paste0("mergeEdge value ", paste(mergeEdge, collapse="|"),
                 " invalid; must be NULL or a vector of length 1\n")
        ))
    } else if (nearBrokenEdge[[mergeEdge]]) {
      return(SPRWarning(parent, child,
                        "Selected mergeEdge will not change tree topology."))
    } else if (DescendantEdges(edge = edgeToBreak, parent, child,
                               nEdge = nEdge)[[mergeEdge]]) {
      stop("mergeEdge is within pruned subtree")
    }
  } else {
    edgesOnThisSide <- if (rightSide[edgeToBreak]) rightSide else leftSide
    mergeEdge <- which(edgesOnThisSide & !nearBrokenEdge & !edgesOnAdriftSegment)
    nCandidates <- length(mergeEdge)
    if (nCandidates > 1) mergeEdge <- SampleOne(mergeEdge, len=nCandidates)
  }
  
  parent[brokenEdgeSister] <- parent[brokenEdgeParent]
  parent[brokenEdgeParent] <- parent[mergeEdge]
  parent[mergeEdge] <- brokenEdge.parentNode
  
  #####Assert(identical(unique(table(parent)), 2L))
  #####Assert(identical(unique(table(child)),  1L))
  
  # Return:
  RenumberEdges(parent, child)
}
