#' TBR Warning
#' Print a warning and return given tree
#'
#' @template treeParent
#' @template treeChild
#' @param error error message to report
#'
#' @return A list with the entries `parent`, `child`.
#' @examples
#' suppressWarnings(TBRWarning(0, 0, "Message text")) # will trigger warning
#' 
#' @author Martin R. Smith
#' @keywords internal
#' @export
TBRWarning <- function (parent, child, error) {
  warning ("No TBR operation performed.\n  > ", error)
  # Return:
  list(parent, child)
}

#' Tree bisection and reconnection (TBR)
#'
#' \code{TBR} performs a single random \acronym{TBR} iteration.
#'
#' Branch lengths are not (yet) supported.
#' 
#' All nodes in a tree must be bifurcating; [ape::collapse.singles] and
#' [ape::multi2di] may help.
#' 
#' @param tree A bifurcating tree of class \code{\link[ape]{phylo}}, with all nodes resolved;
#' @param edgeToBreak (optional) integer specifying the index of an edge to bisect/prune,
#' generated randomly if not specified.  
#' Alternatively, set to \code{-1} to return a complete list
#' of all trees one step from the input tree.
#' @param mergeEdges (optional) vector of length 1 or 2, listing edge(s) to be joined:
#'                   In SPR, this is where the pruned subtree will be reconnected.
#'                   In TBR, these edges will be reconnected (so must be on opposite
#'                   sides of \code{edgeToBreak}); if only a single edge is specified,
#'                   the second will be chosen at random
#' 
#' @return `TBR()` returns a tree in \code{phyDat} format that has undergone one
#' \acronym{TBR} iteration.
#' @references The \acronym{TBR} algorithm is summarized in
#' \insertRef{Felsenstein2004}{TreeSearch}
#' 
#' @examples
#' library("ape")
#' tree <- rtree(20, br=NULL)
#' TBR(tree)
#' @template MRS
#' 
#' @family tree rearrangement functions
#' @seealso [`RootedTBR()`]: useful when the position of the root node should be retained.
#' @importFrom ape root
#' @importFrom TreeTools DescendantEdges Preorder
#' @export
TBR <- function(tree, edgeToBreak = NULL, mergeEdges = NULL) {
  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") {
    tree <- Preorder(tree)
    if (!is.null(edgeToBreak)) {
      warning("Edge numbering modified as tree not in preorder;
               edgeToBreak and mergeEdges ignored.")
      edgeToBreak <- mergeEdges <- NULL
    }
  }
  
  edge <- tree[["edge"]]
  StopUnlessBifurcating(edge[, 1])
  newEdge <- TBRSwap(
    parent = edge[, 1],
    child = edge[, 2],
    edgeToBreak = edgeToBreak,
    mergeEdges = mergeEdges
  )
  tree[["edge"]] <- cbind(newEdge[[1]], newEdge[[2]])
  tree
}

#' @rdname TBR 
#' @return `TBRMoves()` returns a `multiPhylo` object listing all trees one
#'  \acronym{TBR} move away from `tree`, with edges and nodes in preorder,
#'  rooted on the first-labelled tip.
#' @export
TBRMoves <- function (tree, edgeToBreak = integer(0)) UseMethod("TBRMoves")

#' @rdname TBR 
#' @importFrom TreeTools Preorder RootTree
#' @export
TBRMoves.phylo <- function (tree, edgeToBreak = integer(0)) {
  tree <- Preorder(RootTree(tree, 1))
  edges <- unique(all_tbr(tree[["edge"]], edgeToBreak))
  structure(lapply(edges, function (edg) {
    tree[["edge"]] <- edg
    tree
  }), class = "multiPhylo", tip.label = tree[["tip.label"]])
}

#' @rdname TBR
#' @export
TBRMoves.matrix <- function (tree, edgeToBreak = integer(0)) {
  tree <- Preorder(RootTree(tree, 1))
  allMoves <- all_tbr(tree, edgeToBreak)
  unique(allMoves)
}

## TODO Do edges need to be pre-ordered before coming here?
#' @describeIn TBR faster version that takes and returns parent and child
#'  parameters
#' @template treeParent
#' @template treeChild
#' @param nEdge (optional) Number of edges.
#' @return `TBRSwap()` returns a list containing two elements corresponding
#' to the rearranged `parent` and `child` parameters.
#'  
#' @importFrom TreeTools EdgeAncestry
#' @export
TBRSwap <- function(parent, child, nEdge = length(parent),
                    edgeToBreak = NULL,
                    mergeEdges = NULL) {
  if (nEdge < 5) {
    return (list(parent, child)) #TODO do we need to re-root this tree?
  }
  
  # Pick an edge at random
  allEdges <- seq_len(nEdge - 1L) + 1L # Only include one root edge
  not1 <- !logical(nEdge)
  not1[[1]] <- FALSE
  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(allEdges, len = nEdge - 1L)
  } else {
    if (edgeToBreak > nEdge) {
      return(TBRWarning(parent, child, "edgeToBreak > nEdge"))
    }
    if (edgeToBreak < 1) {
      return(TBRWarning(parent, child, "edgeToBreak < 1"))
    }
    if (edgeToBreak == 1) {
      edgeToBreak <- which(parent == parent[[1]])[-1] # Use other side of root
    }
  }
  # More efficient than tabulate or c(logical(), TRUE, logical())
  brokenEdge <- seq_along(parent) == edgeToBreak
  
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) {
      return(TBRWarning(parent, child, "mergeEdges value > number of edges"))
    } else if (length(mergeEdges) > 2 || length(mergeEdges) == 0) {
      return(TBRWarning(
        parent, child,
        paste0("mergeEdges value ", paste(mergeEdges, collapse = "|"),
               " invalid; must be NULL or a vector of length 1 or 2\n  ")
      ))
    } else if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) {
      return(TBRWarning(parent, child, "mergeEdges values must differ"))
    }
  }
  
  edgesCutAdrift <- DescendantEdges(edge = edgeToBreak, parent = parent,
                                    child = child, nEdge = nEdge)
  edgesRemaining <- !edgesCutAdrift & !brokenEdge
  
  brokenEdgeParent <- child == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdge | 
    brokenEdgeSister | 
    brokenEdgeParent | 
    brokenEdgeDaughters
  
  if (breakingRootEdge <- !any(brokenEdgeParent)) { 
    # Edge to break is the Root Node.
    brokenRootDaughters <- parent == child[brokenEdgeSister]
    nearBrokenEdge <- nearBrokenEdge | brokenRootDaughters
  }
  
  if (is.null(mergeEdges)) {
    candidateEdges <- which(!nearBrokenEdge & not1)
    nCandidates <- length(candidateEdges)
    mergeEdges <- if (nCandidates > 1) {
      SampleOne(candidateEdges, len = nCandidates) 
    } else {
      candidateEdges
    }
  }
  if (length(mergeEdges) == 1) {
    if (edgesCutAdrift[[mergeEdges]]) {
      adriftReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[[mergeEdges]]) {
        samplable <- which(!edgesCutAdrift & !nearBrokenEdge & not1)
      } else {
        samplable <- which(!edgesCutAdrift & not1)
        if (all(edgesCutAdrift == not1) && breakingRootEdge) {
          samplable <- 1
        }
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) {
        return(TBRWarning(
          parent, child,
          "No reconnection site would modify the tree; check mergeEdge"
        ))
      }
      rootedReconnectionEdge <- if (nSamplable == 1) {
        samplable
      } else {
        SampleOne(samplable, len = nSamplable)
      }
      #### message(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(edgesCutAdrift & !nearBrokenEdge & not1)
      } else {
        samplable <- which(edgesCutAdrift & not1)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) {
        return(TBRWarning(
          parent, child,
          "No reconnection site would modify the tree; check mergeEdge"
        ))
      }
      adriftReconnectionEdge <- if (nSamplable == 1) {
        samplable
      } else {
        SampleOne(samplable)
      }
      #### message(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesCutAdrift[mergeEdges]
    if (sum(whichAdrift) != 1) {
      return(TBRWarning(
        parent, child,
        paste("Invalid edges selected to merge:",
              mergeEdges[[1]], mergeEdges[[2]])
      ))
    }
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] &&
     nearBrokenEdge[adriftReconnectionEdge]) {
    return(TBRWarning(
      parent, child, "Selected mergeEdges will not change tree topology."
    ))
  }
  #### edgelabels(edge = edgeToBreak, bg="orange", cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg="cyan")    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg="magenta") #### DEBUGGING AID
  
  if (!nearBrokenEdge[[adriftReconnectionEdge]]) {
    edgesToInvert <- EdgeAncestry(adriftReconnectionEdge, parent, child,
                                  stopAt = edgeToBreak) & !brokenEdge
    #### which(edgesToInvert)
    if (any(edgesToInvert)) {
      tmp <- parent[edgesToInvert]
      parent[edgesToInvert] <- child[edgesToInvert]
      child[edgesToInvert] <- tmp
    }
    reconnectionSideEdges <- edgesToInvert
    reconnectionSideEdges[adriftReconnectionEdge] <- TRUE
    
    repurposedDaughterEdge <- brokenEdgeDaughters & reconnectionSideEdges
    spareDaughterEdge      <- brokenEdgeDaughters & !reconnectionSideEdges
    #########Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    #########Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[[adriftReconnectionEdge]] <- brokenEdge.childNode
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    if (breakingRootEdge) {
      parent[brokenRootDaughters] <- brokenEdge.parentNode
      spareNode <- child[brokenEdgeSister]
      child [brokenEdgeSister] <- child[rootedReconnectionEdge]
      parent[brokenEdge | brokenEdgeSister] <- spareNode
      child[rootedReconnectionEdge] <- spareNode
    } else {
      parent[brokenEdgeSister] <- parent[brokenEdgeParent]
      parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
      parent[rootedReconnectionEdge] <- brokenEdge.parentNode
    }
  }
  
  #########Assert(identical(unique(table(parent)), 2L))
  #########Assert(identical(unique(table(child)),  1L))
  # Return:
  RenumberEdges(parent, child)
}

#' Rooted TBR 
#' @describeIn TBR Perform \acronym{TBR} rearrangement, retaining position of root
#' @importFrom TreeTools Preorder
#' @export
RootedTBR <- function(tree, edgeToBreak = NULL, mergeEdges = NULL) {
  if (is.null(treeOrder <- attr(tree, "order")) || treeOrder != "preorder") {
    tree <- Preorder(tree)
  }
  edge   <- tree[["edge"]]
  edgeList <- RootedTBRSwap(edge[, 1], edge[, 2], 
                            edgeToBreak=edgeToBreak, mergeEdges=mergeEdges)
  tree[["edge"]] <- cbind(edgeList[[1]], edgeList[[2]])
  tree
}

#' @describeIn TBR faster version that takes and returns parent and child parameters
#' @importFrom TreeTools EdgeAncestry
#' @export
RootedTBRSwap <- function (parent, child, nEdge=length(parent), 
                           edgeToBreak = NULL, mergeEdges = NULL) {
  if (nEdge < 5) return (TBRWarning(parent, child, "Fewer than 4 tips"))
  nTips <- (nEdge / 2L) + 1L
  rootNode <- parent[1]
  rootEdges <- parent == rootNode
  rightTree <- DescendantEdges(parent, child, edge = 1, nEdge = nEdge)
  selectableEdges <- !rootEdges
  if (sum( rightTree) < 4) {
    selectableEdges[ rightTree] <- FALSE
  } else if (sum( rightTree) < 6) {
    rightChild <- child[1]
    rightGrandchildEdges   <- parent==rightChild
    rightGrandchildren     <- child[rightGrandchildEdges]
    rightGrandchildrenTips <- rightGrandchildren <= nTips
    selectableEdges[which(rightGrandchildEdges)[!rightGrandchildrenTips]] <- FALSE  
  }
  if (sum(!rightTree) < 4) {
    selectableEdges[!rightTree] <- FALSE
  } else if (sum(!rightTree) < 6) {
     leftChild <- child[rootEdges][2]
     leftGrandchildEdges   <- parent==leftChild
     leftGrandchildren     <- child[ leftGrandchildEdges]
     leftGrandchildrenTips <-  leftGrandchildren <= nTips
     selectableEdges[which( leftGrandchildEdges)[! leftGrandchildrenTips]] <- FALSE  
  }
  
  if (!any(selectableEdges)) return(TBRWarning(parent, child, "No opportunity to rearrange tree due to root position"))

  if (is.null(edgeToBreak)) {
    edgeToBreak <- SampleOne(which(selectableEdges)) # Pick an edge at random
  } else {
    if (edgeToBreak > nEdge) return(TBRWarning(parent, child, "edgeToBreak > nEdge"))
    if (edgeToBreak < 1) return(TBRWarning(parent, child, "edgeToBreak < 1"))
    if (rootEdges[edgeToBreak]) return(TBRWarning(parent, child, "RootedTBR cannot break root edge; try TBR"))
    if (!selectableEdges[edgeToBreak]) return(TBRWarning(parent, child, paste("Breaking edge", edgeToBreak,
                                              "does not allow a changing reconnection")))
  }
  repeat {
    edgeInRight <- rightTree[edgeToBreak]
    subtreeWithRoot <- if (edgeInRight) rightTree else !rightTree
    subtreeEdges <- !rootEdges & subtreeWithRoot
    edgesCutAdrift <- DescendantEdges(parent, child, edge = edgeToBreak,
                                      nEdge = nEdge)
    if (sum(edgesCutAdrift) > 2) {
      break;
    }
    if (sum(subtreeEdges, -edgesCutAdrift) > 2) {
      break; # the edge itself, and somewheres else
    }
    # TODO check that all expected selections are valid
    selectableEdges[edgeToBreak] <- FALSE
    ###Assert(any(selectableEdges))
    edgeToBreak <- SampleOne(which(selectableEdges))
  }
  brokenEdge <- seq_along(parent) == edgeToBreak
  brokenEdge.parentNode <- parent[edgeToBreak]
  brokenEdge.childNode  <-  child[edgeToBreak]
  
  edgesRemaining <- !edgesCutAdrift & subtreeEdges
  edgesOnAdriftSegment <- edgesCutAdrift | brokenEdge
  
  if (!is.null(mergeEdges)) { # Quick sanity checks
    if (any(mergeEdges > nEdge)) return(TBRWarning(parent, child, "mergeEdges value > number of edges"))
    if (length(mergeEdges) > 2 || length(mergeEdges) == 0) 
        return(TBRWarning(parent, child, paste0("mergeEdges value ", paste(mergeEdges, collapse="|"),  
               " invalid; must be NULL or a vector of length 1 or 2\n  ")))
    if (length(mergeEdges) == 2 && mergeEdges[1] == mergeEdges[2]) 
      return(TBRWarning(parent, child, "mergeEdges values must differ"))
    if (!all(subtreeWithRoot[mergeEdges])) return(TBRWarning(parent, child, paste("mergeEdges", 
          mergeEdges[1], mergeEdges[2], "not on same side of root as edgeToBreak", edgeToBreak)))
  }  
  
  brokenEdgeParent <- child  == brokenEdge.parentNode
  brokenEdgeSister <- parent == brokenEdge.parentNode & !brokenEdge
  
  brokenEdgeDaughters <- parent == brokenEdge.childNode
  nearBrokenEdge <- brokenEdgeSister | brokenEdgeParent | brokenEdgeDaughters | brokenEdge
  ###Assert(any(brokenEdgeParent))
  
  if (is.null(mergeEdges)) {
    mergeEdges  <- which(subtreeEdges & !nearBrokenEdge)
    nCandidates <- length(mergeEdges)
    if (nCandidates > 1) mergeEdges <- SampleOne(mergeEdges, len=nCandidates)
  }
  if (length(mergeEdges) == 0) {
    return(TBRWarning(parent, child, paste("Breaking edge", edgeToBreak, "does not allow any new reconnections whilst preserving root position.")))
  }
  if (length(mergeEdges) == 1) {
    if (edgesOnAdriftSegment[mergeEdges]) {
      adriftReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & !edgesOnAdriftSegment)
        ###Assert(length(samplable) > 0)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
      rootedReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### message(" - Selected rooted Reconnection Edge: ", rootedReconnectionEdge, "\n")  #### DEBUGGING AID
    } else {
      rootedReconnectionEdge <- mergeEdges
      if (nearBrokenEdge[mergeEdges]) {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment & !nearBrokenEdge)
      } else {
        samplable <- which(subtreeEdges & edgesOnAdriftSegment)
      }
      nSamplable <- length(samplable)
      if (nSamplable == 0) return(TBRWarning(parent, child, "No reconnection site would modify the tree; check mergeEdge"))
      adriftReconnectionEdge <- if (nSamplable == 1) samplable else SampleOne(samplable, len=nSamplable)
      #### message(" - Selected adrift Reconnection Edge: ", adriftReconnectionEdge, "\n") #### DEBUGGING AID
    }
  } else {
    whichAdrift <- edgesOnAdriftSegment[mergeEdges]
    if (sum(whichAdrift) != 1) return(TBRWarning(parent, child, paste("Invalid edges selected to merge:",
            mergeEdges[1], mergeEdges[2], " - etb= ", edgeToBreak)))
    adriftReconnectionEdge <- mergeEdges[whichAdrift]
    rootedReconnectionEdge <- mergeEdges[!whichAdrift]
  }
  if(nearBrokenEdge[rootedReconnectionEdge] && nearBrokenEdge[adriftReconnectionEdge]) 
    return(TBRWarning(parent, child, "Selected mergeEdges will not change tree topology."))
  #### edgelabels(edge = edgeToBreak, bg="orange", cex=1.8)  #### DEBUGGING AID
  #### edgelabels(edge=adriftReconnectionEdge, bg="cyan")    #### DEBUGGING AID
  #### edgelabels(edge=rootedReconnectionEdge, bg="magenta") #### DEBUGGING AID
  
  ###Assert(edgesOnAdriftSegment[adriftReconnectionEdge])
  ###Assert(!edgesOnAdriftSegment[rootedReconnectionEdge])
  
  if (!nearBrokenEdge[adriftReconnectionEdge]) {
    edgesToInvert <- EdgeAncestry(adriftReconnectionEdge, parent, child, stopAt = edgeToBreak) & !brokenEdge
    if (any(edgesToInvert)) {
      tmp <- parent[edgesToInvert]
      parent[edgesToInvert] <- child[edgesToInvert]
      child[edgesToInvert] <- tmp
    }
    reconnectionSideEdges <- edgesToInvert
    reconnectionSideEdges[adriftReconnectionEdge] <- TRUE
    
    repurposedDaughterEdge <- brokenEdgeDaughters & reconnectionSideEdges
    spareDaughterEdge      <- brokenEdgeDaughters & !reconnectionSideEdges
    ###Assert(identical(sum(repurposedDaughterEdge), sum(spareDaughterEdge), 1))
    #### which(repurposedDaughterEdge)
    #### which(spareDaughterEdge)
    child[repurposedDaughterEdge] <- child[spareDaughterEdge]
    child[spareDaughterEdge] <- parent[adriftReconnectionEdge]
    ###Assert(parent[spareDaughterEdge] == brokenEdge.childNode)
    parent[adriftReconnectionEdge] <- child[edgeToBreak]
  }
  if (!nearBrokenEdge[rootedReconnectionEdge]) {
    parent[brokenEdgeSister] <- parent[brokenEdgeParent]
    parent[brokenEdgeParent] <- parent[rootedReconnectionEdge]
    parent[rootedReconnectionEdge] <- brokenEdge.parentNode
  }
  
  ###Assert(identical(unique(table(parent)), 2L))
  ###Assert(identical(unique(table(child)),  1L))
  return (RenumberEdges(parent, child))
}
