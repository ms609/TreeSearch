#' Consensus tree minimizing Clustering Information Distance
#'
#' Find a consensus tree that minimizes the mean Clustering Information
#' Distance (CID) to a set of input trees, using tree rearrangement
#' heuristics.
#'
#' Unlike the majority-rule consensus, which minimizes Robinson-Foulds
#' distance and can be highly unresolved when phylogenetic signal is low,
#' `CIDConsensus()` uses tree search (SPR, TBR, NNI, or the parsimony
#' ratchet) to find a more resolved tree that minimizes a finer-grained
#' information-theoretic distance to the input trees.
#'
#' The search has two phases:
#' 1. **Binary rearrangement** using SPR/TBR/NNI via [`TreeSearch()`] or
#'    [`Ratchet()`].
#' 2. **Collapse/resolve refinement** (when `collapse = TRUE`): greedily
#'    collapse internal edges whose removal reduces CID, then try resolving
#'    remaining polytomies.  This allows the result to be partially
#'    unresolved when that better represents the input trees.
#'
#' The ratchet variant perturbs the objective by resampling input trees
#' with replacement, analogous to character bootstrapping in parsimony.
#'
#' @param trees An object of class `multiPhylo`: the input trees.
#'   All trees must share the same tip labels.
#' @param metric Distance function with signature `f(tree1, tree2)` returning
#'   a numeric vector of distances.
#'   Default: [`ClusteringInfoDistance`][TreeDist::ClusteringInfoDistance].
#'   Other options include
#'   [`MutualClusteringInfo`][TreeDist::MutualClusteringInfo] or
#'   [`SharedPhylogeneticInfo`][TreeDist::SharedPhylogeneticInfo].
#' @param start A `phylo` tree to start the search from, or `NULL`
#'   (default) to start from the majority-rule consensus.
#'   Any non-binary starting tree is resolved with [`ape::multi2di()`]
#'   for the binary search phase.
#' @param method Character: search strategy.
#'   `"ratchet"` (default) uses the parsimony ratchet with TBR + SPR + NNI.
#'   `"spr"`, `"tbr"`, `"nni"` use the corresponding single swapper.
#' @param ratchIter Integer: maximum ratchet iterations (ignored unless
#'   `method = "ratchet"`).
#' @param ratchHits Integer: stop ratchet after hitting the same best score
#'   this many times.
#' @param searchIter Integer: maximum rearrangements per search iteration.
#' @param searchHits Integer: maximum times to hit best score before
#'   stopping a search iteration.
#' @param collapse Logical: if `TRUE` (default), run a collapse/resolve
#'   refinement phase after the binary search.  This can produce a
#'   non-binary result when collapsing a split reduces the mean CID.
#' @param verbosity Integer controlling console output (0 = silent).
#' @param \dots Additional arguments passed to [`Ratchet()`] or
#'   [`TreeSearch()`].
#'
#' @return A tree of class `phylo` with attributes:
#' - `"score"`: the mean distance to input trees under `metric`.
#' - `"hits"`: the number of times this score was found.
#'
#' @examples
#' library(TreeTools)
#' # Generate some trees
#' trees <- as.phylo(1:30, nTip = 12)
#'
#' # Quick search (increase ratchIter for real analyses)
#' result <- CIDConsensus(trees, ratchIter = 2, searchHits = 5,
#'                        verbosity = 0)
#' plot(result)
#' attr(result, "score")
#'
#' # Compare with majority-rule consensus
#' mr <- Consensus(trees, p = 0.5)
#' mean(TreeDist::ClusteringInfoDistance(mr, trees))
#' mean(TreeDist::ClusteringInfoDistance(result, trees))
#'
#' @seealso [Ratchet()] for the underlying search algorithm.
#'
#' @references
#' \insertRef{Smith2020}{TreeSearch}
#'
#' @importFrom TreeDist ClusteringInfoDistance ClusteringEntropy
#'   MutualClusteringInfoSplits
#' @importFrom TreeTools as.Splits Consensus RenumberEdges
#' @importFrom ape multi2di
#' @family custom search functions
#' @export
CIDConsensus <- function(trees,
                         metric = ClusteringInfoDistance,
                         start = NULL,
                         method = c("ratchet", "spr", "tbr", "nni"),
                         ratchIter = 100L,
                         ratchHits = 10L,
                         searchIter = 500L,
                         searchHits = 20L,
                         collapse = TRUE,
                         verbosity = 1L,
                         ...) {
  method <- match.arg(method)
  if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be an object of class 'multiPhylo'.")
  }
  if (length(trees) < 2L) {
    stop("Need at least 2 input trees.")
  }

  tipLabels <- trees[[1]][["tip.label"]]

  # Validate consistent tip labels
  nTip <- length(tipLabels)
  for (i in seq_along(trees)) {
    if (length(trees[[i]][["tip.label"]]) != nTip) {
      stop("All input trees must have the same number of tips. ",
           "Tree ", i, " has ", length(trees[[i]][["tip.label"]]),
           " tips; expected ", nTip, ".")
    }
  }

  # Build starting tree
  if (is.null(start)) {
    start <- Consensus(trees, p = 0.5)
  }
  start <- multi2di(start)

  # Prepare CID dataset (environment for reference semantics)
  cidData <- .MakeCIDData(trees, metric, tipLabels)

  if (method == "ratchet") {
    result <- Ratchet(start, cidData,
                      InitializeData = identity,
                      CleanUpData    = .NoOp,
                      TreeScorer     = .CIDScorer,
                      Bootstrapper   = .CIDBootstrap,
                      swappers = list(RootedTBRSwap, RootedSPRSwap,
                                      RootedNNISwap),
                      ratchIter = ratchIter,
                      ratchHits = ratchHits,
                      searchIter = searchIter,
                      searchHits = searchHits,
                      verbosity = verbosity,
                      ...)
  } else {
    edgeSwapper <- switch(method,
      spr = RootedSPRSwap,
      tbr = RootedTBRSwap,
      nni = RootedNNISwap
    )
    result <- TreeSearch(start, cidData,
                         InitializeData = identity,
                         CleanUpData    = .NoOp,
                         TreeScorer     = .CIDScorer,
                         EdgeSwapper    = edgeSwapper,
                         maxIter = searchIter,
                         maxHits = searchHits,
                         verbosity = verbosity,
                         ...)
  }

  # Phase 2: Collapse/resolve refinement
  if (collapse) {
    result <- .CollapseRefine(result, cidData, verbosity)
  }

  result
}


.NoOp <- function(x) invisible(NULL)


# Build CID dataset.
# Uses an environment for reference semantics (CIDBootstrap needs to swap
# the tree list temporarily). S3 class "cidData" with a names() method
# so that TreeSearch/Ratchet can call names(dataset) to get tip labels.
#
# For CID (the default metric), precomputes input tree splits and
# clustering entropies to avoid redundant O(N) work per candidate.
.MakeCIDData <- function(trees, metric, tipLabels) {
  env <- new.env(parent = emptyenv())
  env$trees <- trees
  env$metric <- metric
  env$tipLabels <- tipLabels
  env$nTip <- length(tipLabels)

  # Precompute splits and entropies for the default CID path
  isCID <- identical(metric, ClusteringInfoDistance)
  env$isCID <- isCID
  if (isCID) {
    inputSplits <- lapply(trees, as.Splits, tipLabels)
    env$inputCE <- vapply(inputSplits, ClusteringEntropy, double(1))
    env$meanInputCE <- mean(env$inputCE)
    # Store raw (unclass'd) split matrices for direct C++ access
    env$inputSplitsRaw <- lapply(inputSplits, unclass)
  }

  class(env) <- "cidData"
  env
}

#' @export
names.cidData <- function(x) x$tipLabels


# CID-based TreeScorer: mean distance from candidate to all input trees.
# For ClusteringInfoDistance, uses precomputed raw splits for ~9x speedup
# over the naive ClusteringInfoDistance(phylo, multiPhylo) approach.
.CIDScorer <- function(parent, child, dataset, ...) {
  if (dataset$isCID) {
    .CIDScoreFast(parent, child, dataset)
  } else {
    candidate <- .EdgeListToPhylo(parent, child, dataset$tipLabels)
    mean(dataset$metric(candidate, dataset$trees))
  }
}


# Fast CID scorer using precomputed input splits.
# CID(cand, tree_i) = CE(cand) + CE(tree_i) - 2*MCI(cand, tree_i)
# mean(CID) = CE(cand) + mean(CE(inputs)) - 2*mean(MCI)
#
# Uses a for loop instead of vapply to avoid function-call overhead per
# input tree, and operates on unclass'd raw split matrices.
.CIDScoreFast <- function(parent, child, dataset) {
  nTip <- dataset$nTip
  candidate <- .EdgeListToPhylo(parent, child, dataset$tipLabels)
  candSp <- as.Splits(candidate, dataset$tipLabels)
  candCE <- ClusteringEntropy(candSp)
  inputSplitsRaw <- dataset$inputSplitsRaw
  nTree <- length(inputSplitsRaw)
  candRaw <- unclass(candSp)
  mciSum <- 0
  for (i in seq_len(nTree)) {
    mciSum <- mciSum + MutualClusteringInfoSplits(candSp, inputSplitsRaw[[i]],
                                                  nTip)
  }
  candCE + dataset$meanInputCE - 2 * mciSum / nTree
}


# Convert parent/child edge vectors to a phylo object.
.EdgeListToPhylo <- function(parent, child, tipLabels) {
  nTip <- length(tipLabels)
  structure(list(
    edge = cbind(parent, child),
    tip.label = tipLabels,
    Nnode = length(unique(parent[parent > nTip]))
  ), class = "phylo")
}


# CID bootstrapper: resample input trees with replacement, then search.
# Also resamples precomputed splits/CEs for the fast CID path.
.CIDBootstrap <- function(edgeList, cidData,
                          EdgeSwapper = RootedNNISwap,
                          maxIter, maxHits,
                          verbosity = 1L,
                          stopAtPeak = FALSE,
                          stopAtPlateau = 0L, ...) {
  origTrees <- cidData$trees
  nTree <- length(origTrees)
  idx <- sample.int(nTree, replace = TRUE)
  cidData$trees <- origTrees[idx]

  # Also resample precomputed data for fast CID path
  if (cidData$isCID) {
    origSplitsRaw <- cidData$inputSplitsRaw
    origCE <- cidData$inputCE
    origMeanCE <- cidData$meanInputCE
    cidData$inputSplitsRaw <- origSplitsRaw[idx]
    cidData$inputCE <- origCE[idx]
    cidData$meanInputCE <- mean(origCE[idx])
    on.exit({
      cidData$trees <- origTrees
      cidData$inputSplitsRaw <- origSplitsRaw
      cidData$inputCE <- origCE
      cidData$meanInputCE <- origMeanCE
    })
  } else {
    on.exit(cidData$trees <- origTrees)
  }

  res <- EdgeListSearch(edgeList[1:2], cidData,
                        TreeScorer = .CIDScorer,
                        EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        verbosity = verbosity,
                        stopAtPeak = stopAtPeak,
                        stopAtPlateau = stopAtPlateau, ...)
  res[1:2]
}


# Greedy collapse/resolve refinement.
# Iteratively tries collapsing each internal edge; accepts if score improves.
# Then tries resolving each polytomy; accepts if score improves.
# Repeats until no improvement.
.CollapseRefine <- function(tree, cidData, verbosity = 0L) {
  tipLabels <- cidData$tipLabels
  nTip <- cidData$nTip
  edge <- tree[["edge"]]
  parent <- edge[, 1]
  child <- edge[, 2]

  bestScore <- .CIDScorer(parent, child, cidData)

  if (verbosity > 0L) {
    message("  - Collapse/resolve refinement. Starting score: ",
            signif(bestScore, 6))
  }

  improved <- TRUE
  while (improved) {
    improved <- FALSE

    # --- Collapse pass: try removing each internal edge ---
    internalEdges <- which(child > nTip)
    if (length(internalEdges) > 0L) {
      for (i in rev(seq_along(internalEdges))) {
        edgeIdx <- internalEdges[i]
        candidate <- .CollapseSpecificEdge(parent, child, edgeIdx, nTip)
        candScore <- .CIDScorer(candidate[[1]], candidate[[2]], cidData)
        if (candScore < bestScore - sqrt(.Machine[["double.eps"]])) {
          parent <- candidate[[1]]
          child <- candidate[[2]]
          bestScore <- candScore
          improved <- TRUE
          if (verbosity > 1L) {
            message("    * Collapsed edge → score ", signif(bestScore, 6))
          }
          # Recompute internal edges since topology changed
          break
        }
      }
      if (improved) next
    }

    # --- Resolve pass: try resolving each polytomy ---
    degrees <- tabulate(parent, nbins = max(parent))
    polyNodes <- which(degrees > 2L)
    if (length(polyNodes) > 0L) {
      for (node in polyNodes) {
        childEdges <- which(parent == node)
        nChildren <- length(childEdges)
        if (nChildren <= 2L) next

        # Try all pairs of children as candidates for a new clade
        bestResolve <- NULL
        bestResolveScore <- bestScore
        for (a in 1:(nChildren - 1L)) {
          for (b in (a + 1L):nChildren) {
            candidate <- .ResolveSpecificPair(
              parent, child, node,
              childEdges[c(a, b)], nTip
            )
            candScore <- .CIDScorer(candidate[[1]], candidate[[2]], cidData)
            if (candScore < bestResolveScore - sqrt(.Machine[["double.eps"]])) {
              bestResolve <- candidate
              bestResolveScore <- candScore
            }
          }
        }
        if (!is.null(bestResolve)) {
          parent <- bestResolve[[1]]
          child <- bestResolve[[2]]
          bestScore <- bestResolveScore
          improved <- TRUE
          if (verbosity > 1L) {
            message("    * Resolved polytomy → score ", signif(bestScore, 6))
          }
          break
        }
      }
    }
  }

  if (verbosity > 0L) {
    message("  - Collapse/resolve complete. Final score: ",
            signif(bestScore, 6))
  }

  result <- .EdgeListToPhylo(parent, child, tipLabels)
  attr(result, "score") <- bestScore
  result
}


# Collapse a specific edge (by index), returning renumbered parent/child.
.CollapseSpecificEdge <- function(parent, child, edgeIdx, nTip) {
  nEdge <- length(parent)
  collapseParent <- parent[edgeIdx]
  collapseChild <- child[edgeIdx]

  # Reparent all children of collapseChild to collapseParent
  childOfCollapsed <- which(parent == collapseChild)
  parent[childOfCollapsed] <- collapseParent

  # Remove the collapsed edge
  keep <- seq_len(nEdge) != edgeIdx
  parent <- parent[keep]
  child <- child[keep]

  .RenumberNodes(parent, child, nTip)
}


# Resolve a specific pair of child edges under a polytomy node.
.ResolveSpecificPair <- function(parent, child, node, moveEdges, nTip) {
  newNode <- max(c(parent, child)) + 1L

  # New edge: node -> newNode
  parent <- c(parent, node)
  child <- c(child, newNode)

  # Reparent selected children
  parent[moveEdges] <- newNode

  list(parent, child)
}


# Collapse a random internal edge, creating a polytomy.
.CollapseEdge <- function(parent, child, nTip) {
  nEdge <- length(parent)
  # Internal edges: both parent and child are internal nodes
  internalEdges <- which(child > nTip)

  if (length(internalEdges) == 0L) {
    # Star tree — nothing to collapse
    return(list(parent, child))
  }

  # Pick a random internal edge to collapse
  edgeIdx <- internalEdges[sample.int(length(internalEdges), 1L)]
  collapseParent <- parent[edgeIdx]
  collapseChild <- child[edgeIdx]

  # Reparent all children of collapseChild to collapseParent
  childOfCollapsed <- which(parent == collapseChild)
  parent[childOfCollapsed] <- collapseParent

  # Remove the collapsed edge
  keep <- seq_len(nEdge) != edgeIdx
  parent <- parent[keep]
  child <- child[keep]

  # Renumber internal nodes to fill the gap
  .RenumberNodes(parent, child, nTip)
}


# Resolve a random polytomy by inserting a new internal node.
.ResolveNode <- function(parent, child, nTip) {
  # Find polytomy nodes (> 2 children)
  degrees <- tabulate(parent)
  polyNodes <- which(degrees > 2L)

  if (length(polyNodes) == 0L) {
    return(list(parent, child))
  }

  # Pick a random polytomy
  node <- polyNodes[sample.int(length(polyNodes), 1L)]
  childEdges <- which(parent == node)
  nChildren <- length(childEdges)

  # Pick 2+ children to move to a new internal node
  # (pick exactly 2 for a single-step resolution)
  nToMove <- 2L
  if (nChildren <= nToMove) {
    # Can't resolve a node with <= 2 children
    return(list(parent, child))
  }
  moveEdges <- childEdges[sample.int(nChildren, nToMove)]

  # Insert new node
  nNode <- max(parent)
  newNode <- nNode + 1L

  # New edge: node -> newNode
  parent <- c(parent, node)
  child <- c(child, newNode)

  # Reparent selected children
  parent[moveEdges] <- newNode

  list(parent, child)
}


# EdgeSwapper: collapse a random internal edge.
.CollapseSwap <- function(parent, child, nTip = min(parent) - 1L, ...) {
  .CollapseEdge(parent, child, nTip)
}


# EdgeSwapper: resolve a random polytomy.
.ResolveSwap <- function(parent, child, nTip = min(parent) - 1L, ...) {
  .ResolveNode(parent, child, nTip)
}


# Renumber internal nodes to be contiguous from nTip+1.
.RenumberNodes <- function(parent, child, nTip) {
  allNodes <- sort(unique(c(parent, child)))
  internalNodes <- allNodes[allNodes > nTip]
  # Map old node numbers to new contiguous range
  newNumbers <- seq_along(internalNodes) + nTip
  nodeMap <- integer(max(internalNodes))
  # Tips map to themselves
  nodeMap[seq_len(nTip)] <- seq_len(nTip)
  nodeMap[internalNodes] <- newNumbers

  list(nodeMap[parent], nodeMap[child])
}
