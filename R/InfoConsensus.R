#' Information-theoretic consensus tree
#'
#' Find a consensus tree that maximizes the mean Mutual Clustering
#' Information (MCI) with a set of input trees, using a driven search with
#' TBR, ratchet, drift, sectorial search, and tree fusing.
#'
#' Unlike the majority-rule consensus, which minimizes Robinson-Foulds
#' distance and can be highly unresolved when phylogenetic signal is low,
#' `InfoConsensus()` finds a more resolved tree that maximizes a finer-grained
#' information-theoretic measure of agreement with the input trees.
#'
#' The search uses MRP (Matrix Representation with Parsimony) characters
#' for fast incremental screening during TBR, with MCI verification for
#' move acceptance.  This provides the full power of the driven search
#' pipeline (ratchet, drift, sectorial search, tree fusing, multi-replicate
#' parallelism) with MCI scoring.
#'
#' The search proceeds in up to three phases:
#' 1. **Driven search** using the C++ engine with MCI scoring.
#' 2. **Collapse/resolve refinement** (when `collapse = TRUE`): greedily
#'    collapse internal edges whose removal improves the score, then try
#'    resolving remaining polytomies.
#' 3. **Rogue taxon dropping** (when `neverDrop != TRUE`): iteratively
#'    identify and remove taxa whose absence improves consensus quality,
#'    then attempt to restore previously dropped taxa.
#'
#' @param trees An object of class `multiPhylo`: the input trees.
#'   All trees must share the same tip labels.
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates.  Each replicate starts from a different random Wagner tree.
#' @param targetHits Integer: stop after finding the best score this many
#'   times independently.
#' @param maxSeconds Numeric: timeout in seconds (0 = no timeout).
#' @param nThreads Integer: number of threads for inter-replicate parallelism.
#'   Defaults to `getOption("mc.cores", 1L)`.
#' @param collapse Logical: if `TRUE` (default), run a collapse/resolve
#'   refinement phase after the binary search.  This can produce a
#'   non-binary result when collapsing a split improves the mean MCI.
#' @param neverDrop Controls rogue taxon dropping (Phase 3).
#'   `TRUE` disables dropping entirely.
#'   `FALSE` (default) allows all taxa to be dropped if doing so improves
#'   the consensus.
#'   A character or integer vector specifies tips that must never be dropped;
#'   all others are candidates.
#' @param maxDrop Integer: maximum number of tips that may be dropped during
#'   rogue screening.  Default `ceiling(nTip / 10)` (10 percent of tips).
#' @param control A [`SearchControl()`] object for expert tuning of the
#'   driven search strategy.
#' @param screeningK Numeric: implied-weights concavity constant for MRP
#'   character screening during TBR.  The search uses MRP (Matrix
#'   Representation with Parsimony) characters as a fast proxy for MCI
#'   scoring; `screeningK` controls how these characters are weighted.
#'   Default `7` (IW with `k = 7`), which empirically maximizes rank
#'   correlation with MCI scores.  Use `Inf` for equal-weight screening.
#' @param screeningTolerance Numeric (>= 0): controls how generously
#'   candidate moves are sent to full MCI evaluation.  A value of `0`
#'   (default) only evaluates the single best MRP-screened candidate per
#'   clip.  Values > 0 relax the screening threshold, allowing candidates
#'   whose MRP score exceeds the current best by up to this fraction (e.g.,
#'   `0.02` = $2%$ tolerance).  Higher values improve search quality at the
#'   cost of more MCI evaluations per step.
#' @param verbosity Integer controlling console output (0 = silent).
#'
#' @return A tree of class `phylo` with attributes:
#' - `"score"`: mean MCI between the consensus and the input trees
#'   (higher is better).
#' - `"hits"`: the number of times this score was found.
#' - `"droppedTips"`: character vector of dropped taxa (if any), or `NULL`.
#'
#' @examples
#' library(TreeTools)
#' # Generate some trees
#' trees <- as.phylo(1:30, nTip = 12)
#'
#' # Quick search
#' result <- InfoConsensus(trees, maxReplicates = 3L, targetHits = 2L,
#'                         neverDrop = TRUE, verbosity = 0)
#' plot(result)
#' attr(result, "score")
#'
#' # Compare with majority-rule consensus
#' mr <- Consensus(trees, p = 0.5)
#' mean(TreeDist::ClusteringInfoDistance(mr, trees))
#' mean(TreeDist::ClusteringInfoDistance(result, trees))
#'
#' @seealso [MaximizeParsimony()] uses the same driven search engine for
#'   parsimony.
#' @seealso [TreeDist::TransferDistance()] and [Quartet::QuartetConsensus()]
#'   for alternative consensus methods.
#' % TODO: add \insertCite reference for Smith 2026 when available in TreeDist
#'
#' @references
#' \insertRef{Smith2020}{TreeSearch}
#'
#' @importFrom TreeDist ClusteringInfoDistance ClusteringEntropy
#'   MutualClusteringInfoSplits
#' @importFrom TreeTools as.Splits Consensus MakeTreeBinary NTip RenumberEdges
#' @family custom search functions
#' @export
InfoConsensus <- function(trees,
                          maxReplicates = 100L,
                          targetHits = 10L,
                          maxSeconds = 0,
                          nThreads = getOption("mc.cores", 1L),
                          collapse = TRUE,
                          neverDrop = FALSE,
                          maxDrop = ceiling(NTip(trees[[1]]) / 10),
                          control = SearchControl(),
                          screeningK = 7,
                          screeningTolerance = 0,
                          verbosity = 1L) {
  if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be an object of class 'multiPhylo'.")
  }
  if (length(trees) < 2L) {
    stop("Need at least 2 input trees.")
  }
  
  tipLabels <- trees[[1]][["tip.label"]]
  nTip <- length(tipLabels)
  
  # Validate consistent tip labels
  for (i in seq_along(trees)) {
    if (length(trees[[i]][["tip.label"]]) != nTip) {
      stop("All input trees must have the same number of tips. ",
           "Tree ", i, " has ", length(trees[[i]][["tip.label"]]),
           " tips; expected ", nTip, ".")
    }
  }
  
  # Phase 1: C++ driven search
  result <- .CIDDrivenSearch(trees, tipLabels, nTip,
                              maxReplicates, targetHits, maxSeconds,
                              nThreads, control,
                              screeningK, screeningTolerance,
                              verbosity)
  
  # Phase 2: Collapse/resolve refinement
  if (collapse) {
    cidData <- .MakeCIDData(trees, tipLabels)
    result <- .CollapseRefine(result, cidData, verbosity)
  }
  
  # Phase 3: rogue taxon dropping
  if (!isTRUE(neverDrop)) {
    cidData <- if (exists("cidData", inherits = FALSE)) {
      cidData
    } else {
      .MakeCIDData(trees, tipLabels)
    }
    result <- .RogueRefine(result, cidData, neverDrop, maxDrop,
                           "ratchet",       # method for re-optimization
                           5L, 3L,          # light ratchet for re-opt
                           100L, 10L, collapse,
                           verbosity)
  }
  
  # Convert internal score (negated MCI) to user-facing positive MCI
  internalScore <- attr(result, "score")
  if (!is.null(internalScore)) {
    attr(result, "score") <- -internalScore
  }
  
  result
}





.NoOp <- function(x) invisible(NULL)

# Null-coalesce (base R %||% requires R >= 4.4)
.NullOr <- function(x, default) if (is.null(x)) default else x


# Light re-optimization via the C++ driven search.
# Used by .RogueRefine() after dropping a rogue taxon.
# Arguments mirror the old TopologySearch interface for backward compatibility
# with .RogueRefine(); the actual search delegates to ts_cid_consensus.
.TopologySearch <- function(tree, cidData, method,
                            ratchIter, ratchHits,
                            searchIter, searchHits, collapse,
                            verbosity) {
  tipLabels <- cidData$tipLabels
  nTip <- cidData$nTip

  splitMats <- cidData$inputSplitsRaw

  # C++ engine expects a binary tree; resolve any polytomies
  startTree <- if (ape::is.binary(tree)) tree else MakeTreeBinary(tree)
  startEdge <- startTree[["edge"]]

  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = nTip,
    normalize = FALSE,
    maxReplicates = max(1L, as.integer(ratchIter)),
    targetHits = max(1L, as.integer(ratchHits)),
    verbosity = max(0L, as.integer(verbosity) - 1L),
    nThreads = 1L,
    startEdge = startEdge
  )

  if (result[["pool_size"]] == 0L) {
    attr(tree, "score") <- .ScoreTree(tree, cidData)
    return(tree)
  }

  bestEdge <- result[["trees"]][[1]]
  bestTree <- .EdgeListToPhylo(bestEdge[, 1], bestEdge[, 2], tipLabels)
  attr(bestTree, "score") <- result[["best_score"]]

  if (collapse) {
    bestTree <- .CollapseRefine(bestTree, cidData,
                                verbosity = max(0L, verbosity - 1L))
  }

  bestTree
}


# Phase 1: C++ driven search with MCI scoring.
# Converts input trees to split matrices and calls the C++ engine.
# Returns tree with attr("score") = positive MCI (higher = better).
.CIDDrivenSearch <- function(trees, tipLabels, nTip,
                              maxReplicates, targetHits,
                              maxSeconds, nThreads, control,
                              screeningK, screeningTolerance,
                              verbosity) {
  # Convert input trees to split matrices (RawMatrix format)
  splitMats <- lapply(trees, function(tr) {
    unclass(as.Splits(tr, tipLabels))
  })
  
  # Extract SearchControl parameters
  ctrl <- control
  
  # Sectors only benefit trees large enough for meaningful partitioning.
  # CID scoring is expensive per evaluation, so skip sectors on small trees
  # where TBR alone converges quickly.
  useSectors <- nTip >= 20L
  
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = nTip,
    normalize = FALSE,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    tbrMaxHits = .NullOr(ctrl[["tbrMaxHits"]], 1L),
    ratchetCycles = .NullOr(ctrl[["ratchetCycles"]], 10L),
    ratchetPerturbProb = .NullOr(ctrl[["ratchetPerturbProb"]], 0.04),
    ratchetPerturbMode = .NullOr(ctrl[["ratchetPerturbMode"]], 0L),
    ratchetAdaptive = .NullOr(ctrl[["ratchetAdaptive"]], FALSE),
    driftCycles = .NullOr(ctrl[["driftCycles"]], 6L),
    driftAfdLimit = .NullOr(ctrl[["driftAfdLimit"]], 3L),
    driftRfdLimit = .NullOr(ctrl[["driftRfdLimit"]], 0.1),
    xssRounds = .NullOr(ctrl[["xssRounds"]], if (useSectors) 3L else 0L),
    xssPartitions = .NullOr(ctrl[["xssPartitions"]], 4L),
    rssRounds = .NullOr(ctrl[["rssRounds"]], if (useSectors) 3L else 0L),
    cssRounds = .NullOr(ctrl[["cssRounds"]], if (useSectors) 2L else 0L),
    cssPartitions = .NullOr(ctrl[["cssPartitions"]], 4L),
    sectorMinSize = .NullOr(ctrl[["sectorMinSize"]], 6L),
    sectorMaxSize = .NullOr(ctrl[["sectorMaxSize"]], 50L),
    fuseInterval = .NullOr(ctrl[["fuseInterval"]], 3L),
    fuseAcceptEqual = .NullOr(ctrl[["fuseAcceptEqual"]], FALSE),
    poolMaxSize = 100L,
    poolSuboptimal = 0.0,
    maxSeconds = maxSeconds,
    verbosity = verbosity,
    tabuSize = .NullOr(ctrl[["tabuSize"]], 100L),
    wagnerStarts = .NullOr(ctrl[["wagnerStarts"]], 1L),
    nThreads = nThreads,
    screeningK = screeningK,
    screeningTolerance = screeningTolerance
  )
  
  # Convert best tree from edge matrix to phylo
  if (result[["pool_size"]] == 0L) {
    warning("MCI search found no trees; returning majority-rule consensus.",
            call. = FALSE)
    tree <- Consensus(trees, p = 0.5)
    # Internal: 0 is worst possible negated MCI
    attr(tree, "score") <- 0
    attr(tree, "hits") <- 0L
    return(tree)
  }
  
  bestEdge <- result[["trees"]][[1]]
  tree <- .EdgeListToPhylo(bestEdge[, 1], bestEdge[, 2], tipLabels)
  # Internal score: negated MCI from C++ (lower = better)
  attr(tree, "score") <- result[["best_score"]]
  attr(tree, "hits") <- result[["hits_to_best"]]
  tree
}


# Build CID dataset.
# Uses an environment for reference semantics. S3 class "cidData" with a
# names() method so that TreeSearch/Ratchet can call names(dataset) to get
# tip labels.  Precomputes input tree splits and clustering entropies.
.MakeCIDData <- function(trees, tipLabels) {
  env <- new.env(parent = emptyenv())
  env$trees <- trees
  env$tipLabels <- tipLabels
  env$nTip <- length(tipLabels)
  
  inputSplits <- lapply(trees, as.Splits, tipLabels)
  env$inputCE <- vapply(inputSplits, ClusteringEntropy, double(1))
  env$meanInputCE <- mean(env$inputCE)
  # Store raw (unclass'd) split matrices for direct C++ access
  env$inputSplitsRaw <- lapply(inputSplits, unclass)
  
  class(env) <- "cidData"
  env
}

#' @export
names.cidData <- function(x) x$tipLabels


# Score candidate tree against all input trees using MCI.
# Returns negated mean MCI (lower = better).
.CIDScorer <- function(parent, child, dataset) {
  .CIDScoreFast(parent, child, dataset)
}


# Fast MCI scorer using precomputed input splits.
# Returns negated mean MCI (lower = better), consistent with C++ cid_score().
.CIDScoreFast <- function(parent, child, dataset) {
  nTip <- dataset$nTip
  candidate <- .EdgeListToPhylo(parent, child, dataset$tipLabels)
  candSp <- as.Splits(candidate, dataset$tipLabels)
  inputSplitsRaw <- dataset$inputSplitsRaw
  nTree <- length(inputSplitsRaw)
  
  mciSum <- 0
  for (i in seq_len(nTree)) {
    mciSum <- mciSum + MutualClusteringInfoSplits(candSp, inputSplitsRaw[[i]],
                                                  nTip)
  }
  -mciSum / nTree
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
  
  # Also resample precomputed splits/entropies
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
    message("  - Collapse/resolve refinement. Starting MCI: ",
            signif(-bestScore, 6))
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
            message("    * Collapsed edge -> MCI ", signif(-bestScore, 6))
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
            message("    * Resolved polytomy -> MCI ", signif(-bestScore, 6))
          }
          break
        }
      }
    }
  }
  
  if (verbosity > 0L) {
    message("  - Collapse/resolve complete. Final MCI: ",
            signif(-bestScore, 6))
  }
  
  result <- .EdgeListToPhylo(parent, child, tipLabels)
  attr(result, "score") <- bestScore
  result
}



# --- Phase 3: Rogue taxon dropping and restoration -------------------------

# Greedy rogue dropping followed by greedy restoration.
.RogueRefine <- function(tree, cidData, neverDrop, maxDrop,
                         method, ratchIter, ratchHits,
                         searchIter, searchHits, collapse,
                         verbosity) {
  originalTrees <- cidData$trees
  allTipLabels <- cidData$tipLabels
  bestScore <- .ScoreTree(tree, cidData)
  currentTips <- tree[["tip.label"]]
  nTip <- length(currentTips)
  if (nTip < 5L) return(tree)
  protected <- if (isFALSE(neverDrop)) {
    character(0)
  } else if (is.character(neverDrop)) {
    bad <- setdiff(neverDrop, allTipLabels)
    if (length(bad)) {
      warning("neverDrop tips not found in trees: ",
              paste(bad, collapse = ", "), call. = FALSE)
    }
    intersect(neverDrop, allTipLabels)
  } else if (is.numeric(neverDrop)) {
    allTipLabels[as.integer(neverDrop)]
  } else {
    character(0)
  }
  droppedTips <- character(0)
  maxDrop <- as.integer(min(maxDrop, nTip - max(4L, length(protected) + 1L)))
  if (verbosity > 0L) {
    message("  - Rogue screening phase. Starting MCI: ",
            signif(-bestScore, 6), " (", nTip, " tips, maxDrop = ",
            maxDrop, ")")
  }
  improved <- TRUE
  while (improved &&
         length(droppedTips) < maxDrop &&
         length(currentTips) > max(4L, length(protected) + 1L)) {
    improved <- FALSE
    droppable <- setdiff(currentTips, protected)
    if (length(droppable) == 0L) break
    prescreenScores <- .PrescreenMarginalNID(
      tree, cidData, droppable, originalTrees, allTipLabels, droppedTips
    )
    candidates <- names(sort(prescreenScores))
    for (tip in candidates) {
      if (prescreenScores[tip] >= bestScore - sqrt(.Machine[["double.eps"]])) {
        break
      }
      reducedTips <- setdiff(currentTips, tip)
      allDropped <- c(droppedTips, tip)
      reducedInputTrees <- .PruneTrees(originalTrees, allDropped)
      reducedCidData <- .MakeCIDData(reducedInputTrees, reducedTips)
      reducedTree <- DropTip(tree, tip)
      nReducedTips <- length(reducedTips)
      if (nReducedTips < 5L) {
        reoptResult <- reducedTree
        attr(reoptResult, "score") <- .ScoreTree(reducedTree, reducedCidData)
      } else if (nReducedTips >= 8L && method == "ratchet") {
        reoptResult <- .TopologySearch(
          reducedTree, reducedCidData, method,
          max(1L, ratchIter %/% 2L), ratchHits,
          searchIter, searchHits, collapse,
          max(0L, verbosity - 1L))
      } else {
        reoptResult <- .TopologySearch(
          reducedTree, reducedCidData, "nni",
          1L, 1L, searchIter, searchHits, collapse,
          max(0L, verbosity - 1L))
      }
      newScore <- attr(reoptResult, "score")
      if (is.null(newScore)) newScore <- .ScoreTree(reoptResult, reducedCidData)
      if (newScore < bestScore - sqrt(.Machine[["double.eps"]])) {
        tree <- reoptResult
        bestScore <- newScore
        currentTips <- reducedTips
        cidData <- reducedCidData
        droppedTips <- allDropped
        improved <- TRUE
        if (verbosity > 0L) {
          message("    * Dropped '", tip, "' -> MCI ",
                  signif(-bestScore, 6),
                  " (", length(currentTips), " tips)")
        }
        break
      }
    }
  }
  if (length(droppedTips) > 0L) {
    if (verbosity > 0L) {
      message("  - Restore phase: trying to re-insert ",
              length(droppedTips), " dropped tip(s)")
    }
    restoredAny <- TRUE
    while (restoredAny) {
      restoredAny <- FALSE
      for (idx in rev(seq_along(droppedTips))) {
        tip <- droppedTips[idx]
        insertion <- .BestInsertion(
          tree, tip, originalTrees,
          droppedTips[-idx]
        )
        if (insertion$score < bestScore - sqrt(.Machine[["double.eps"]])) {
          tree <- insertion$tree
          bestScore <- insertion$score
          cidData <- insertion$cidData
          currentTips <- tree[["tip.label"]]
          droppedTips <- droppedTips[-idx]
          restoredAny <- TRUE
          if (verbosity > 0L) {
            message("    * Restored '", tip, "' -> MCI ",
                    signif(-bestScore, 6),
                    " (", length(currentTips), " tips)")
          }
          break
        }
      }
    }
  }
  if (verbosity > 0L) {
    message("  - Rogue screening complete. Final MCI: ",
            signif(-bestScore, 6),
            if (length(droppedTips)) paste0(
              " (dropped: ", paste(droppedTips, collapse = ", "), ")"
            ) else " (no tips dropped)")
  }
  attr(tree, "score") <- bestScore
  attr(tree, "droppedTips") <- if (length(droppedTips)) droppedTips else NULL
  tree
}


.PrescreenMarginalNID <- function(tree, cidData, droppable,
                                  originalTrees, allTipLabels,
                                  alreadyDropped) {
  vapply(droppable, function(tip) {
    reducedTips <- setdiff(tree[["tip.label"]], tip)
    allDropped <- c(alreadyDropped, tip)
    reducedInputTrees <- .PruneTrees(originalTrees, allDropped)
    reducedCidData <- .MakeCIDData(reducedInputTrees, reducedTips)
    reducedTree <- DropTip(tree, tip)
    .ScoreTree(reducedTree, reducedCidData)
  }, double(1))
}


.PruneTrees <- function(trees, tipsToDrop) {
  if (length(tipsToDrop) == 0L) return(trees)
  pruned <- lapply(trees, DropTip, tip = tipsToDrop)
  class(pruned) <- "multiPhylo"
  pruned
}


.ScoreTree <- function(tree, cidData) {
  edge <- tree[["edge"]]
  .CIDScorer(edge[, 1], edge[, 2], cidData)
}


.EdgeListToPhylo <- function(parent, child, tipLabels) {
  nTip <- length(tipLabels)
  result <- structure(list(
    edge = cbind(parent, child),
    tip.label = tipLabels,
    Nnode = length(unique(parent[parent > nTip]))
  ), class = "phylo")
  reorder(result)
}


.BestInsertion <- function(tree, tipLabel, originalTrees,
                           otherDropped) {
  currentTips <- c(tree[["tip.label"]], tipLabel)
  if (length(otherDropped) > 0L) {
    inputTrees <- .PruneTrees(originalTrees, otherDropped)
  } else {
    inputTrees <- originalTrees
  }
  testCidData <- .MakeCIDData(inputTrees, currentTips)
  nEdge <- nrow(tree[["edge"]])
  bestScore <- Inf
  bestTree <- NULL
  for (i in seq_len(nEdge)) {
    candidate <- .InsertTipAtEdge(tree, tipLabel, i)
    score <- .ScoreTree(candidate, testCidData)
    if (score < bestScore) {
      bestScore <- score
      bestTree <- candidate
    }
  }
  list(tree = bestTree, score = bestScore, cidData = testCidData)
}


.InsertTipAtEdge <- function(tree, tipLabel, edgeIdx) {
  edge <- tree[["edge"]]
  tipLabels <- tree[["tip.label"]]
  nTip <- length(tipLabels)
  nNode <- tree[["Nnode"]]
  pNode <- edge[edgeIdx, 1]
  cNode <- edge[edgeIdx, 2]
  newTipIdx <- nTip + 1L
  newEdge <- edge
  internals <- newEdge > nTip
  newEdge[internals] <- newEdge[internals] + 1L
  newPNode <- pNode + 1L
  newCNode <- if (cNode > nTip) cNode + 1L else cNode
  newInternalIdx <- as.integer(nTip + nNode + 2L)
  newEdge[edgeIdx, 2] <- newInternalIdx
  newEdge <- rbind(newEdge,
                   c(newInternalIdx, newCNode),
                   c(newInternalIdx, newTipIdx))
  result <- structure(list(
    edge = newEdge,
    tip.label = c(tipLabels, tipLabel),
    Nnode = nNode + 1L
  ), class = "phylo")
  reorder(result)
}

# --- Topology manipulation helpers (unchanged) ----------------------------

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
    # Star tree -- nothing to collapse
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
