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
#' @param treeSample Controls how many input trees are used during Phase 1
#'   (driven search).  CID verification cost scales linearly with the number
#'   of input trees, so with large tree sets (hundreds or thousands of trees),
#'   using a representative subsample for the search phase and verifying
#'   against the full set afterwards can be faster without sacrificing quality.
#'
#'   - `"auto"` (default): automatically selects a subsample size based on
#'     the number of tips.  Benchmarking on a 1449-taxon mammal bootstrap
#'     dataset (Lemoine _et al._ 2018) showed that consensus quality
#'     (measured as full-set MCI) saturates at 50--100 input trees and can
#'     _degrade_ with larger subsamples under fixed time budgets, because
#'     the additional CID cost displaces search replicates.  The auto
#'     heuristic is: `min(T, max(50, 2 * n_tip))`, capped at `T`.
#'   - An integer: use exactly this many trees (sampled without replacement).
#'   - `Inf` or `NULL`: use all trees (no subsampling).
#'
#'   Subsampling applies only to Phase 1.  Phases 2 (collapse/resolve) and 3
#'   (rogue taxon dropping) always score against the **full** input tree set,
#'   and the returned `score` attribute reflects the full-set MCI.
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
#' @param screeningTopK Integer: number of top MRP-screened candidates to
#'   evaluate via full MCI scoring per TBR clip.  Default `1` uses the
#'   single best MRP candidate (original behaviour).  Values > 1 (e.g., `5`)
#'   score the top-k MRP candidates and accept the one with the best MCI,
#'   catching moves where MRP and MCI rankings disagree.  Cost scales
#'   linearly with `screeningTopK`.
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
                          maxReplicates = 5L,
                          targetHits = 2L,
                          maxSeconds = 0,
                          nThreads = getOption("mc.cores", 1L),
                          treeSample = "auto",
                          collapse = TRUE,
                          neverDrop = FALSE,
                          maxDrop = min(5L, ceiling(NTip(trees[[1]]) / 10)),
                          control = SearchControl(),
                          screeningK = 7,
                          screeningTolerance = 0,
                          screeningTopK = 1L,
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
  
  nTree <- length(trees)
  
  # --- Tree subsampling for Phase 1 ---
  # CID verification cost is O(T * n_tip) per candidate. With large tree
  # sets, this dominates wall time and displaces search replicates.
  # Benchmarking (mammals 1449-taxon bootstrap, 1000 trees, 50-/100-tip
  # subsets, 60-/120-s budgets) showed:
  #   - MCI quality saturates at T_sub ~ 50-100, regardless of tip count
  #   - At 100 tips, quality *degrades* beyond T_sub ~ 100 under fixed
  #     time budgets (CID cost displaces exploration)
  #   - MRP split deduplication makes the Fitch screening layer insensitive
  #     to T (97% dedup at T=1000), so the bottleneck is CID verification
  # Heuristic: use min(T, max(50, 2 * n_tip)) trees for Phase 1.
  # Phases 2-3 always use the full set.
  searchTrees <- .SubsampleTrees(trees, nTree, nTip, treeSample, verbosity)
  
  # Phase 1: C++ driven search (on subsample if applicable)
  result <- .CIDDrivenSearch(searchTrees, tipLabels, nTip,
                              maxReplicates, targetHits, maxSeconds,
                              nThreads, control,
                              screeningK, screeningTolerance, screeningTopK,
                              verbosity)
  
  # Phases 2-3 use the full tree set for accurate scoring.
  # Re-score Phase 1 result against full set if we subsampled.
  subsampled <- length(searchTrees) < nTree
  if (subsampled) {
    cidData <- .MakeCIDData(trees, tipLabels)
    subScore <- attr(result, "score")
    fullScore <- .ScoreTree(result, cidData)
    attr(result, "score") <- fullScore
    if (verbosity > 0L) {
      message("  Subsample MCI: ", signif(-subScore, 6),
              " -> full-set MCI: ", signif(-fullScore, 6),
              " (", nTree, " trees)")
    }
  }
  
  # Phase 2: Collapse/resolve refinement (full tree set)
  if (collapse) {
    if (!exists("cidData", inherits = FALSE)) {
      cidData <- .MakeCIDData(trees, tipLabels)
    }
    result <- .CollapseRefine(result, cidData, verbosity)
  }
  
  # Phase 3: rogue taxon dropping (full tree set)
  if (!isTRUE(neverDrop)) {
    if (!exists("cidData", inherits = FALSE)) {
      cidData <- .MakeCIDData(trees, tipLabels)
    }
    result <- .RogueRefine(result, cidData, neverDrop, maxDrop,
                           "ratchet",       # method for re-optimization
                           2L, 2L,          # light ratchet for re-opt
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


# Select a subsample of input trees for Phase 1 search.
#
# With large tree sets, CID verification O(T * n_tip) per candidate
# dominates wall time. Subsampling to T_sub trees for the search phase
# preserves quality (the unique split landscape saturates early) while
# freeing time for more search replicates.
#
# Returns: multiPhylo of at most tSub trees (random without replacement).
.SubsampleTrees <- function(trees, nTree, nTip, treeSample, verbosity) {
  # Resolve treeSample to an integer target
  if (is.null(treeSample) || (is.numeric(treeSample) && is.infinite(treeSample))) {
    return(trees) # No subsampling
  }
  
  if (identical(treeSample, "auto")) {
    # Heuristic: min(T, max(50, 2 * nTip)).
    # At 50 tips, T_sub = 100; at 100 tips, T_sub = 200; at 25 tips, T_sub = 50.
    # Benchmark showed quality saturates by T_sub ~ 50-100 regardless of
    # tip count; the 2*nTip term provides a safety margin for larger trees
    # where the split landscape is richer.
    tSub <- min(nTree, max(50L, 2L * nTip))
  } else if (is.numeric(treeSample) && length(treeSample) == 1L) {
    tSub <- min(nTree, max(2L, as.integer(treeSample)))
  } else {
    stop("`treeSample` must be \"auto\", a positive integer, Inf, or NULL.",
         call. = FALSE)
  }
  
  if (tSub >= nTree) return(trees)
  
  if (verbosity > 0L) {
    message("  Subsampling ", tSub, " of ", nTree,
            " input trees for Phase 1 search")
  }
  
  idx <- sample.int(nTree, tSub)
  trees[idx]
}


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
                              screeningTopK,
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
  
  # CID scoring is ~50-100x more expensive per evaluation than parsimony,

  # so use lighter search parameters than the parsimony defaults.
  # Each replicate takes ~10-25s at 50 tips; 10 replicates is usually
  # sufficient for CID landscape convergence.
  result <- ts_cid_consensus(
    splitMatrices = splitMats,
    nTip = nTip,
    normalize = FALSE,
    maxReplicates = maxReplicates,
    targetHits = targetHits,
    tbrMaxHits = .NullOr(ctrl[["tbrMaxHits"]], 1L),
    ratchetCycles = .NullOr(ctrl[["ratchetCycles"]], 2L),
    ratchetPerturbProb = .NullOr(ctrl[["ratchetPerturbProb"]], 0.25),
    ratchetPerturbMode = .NullOr(ctrl[["ratchetPerturbMode"]], 0L),
    ratchetAdaptive = .NullOr(ctrl[["ratchetAdaptive"]], FALSE),
    driftCycles = .NullOr(ctrl[["driftCycles"]], 1L),
    driftAfdLimit = .NullOr(ctrl[["driftAfdLimit"]], 3L),
    driftRfdLimit = .NullOr(ctrl[["driftRfdLimit"]], 0.1),
    xssRounds = .NullOr(ctrl[["xssRounds"]], if (useSectors) 1L else 0L),
    xssPartitions = .NullOr(ctrl[["xssPartitions"]], 4L),
    rssRounds = .NullOr(ctrl[["rssRounds"]], if (useSectors) 1L else 0L),
    cssRounds = .NullOr(ctrl[["cssRounds"]], 0L),
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
    screeningTolerance = screeningTolerance,
    screeningTopK = screeningTopK,
    scoreTol = .NullOr(ctrl[["scoreTol"]], 1e-5),
    plateauReps = .NullOr(ctrl[["plateauReps"]], 3L)
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
# Returns negated mean MCI (lower = better), consistent with C++ cid_score().
# Delegates to ts_cid_score_trees for the optimised C++ path (precomputed
# hash index, log2 values, bounded early exit, persistent scratch buffers).
.CIDScorer <- function(parent, child, dataset) {
  ts_cid_score_trees(dataset$inputSplitsRaw, dataset$nTip,
                     list(cbind(parent, child)))[1L]
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
#
# All collapse/resolve candidates per pass are batch-scored in a single
# ts_cid_score_trees() call, amortising CidData construction (hash-index
# build, log2 precomputation) across the candidate set.
.CollapseRefine <- function(tree, cidData, verbosity = 0L) {
  tipLabels <- cidData$tipLabels
  nTip      <- cidData$nTip
  splitMats <- cidData$inputSplitsRaw
  edge   <- tree[["edge"]]
  parent <- edge[, 1L]
  child  <- edge[, 2L]

  bestScore <- ts_cid_score_trees(splitMats, nTip, list(edge))[1L]

  if (verbosity > 0L) {
    message("  - Collapse/resolve refinement. Starting MCI: ",
            signif(-bestScore, 6))
  }

  improved <- TRUE
  while (improved) {
    improved <- FALSE

    # --- Collapse pass: batch-score all collapse candidates ---
    internalEdges <- which(child > nTip)
    if (length(internalEdges) > 0L) {
      candidates <- lapply(internalEdges, function(idx) {
        cand <- .CollapseSpecificEdge(parent, child, idx, nTip)
        cbind(cand[[1L]], cand[[2L]])
      })
      scores  <- ts_cid_score_trees(splitMats, nTip, candidates)
      bestIdx <- which.min(scores)
      if (scores[[bestIdx]] < bestScore - sqrt(.Machine[["double.eps"]])) {
        bestEdge  <- candidates[[bestIdx]]
        parent    <- bestEdge[, 1L]
        child     <- bestEdge[, 2L]
        bestScore <- scores[[bestIdx]]
        improved  <- TRUE
        if (verbosity > 1L) {
          message("    * Collapsed edge -> MCI ", signif(-bestScore, 6))
        }
        next
      }
    }

    # --- Resolve pass: batch-score all resolve candidates ---
    degrees   <- tabulate(parent, nbins = max(parent))
    polyNodes <- which(degrees > 2L)
    if (length(polyNodes) > 0L) {
      candidates <- list()
      for (node in polyNodes) {
        childEdges <- which(parent == node)
        nChildren  <- length(childEdges)
        if (nChildren <= 2L) next
        for (a in seq_len(nChildren - 1L)) {
          for (b in seq(a + 1L, nChildren)) {
            cand <- .ResolveSpecificPair(parent, child, node,
                                         childEdges[c(a, b)], nTip)
            candidates <- c(candidates, list(cbind(cand[[1L]], cand[[2L]])))
          }
        }
      }
      if (length(candidates) > 0L) {
        scores  <- ts_cid_score_trees(splitMats, nTip, candidates)
        bestIdx <- which.min(scores)
        if (scores[[bestIdx]] < bestScore - sqrt(.Machine[["double.eps"]])) {
          bestEdge  <- candidates[[bestIdx]]
          parent    <- bestEdge[, 1L]
          child     <- bestEdge[, 2L]
          bestScore <- scores[[bestIdx]]
          improved  <- TRUE
          if (verbosity > 1L) {
            message("    * Resolved polytomy -> MCI ", signif(-bestScore, 6))
          }
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
# Uses C++ prescreen (ts_cid_prescreen_rogue) for the first iteration
# (no tips yet dropped); subsequent iterations use R-level pruning.
# Accepts drops based on prescreen score alone (no per-candidate
# re-optimization), then does one light re-optimization at the end.
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
    # Accept the single best drop directly (no per-candidate re-optimization)
    bestIdx <- which.min(prescreenScores)
    if (prescreenScores[bestIdx] < bestScore - sqrt(.Machine[["double.eps"]])) {
      tip <- names(prescreenScores)[bestIdx]
      reducedTips <- setdiff(currentTips, tip)
      allDropped <- c(droppedTips, tip)
      reducedInputTrees <- .PruneTrees(originalTrees, allDropped)
      reducedCidData <- .MakeCIDData(reducedInputTrees, reducedTips)
      tree <- DropTip(tree, tip)
      # Score the pruned tree on the properly pruned input trees
      bestScore <- .ScoreTree(tree, reducedCidData)
      currentTips <- reducedTips
      cidData <- reducedCidData
      droppedTips <- allDropped
      improved <- TRUE
      if (verbosity > 0L) {
        message("    * Dropped '", tip, "' -> MCI ",
                signif(-bestScore, 6),
                " (", length(currentTips), " tips)")
      }
    }
  }

  # Light re-optimization after all drops are decided
  if (length(droppedTips) > 0L && length(currentTips) >= 8L) {
    if (verbosity > 0L) {
      message("  - Re-optimizing after rogue removal")
    }
    reoptResult <- .TopologySearch(
      tree, cidData, "ratchet",
      max(1L, ratchIter %/% 2L), ratchHits,
      searchIter, searchHits, collapse,
      max(0L, verbosity - 1L))
    reoptScore <- attr(reoptResult, "score")
    if (!is.null(reoptScore) && reoptScore < bestScore) {
      tree <- reoptResult
      bestScore <- reoptScore
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


# Pre-screen rogue candidates via C++ bit-masking.
# For each tip in `droppable`, masks out the tip from pre-built split data
# and scores in C++ -- no R-level DropTip/as.Splits needed.
# When tips have already been dropped, falls back to the R-level path.
.PrescreenMarginalNID <- function(tree, cidData, droppable,
                                  originalTrees, allTipLabels,
                                  alreadyDropped) {
  if (length(alreadyDropped) == 0L) {
    # Fast C++ path: mask each tip from the full split data
    tipLabels <- cidData$tipLabels
    dropIdx <- match(droppable, tipLabels)
    scores <- ts_cid_prescreen_rogue(
      cidData$inputSplitsRaw, cidData$nTip,
      tree[["edge"]], dropIdx
    )
    names(scores) <- droppable
    return(scores)
  }

  # Fallback: tips already dropped, need pruned input trees
  vapply(droppable, function(tip) {
    reducedTips <- setdiff(tree[["tip.label"]], tip)
    allDropped  <- c(alreadyDropped, tip)
    reducedInputTrees <- .PruneTrees(originalTrees, allDropped)
    splitMats <- lapply(reducedInputTrees,
                        function(tr) unclass(as.Splits(tr, reducedTips)))
    reducedTree <- DropTip(tree, tip)
    ts_cid_score_trees(splitMats, length(reducedTips),
                       list(reducedTree[["edge"]]))[1L]
  }, double(1))
}


.PruneTrees <- function(trees, tipsToDrop) {
  if (length(tipsToDrop) == 0L) return(trees)
  pruned <- lapply(trees, DropTip, tip = tipsToDrop)
  class(pruned) <- "multiPhylo"
  pruned
}


.ScoreTree <- function(tree, cidData) {
  ts_cid_score_trees(cidData$inputSplitsRaw, cidData$nTip,
                     list(tree[["edge"]]))[1L]
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


# Try inserting tipLabel at every edge of tree; return the best position.
# Batch-scores all insertion candidates in a single ts_cid_score_trees() call
# so CidData is built once (not once per insertion position).
.BestInsertion <- function(tree, tipLabel, originalTrees,
                           otherDropped) {
  currentTips <- c(tree[["tip.label"]], tipLabel)
  inputTrees  <- if (length(otherDropped) > 0L) {
    .PruneTrees(originalTrees, otherDropped)
  } else {
    originalTrees
  }
  testCidData <- .MakeCIDData(inputTrees, currentTips)
  splitMats   <- testCidData$inputSplitsRaw
  nTip        <- length(currentTips)

  nEdge      <- nrow(tree[["edge"]])
  candidates <- lapply(seq_len(nEdge),
                       function(i) .InsertTipAtEdge(tree, tipLabel, i)[["edge"]])

  scores  <- ts_cid_score_trees(splitMats, nTip, candidates)
  bestIdx <- which.min(scores)

  bestEdge <- candidates[[bestIdx]]
  bestTree <- .EdgeListToPhylo(bestEdge[, 1L], bestEdge[, 2L], currentTips)
  list(tree = bestTree, score = scores[[bestIdx]], cidData = testCidData)
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
