# Hierarchy-aware resampling: generates hierarchical weights per replicate
# and calls ts_driven_search with HSJ/xform scoring.
# This is an internal helper called from Resample() when inapplicable != "bgs".
.ResampleHierarchy <- function(dataset, hierarchy, inapplicable, hsj_alpha,
                               method_idx, proportion, nReplicates,
                               contrast, tip_data, weight, levels, nTip,
                               concavity, ratchIter, tbrIter,
                               consArgs, profileArgs, tree) {
  bootstrap <- (method_idx == 2L)

  # Prepare full HSJ args (before resampling)
  hsjBase <- list()
  if (identical(inapplicable, "hsj")) {
    # Get flat blocks grouped by top-level block
    .FlattenOneTop <- function(node) {
      block <- list(
        primary = node$controlling - 1L,
        secondaries = node$dependents - 1L
      )
      child_blocks <- lapply(node$children, .FlattenOneTop)
      c(list(block), unlist(child_blocks, recursive = FALSE))
    }
    hsjBase$blocks_per_top <- lapply(hierarchy, .FlattenOneTop)
    hsjBase$hsjTipLabels <- .BuildTipLabels(dataset)
    hsjBase$hsjAlpha <- as.double(hsj_alpha)
    hsjBase$hsjAbsentState <- .HSJAbsentState(dataset)
  }

  # Prepare full xform args (before resampling)
  xformBase <- list()
  if (identical(inapplicable, "xform")) {
    recoded <- RecodeHierarchy(dataset, hierarchy)
    xformBase$all_chars <- recoded$sankoff_chars
  }

  # Driven search params for resampling context (light search per replicate)
  resampleControl <- SearchControl(
    tbrMaxHits = as.integer(max(tbrIter, 1L)),
    ratchetCycles = as.integer(max(ratchIter, 3L)),
    driftCycles = 0L,
    xssRounds = 0L,
    rssRounds = 0L,
    cssRounds = 0L,
    fuseInterval = 0L,
    poolMaxSize = 1L,
    poolSuboptimal = 0.0
  )
  resampleRuntime <- list(
    maxReplicates = as.integer(max(ratchIter, 5L)),
    targetHits = 2L,
    maxSeconds = 0.0,
    verbosity = 0L,
    nThreads = 1L,
    startEdge = NULL,
    progressCallback = NULL
  )
  resampleScoring <- list(
    min_steps = integer(0),
    concavity = as.double(concavity),
    xpiwe = FALSE,
    xpiwe_r = 0.5,
    xpiwe_max_f = 5.0,
    obs_count = integer(0),
    infoAmounts = profileArgs$infoAmounts
  )

  trees <- vector("list", nReplicates)
  for (r in seq_len(nReplicates)) {
    resamp <- .HierarchicalResampleWeights(
      dataset, hierarchy, bootstrap, proportion
    )

    # Build per-replicate hierarchy args based on retained blocks
    repHsj <- list()
    repXform <- list()

    if (identical(inapplicable, "hsj")) {
      # Expand retained flat blocks (supports bootstrap: block sampled >1 time)
      rep_blocks <- list()
      for (bi in seq_along(resamp$blockCounts)) {
        if (resamp$blockCounts[bi] > 0L) {
          top_blocks <- hsjBase$blocks_per_top[[bi]]
          for (k in seq_len(resamp$blockCounts[bi])) {
            rep_blocks <- c(rep_blocks, top_blocks)
          }
        }
      }
      repHsj$hierarchyBlocks <- rep_blocks
      repHsj$hsjTipLabels <- hsjBase$hsjTipLabels
      repHsj$hsjAlpha <- hsjBase$hsjAlpha
      repHsj$hsjAbsentState <- hsjBase$hsjAbsentState
    }

    if (identical(inapplicable, "xform")) {
      rep_xf <- list()
      for (bi in seq_along(resamp$blockCounts)) {
        if (resamp$blockCounts[bi] > 0L) {
          for (k in seq_len(resamp$blockCounts[bi])) {
            rep_xf <- c(rep_xf, list(xformBase$all_chars[[bi]]))
          }
        }
      }
      repXform$xformChars <- rep_xf
    }

    # Call ts_driven_search with resampled weights
    constraintCfg <- if (length(consArgs) > 0L) consArgs
    hsjCfg <- if (length(repHsj) > 0L) repHsj
    xformCfg <- if (length(repXform) > 0L) repXform

    result <- ts_driven_search(
      contrast, tip_data,
      as.integer(resamp$nonHierarchyWeights), levels,
      resampleControl, resampleRuntime, resampleScoring,
      constraintCfg, hsjCfg, xformCfg
    )

    # Extract best tree
    if (result$pool_size > 0L && length(result$trees) > 0L) {
      tr <- structure(
        list(edge = result$trees[[1L]],
             tip.label = names(dataset),
             Nnode = nTip - 1L),
        class = "phylo"
      )
      attr(tr, "score") <- result$best_score
    } else {
      tr <- if (!is.null(tree) && inherits(tree, "phylo")) tree
            else AdditionTree(dataset)
      attr(tr, "score") <- result$best_score
    }
    trees[[r]] <- tr
  }

  structure(trees, class = "multiPhylo")
}


#' Resampling under custom search criteria
#'
#' @inheritParams MaximizeParsimony
#' @param tree (optional) A bifurcating tree of class \code{\link[ape]{phylo}}.
#' Consulted only as a fallback: if a replicate's search reports no
#' best-score tree, `tree` (if supplied) is returned in its place.  Unlike
#' [`MaximizeParsimony()`]'s `tree` argument, it is not used as a warm-start
#' topology.
#' @param method Unambiguous abbreviation of `jackknife` or `bootstrap`
#' specifying how to resample characters.  Note that jackknife is considered
#' to give more meaningful results.
#'
#' @param proportion Numeric between 0 and 1 specifying what proportion of
#' characters to retain under jackknife resampling.
#'
#' @param ratchIter Numeric: governs search effort per replicate, mapped to
#' `max(ratchIter, 5)` independent search replicates and `max(ratchIter, 3)`
#' parsimony-ratchet cycles per replicate.
#' @param tbrIter Numeric: maximum number of times the best score may be hit
#' during a \acronym{TBR} rearrangement pass before it stops (mapped to the
#' underlying search engine's `tbrMaxHits` control).
#' @param finalIter,maxHits,tolerance,verbosity Deprecated and without
#' effect.  These were parameters of the pre-2.0.0, \pkg{MorphyLib}-based
#' implementation of `Resample()`; the native search engine that replaced it
#' has no equivalent controls (in particular, no progress-reporting hook), so
#' supplying a non-`NULL` value now only issues a deprecation warning.
#' @param nThreads Integer: number of parallel threads for search replicates,
#' as for [`MaximizeParsimony()`].  Only takes effect for `inapplicable =
#' "bgs"` (the default) with `nReplicates > 1`; resampling under `"hsj"` or
#' `"xform"`, or a single replicate, always runs on a single thread.
#' @section Resampling:
#' Note that bootstrap support is a measure of the amount of data supporting
#' a split, rather than the amount of confidence that should be afforded the
#' grouping.
#' "Bootstrap support of 100% is not enough, the tree must also be correct" 
#' \insertCite{Phillips2004}{TreeSearch}.
#' See discussion in \insertCite{Egan2006;textual}{TreeSearch};
#' \insertCite{Wagele2009;textual}{TreeSearch};
#' \insertCite{Simmons2011}{TreeSearch};
#' \insertCite{Kumar2012;textual}{TreeSearch}.
#' 
#' For a discussion of suitable search parameters in resampling estimates, see
#' \insertCite{Muller2005;textual}{TreeSearch}.
#' The user should decide whether to start each resampling
#' from the optimal tree (which may be quicker, but result in overestimated 
#' support values as searches get stuck in local optima close to the 
#' optimal tree) or a random tree (which may take longer as more rearrangements
#' are necessary to find an optimal tree on each iteration).
#' 
#' For other ways to estimate clade concordance, see [`SiteConcordance()`].
#' 
#' @param nReplicates Integer specifying how many resample replicates to run.
#' Default `1L` runs a single replicate (original behaviour).
#' When `> 1`, all replicates are run in a single call, optionally in parallel.
#' @param ... Unused; retained for backward compatibility.
#'
#' @return `Resample()` returns a `multiPhylo` object containing one best tree
#' per resample replicate.
#' @family split support functions
#' @encoding UTF-8
#' @export
Resample <- function(dataset, tree, method = "jack", proportion = 2 / 3,
                     ratchIter = 1L, tbrIter = 8L, finalIter = NULL,
                     maxHits = NULL, concavity = Inf,
                     tolerance = NULL,
                     constraint, verbosity = NULL,
                     nReplicates = 1L, nThreads = 1L,
                     hierarchy = NULL, inapplicable = "bgs",
                     hsj_alpha = 1.0,
                     extended_iw = TRUE,
                     xpiwe_r = 0.5,
                     xpiwe_max_f = 5,
                     ...) {

  if (!is.null(finalIter) || !is.null(maxHits) || !is.null(tolerance) ||
      !is.null(verbosity)) {
    .Deprecated(msg = paste(
      "In `Resample()`: `finalIter`, `maxHits`, `tolerance` and `verbosity`",
      "have had no effect since `Resample()` was rewritten to use the",
      "native search engine, and are deprecated.  They will be removed in",
      "a future release."
    ))
  }

  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be of class `phyDat`.")
  }

  method_idx <- pmatch(tolower(method), c("jackknife", "bootstrap"))
  if (is.na(method_idx)) {
    stop("`method` must be either \"jackknife\" or \"bootstrap\".")
  }

  nReplicates <- as.integer(max(nReplicates, 1L))
  nThreads <- as.integer(nThreads)

  # Validate proportion for jackknife
  index <- attr(dataset, "index")
  if (method_idx == 1L) {
    nKept <- ceiling(proportion * length(index))
    if (nKept < 1L) {
      stop("No characters retained. `proportion` must be positive.")
    }
    if (nKept == length(index)) {
      stop("`proportion` too high; no characters deleted.")
    }
  }

  # --- Validate inapplicable-handling parameters ---
  inapplicable <- tolower(inapplicable)
  if (inapplicable == "brazeau") inapplicable <- "bgs"
  inapplicable <- match.arg(inapplicable, c("bgs", "hsj", "xform"))
  if (inapplicable != "bgs") {
    if (is.null(hierarchy)) {
      stop("A `hierarchy` is required when inapplicable = \"", inapplicable,
           "\". See ?CharacterHierarchy.")
    }
    if (!inherits(hierarchy, "CharacterHierarchy")) {
      stop("`hierarchy` must be a CharacterHierarchy object.")
    }
    ValidateHierarchy(hierarchy, dataset)
  }
  if (!is.numeric(hsj_alpha) || length(hsj_alpha) != 1L ||
      hsj_alpha < 0 || hsj_alpha > 1) {
    stop("`hsj_alpha` must be a single number in [0, 1].")
  }

  # Profile parsimony: prepare data
  useProfile <- identical(concavity, "profile")
  if (useProfile) {
    if (inapplicable != "bgs") {
      stop("Profile parsimony is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
    dataset <- PrepareDataProfile(dataset)
    concavity <- Inf
  }
  if (is.finite(concavity) && inapplicable != "bgs") {
    stop("Implied weighting is not currently supported with inapplicable = \"",
         inapplicable, "\".")
  }
  if (is.finite(concavity) && concavity <= 0) {
    stop("`concavity` must be positive (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }

  # C++ engine path
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- .ScaleWeight(at$weight)
  levels <- at$levels
  nTip <- length(dataset)

  # Prepare constraint
  consArgs <- .PrepareConstraint(
    constraint = if (!missing(constraint)) constraint,
    dataset = dataset
  )

  # Profile parsimony: extract info_amounts
  profileArgs <- list()
  if (useProfile) {
    infoAmounts <- attr(dataset, "info.amounts")
    if (!is.null(infoAmounts) && length(infoAmounts) > 0L) {
      profileArgs$infoAmounts <- infoAmounts
    }
  }

  # --- Hierarchy-aware resampling path ---
  # When inapplicable != "bgs", resample at the unit level (free chars +
  # hierarchy blocks) and run driven_search per replicate with HSJ/xform
  # scoring.
  if (inapplicable != "bgs" && !is.null(hierarchy)) {
    return(.ResampleHierarchy(
      dataset = dataset, hierarchy = hierarchy, inapplicable = inapplicable,
      hsj_alpha = hsj_alpha, method_idx = method_idx, proportion = proportion,
      nReplicates = nReplicates,
      contrast = contrast, tip_data = tip_data, weight = weight,
      levels = levels, nTip = nTip, concavity = concavity,
      ratchIter = ratchIter, tbrIter = tbrIter,
      consArgs = consArgs, profileArgs = profileArgs,
      tree = if (!missing(tree)) tree else NULL
    ))
  }

  # XPIWE: compute per-pattern observed-taxa counts
  useXpiwe <- isTRUE(extended_iw) && is.finite(concavity) && !useProfile
  if (useXpiwe) {
    obsCount <- .ObsCount(dataset)
  }

  searchArgs <- list(
    contrast = contrast,
    tip_data = tip_data,
    weight = weight,
    levels = levels,
    bootstrap = (method_idx == 2L),
    jackProportion = proportion,
    maxReplicates = as.integer(max(ratchIter, 5L)),
    targetHits = 2L,
    tbrMaxHits = as.integer(max(tbrIter, 1L)),
    ratchetCycles = as.integer(max(ratchIter, 3L)),
    min_steps = if (is.finite(concavity))
      as.integer(MinimumLength(dataset, compress = TRUE)) else integer(0),
    concavity = as.double(concavity),
    xpiwe = useXpiwe,
    xpiwe_r = as.double(xpiwe_r),
    xpiwe_max_f = as.double(xpiwe_max_f),
    obs_count = if (useXpiwe) obsCount else integer(0)
  )

  if (nReplicates > 1L) {
    # Batch mode: run all replicates at once (optionally in parallel)
    batchArgs <- c(searchArgs,
                   list(nReplicates = nReplicates, nThreads = nThreads),
                   consArgs, profileArgs)
    result <- do.call(ts_parallel_resample, batchArgs)

    trees <- vector("list", nReplicates)
    for (r in seq_len(nReplicates)) {
      em <- result$edges[[r]]
      if (nrow(em) == 0L) {
        tr <- if (!missing(tree) && inherits(tree, "phylo")) tree
              else AdditionTree(dataset)
      } else {
        tr <- structure(
          list(edge = em,
               tip.label = names(dataset),
               Nnode = nTip - 1L),
          class = "phylo"
        )
      }
      attr(tr, "score") <- result$scores[r]
      trees[[r]] <- tr
    }
    return(structure(trees, class = "multiPhylo"))
  }

  # Single-replicate path (original behavior)
  result <- do.call(ts_resample_search, c(searchArgs, consArgs, profileArgs))

  if (nrow(result$edge) == 0L) {
    tr <- if (!missing(tree) && inherits(tree, "phylo")) tree
          else AdditionTree(dataset)
    attr(tr, "score") <- result$score
    return(structure(list(tr), class = "multiPhylo"))
  }

  tr <- structure(
    list(edge = result$edge,
         tip.label = names(dataset),
         Nnode = nTip - 1L),
    class = "phylo"
  )
  attr(tr, "score") <- result$score

  structure(list(tr), class = "multiPhylo")
}
