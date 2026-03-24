#' Native C++ scoring for custom search functions
#'
#' These functions provide a native C++ scoring interface for use with
#' [`TreeSearch()`], [`Ratchet()`], and [`Jackknife()`], replacing the
#' MorphyLib-based defaults (`PhyDat2Morphy`, `UnloadMorphy`,
#' `MorphyLength`, `MorphyBootstrap`).
#'
#' @name NativeSearch
#' @family custom search functions
NULL

#' @describeIn NativeSearch Prepare a phyDat dataset for native C++ scoring.
#' Replaces [`PhyDat2Morphy()`] as the `InitializeData` function.
#'
#' @param dataset A phyDat object.
#' @param concavity Concavity constant for implied weighting, `Inf` for
#'   equal weights, or `"profile"` for profile parsimony.
#' @return `PrepareNativeData()` returns a list containing pre-extracted data
#'   matrices suitable for repeated scoring with [`NativeLength()`].
#' @export
PrepareNativeData <- function(dataset, concavity = Inf) {
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  useProfile <- identical(concavity, "profile")
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
    infoAmounts <- attr(dataset, "info.amounts")
  } else {
    infoAmounts <- NULL
  }

  iw <- !useProfile && is.finite(concavity)
  if (iw) {
    if (concavity <= 0) {
      stop("`concavity` must be positive (or Inf for equal weights, ",
           "or \"profile\" for profile parsimony).")
    }
    if (!("min.length" %in% names(attributes(dataset)))) {
      dataset <- PrepareDataIW(dataset)
    }
    minSteps <- as.integer(attr(dataset, "min.length"))
  } else {
    minSteps <- integer(0)
  }

  at <- attributes(dataset)
  weight <- at[["weight"]]

  list(
    contrast = at[["contrast"]],
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = as.integer(weight),
    levels = at[["levels"]],
    min_steps = minSteps,
    concavity = if (iw) as.double(concavity) else -1.0,
    info_amounts = infoAmounts,
    original_weight = as.integer(weight),
    index = at[["index"]]
  )
}

#' @describeIn NativeSearch No-op cleanup (native data has no external
#'   resources). Replaces [`UnloadMorphy()`] as the `CleanUpData` function.
#' @param dataset Prepared dataset returned by `PrepareNativeData()`.
#' @return `CleanNativeData()` returns `dataset` invisibly.
#' @export
CleanNativeData <- function(dataset) {
  invisible(dataset)
}

#' @describeIn NativeSearch Score a tree using the native C++ Fitch engine.
#' Replaces [`MorphyLength()`] as the `TreeScorer` function.
#'
#' @param parent Integer vector of parent node indices for each edge.
#' @param child Integer vector of child node indices for each edge.
#' @param dataset Prepared dataset from [`PrepareNativeData()`].
#' @param \dots Ignored (present for interface compatibility).
#' @return `NativeLength()` returns a double: the parsimony score of the tree.
#' @export
NativeLength <- function(parent, child, dataset, ...) {
  ts_fitch_score(
    cbind(parent, child),
    dataset[["contrast"]],
    dataset[["tip_data"]],
    dataset[["weight"]],
    dataset[["levels"]],
    dataset[["min_steps"]],
    dataset[["concavity"]],
    dataset[["info_amounts"]]
  )
}

#' @describeIn NativeSearch Bootstrap resampling for use with [`Ratchet()`].
#' Replaces [`MorphyBootstrap()`] as the `Bootstrapper` function.
#'
#' @param edgeList A list containing parent and child integer vectors
#'   (and optionally a score and hit count).
#' @param dataset Prepared dataset from [`PrepareNativeData()`].
#' @param EdgeSwapper A function to rearrange edges, such as
#'   [`NNISwap()`] or [`TBRSwap()`].
#' @param maxIter Maximum rearrangement iterations.
#' @param maxHits Maximum hits of best score before stopping.
#' @param verbosity Verbosity level.
#' @param stopAtPeak,stopAtPlateau Passed to [`EdgeListSearch()`].
#' @param \dots Further parameters passed to the `TreeScorer`.
#' @return `NativeBootstrap()` returns an edge list (list of parent and
#'   child vectors).
#' @export
NativeBootstrap <- function(edgeList, dataset, EdgeSwapper = NNISwap,
                            maxIter, maxHits, verbosity = 1L,
                            stopAtPeak = FALSE, stopAtPlateau = 0L, ...) {
  startWeights <- dataset[["original_weight"]]
  eachChar <- seq_along(startWeights)
  deindexedChars <- rep.int(eachChar, startWeights)
  resampling <- tabulate(sample(deindexedChars, replace = TRUE),
                         length(startWeights))

  # R's copy-on-modify: local copy, caller's dataset unchanged
  dataset[["weight"]] <- as.integer(resampling)

  res <- EdgeListSearch(edgeList[1:2], dataset,
                        TreeScorer = NativeLength,
                        EdgeSwapper = EdgeSwapper,
                        maxIter = maxIter, maxHits = maxHits,
                        stopAtPeak = stopAtPeak,
                        stopAtPlateau = stopAtPlateau,
                        verbosity = verbosity - 1L, ...)

  # Return:
  res[1:2]
}
