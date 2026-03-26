#' Tree search using successive approximations
#'
#' Searches for a tree that is optimal under the Successive Approximations 
#' criterion \insertCite{Farris1969}{TreeSearch}.
#'
#' @inheritParams TreeTools::Renumber
#' @inheritParams MaximizeParsimony
#' @param outgroup if not NULL, taxa on which the tree should be rooted
#' @param k Constant for successive approximations, see Farris 1969 p. 379
#' @param maxSuccIter maximum iterations of successive approximation
#' @param ratchetHits Number of replicates.
#'   Internally capped at 100 and passed to the C++ engine as `maxReplicates`.
#' @param searchHits Convergence criterion: stop after finding this many
#'   trees with the best score.
#'   Internally capped at 10 and passed to the C++ engine as `targetHits`.
#' @param searchIter Unused (retained for backward compatibility).
#' @param ratchetIter Controls ratchet intensity within each replicate.
#'   Converted to `ratchetCycles` (approximately `ratchetIter / 500`,
#'   capped at 10).
#' @param suboptimal Retain trees that are this proportion less optimal
#'   than the optimal tree.
#' 
#' @return `SuccessiveApproximations()` returns a list of class `multiPhylo`
#' containing optimal (and slightly suboptimal, if suboptimal > 0) trees.
#'
#' @references 
#' \insertAllCited{}
#'
#' @importFrom ape consensus root
#' @family custom search functions
#' @export
SuccessiveApproximations <- function (tree, dataset, outgroup = NULL, k = 3,
                                      maxSuccIter = 20, ratchetHits = 100,
                                      searchHits = 50, searchIter = 500,
                                      ratchetIter = 5000, verbosity = 0,
                                      suboptimal = 0.1,
                                      concavity = Inf,
                                      constraint = NULL,
                                      extended_iw = TRUE,
                                      xpiwe_r = 0.5,
                                      xpiwe_max_f = 5) {

  if (k < 1) stop("k should be at least 1, see Farris 1969 p.379")

  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be of class `phyDat`.")
  }

  # Profile parsimony: prepare data
  useProfile <- identical(concavity, "profile")
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
    concavity <- Inf
  }
  if (is.finite(concavity) && concavity <= 0) {
    stop("`concavity` must be positive (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }

  nTip <- length(dataset)

  # Extract data for C++ engine
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = nTip, byrow = TRUE)
  weight <- at$weight
  levels <- at$levels

  # Prepare constraint
  consArgs <- .PrepareConstraint(constraint = constraint, dataset = dataset)

  # Profile parsimony: extract info_amounts
  profileArgs <- list()
  if (useProfile) {
    infoAmounts <- attr(dataset, "info.amounts")
    if (!is.null(infoAmounts) && length(infoAmounts) > 0L) {
      profileArgs$infoAmounts <- infoAmounts
    }
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
    saK = as.double(k),
    maxSAIter = as.integer(maxSuccIter),
    maxReplicates = as.integer(min(ratchetHits, 100L)),
    targetHits = as.integer(min(searchHits, 10L)),
    tbrMaxHits = 1L,
    ratchetCycles = as.integer(min(ceiling(ratchetIter / 500), 10L)),
    min_steps = if (is.finite(concavity))
      as.integer(MinimumLength(dataset, compress = TRUE)) else integer(0),
    concavity = as.double(concavity),
    xpiwe = useXpiwe,
    xpiwe_r = as.double(xpiwe_r),
    xpiwe_max_f = as.double(xpiwe_max_f),
    obs_count = if (useXpiwe) obsCount else integer(0)
  )
  result <- do.call(ts_successive_approx, c(searchArgs, consArgs, profileArgs))

  if (result$converged && verbosity > 0) {
    message("Successive approximations converged after ",
            result$sa_iterations, " iteration(s).")
  } else if (!result$converged) {
    message("Stability not reached after ", result$sa_iterations,
            " iteration(s).")
  }

  # Reconstruct phylo from C++ edge matrix
  if (nrow(result$edge) == 0L) {
    tr <- if (!missing(tree) && inherits(tree, "phylo")) tree
          else AdditionTree(dataset)
    attr(tr, "score") <- result$score
  } else {
    tr <- structure(
      list(edge = result$edge,
           tip.label = names(dataset),
           Nnode = nTip - 1L),
      class = "phylo"
    )
    attr(tr, "score") <- result$score
  }

  if (!is.null(outgroup)) {
    tr <- root(tr, outgroup, resolve.root = TRUE)
  }

  structure(
    list(tr),
    score = result$score,
    sa_iterations = result$sa_iterations,
    converged = result$converged,
    class = "multiPhylo"
  )
}

#' Tree suboptimality
#'
#' How suboptimal is a tree?
#'
#' @param trees list of trees, to include an optimal tree
#' @param proportional logical stating whether to normalise results to lowest
#' score
#' @return `Suboptimality()` returns a vector listing, for each tree, how much
#' its score differs from the optimal (lowest) score.
#' @keywords internal
#' @export
Suboptimality <- function (trees, proportional = FALSE) {
  scores <- vapply(trees, attr, double(1), "score")
  
  # Return:
  if (proportional) {
    (scores - min(scores)) / min(scores)
  } else {
    scores - min(scores)
  }
}

#' @rdname SuccessiveApproximations
#' @return `SuccessiveWeights()` returns the score of a tree, given the
#' weighting instructions specified in the attributes of the dataset.
#' 
#' @keywords internal
#' @export
SuccessiveWeights <- function(tree, dataset) {
  # Data
  if (inherits(dataset, "phyDat")) dataset <- PrepareDataSA(dataset)
  if (!inherits(dataset, "saDat")) {
    stop("Invalid data type; prepare data with PhyDat() or PrepareDataSA().")
  }
  at <- attributes(dataset)
  weight <- at[["weight"]]
  sa.weights <- at[["sa.weights"]]
  if (is.null(sa.weights)) sa.weights <- rep.int(1, length(weight))
  steps <- CharacterLength(tree, dataset, compress = TRUE)
  
  # Return:
  sum(steps * sa.weights * weight)
}

PrepareDataSA <- function (dataset) {
# Written with reference to phangorn:::prepareDataFitch
  at <- attributes(dataset)
  nam <- at[["names"]]
  nLevel <- length(at[["levels"]])
  nChar <- at[["nr"]]
  cont <- attr(dataset, "contrast")
  nTip <- length(dataset)
  
  at[["names"]] <- NULL
  powers.of.2 <- 2L ^ c(0L:(nLevel - 1L))
  tmp <- cont %*% powers.of.2
  tmp <- as.integer(tmp)
  dataset <- unlist(dataset, FALSE, FALSE)
  ret <- tmp[dataset] 
  ret <- as.integer(ret)
  attributes(ret) <- at
  inappLevel <- which(at[["levels"]] == "-")
  attr(ret, "dim") <- c(nChar, nTip)  
  applicableTokens <- setdiff(powers.of.2, 2 ^ (inappLevel - 1))
  attr(ret, "split.sizes") <- t(apply(ret, 1, function(x) vapply(applicableTokens, function (y) sum(x == y), integer(1))))
  dimnames(ret) <- list(NULL, nam)
  attr(ret, "bootstrap") <- list("split.sizes", "sa.weights")
  class(ret) <- "saDat"
  ret
}
