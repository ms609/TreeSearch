#' Bremer (decay) support
#'
#' `Bremer()` calculates the Bremer support (decay index)
#' \insertCite{Bremer1988,Bremer1994}{TreeSearch} of each clade in a reference
#' tree: the number of extra steps required before the clade is no longer
#' present in an optimal tree.  Formally, for a clade _C_,
#' \deqn{Bremer(C) = L(\neg C) - L^*}{Bremer(C) = L(not C) - L*}
#' where \eqn{L(\neg C)}{L(not C)} is the length of the shortest tree that does
#' *not* contain _C_ and \eqn{L^*}{L*} is the length of the most-parsimonious
#' tree.  Larger values indicate better-supported clades.
#'
#' Two engines are available:
#'
#' \describe{
#'   \item{`method = "constraint"` (default, rigorous)}{Runs one *converse-
#'     constraint* search per clade, forcing that clade to be absent and taking
#'     the length of the shortest resulting tree.  This directly targets
#'     \eqn{L(\neg C)}{L(not C)} and needs only bounded memory (one tree per
#'     clade).  It is the reliable choice for reported support values.}
#'   \item{`method = "pool"` (fast, approximate)}{Collects a pool of suboptimal
#'     trees with [`SuboptimalTrees()`] and, for each clade, takes the shortest
#'     retained tree that lacks it.  This is quick but **over-estimates**
#'     support (the minimum over a sampled subset can only exceed the true
#'     minimum) and is **right-censored** at the sampling depth `maxBremer`: a
#'     clade broken by no retained tree is reported as `Inf` and flagged in the
#'     `censored` attribute.  Use it as a preview, not for publication values.}
#' }
#'
#' The reference `tree` supplies the clades to be evaluated.  Pass a single
#' `phylo` (e.g. one most-parsimonious tree) or a `multiPhylo` search result, in
#' which case the strict consensus is used and support is calculated only for
#' its resolved bipartitions.  Scoring options (`concavity`, `inapplicable`,
#' ...) **must match those used to find the trees**, or the extra-step counts
#' will be meaningless; they default to equal-weights Fitch parsimony.
#'
#' @param tree A most-parsimonious tree (`phylo`) whose clades are to be
#' evaluated, or a `multiPhylo` whose strict consensus provides them.  A
#' `multiPhylo` returned by [`MaximizeParsimony()`] additionally supplies the
#' optimal score \eqn{L^*}{L*} via its `score` attribute.
#' @param dataset A phylogenetic data matrix of class `phyDat`.
#' @param method Character: `"constraint"` (default) for rigorous converse-
#' constraint searches, or `"pool"` for the fast suboptimal-pool approximation.
#' @param maxBremer Numeric: the largest decay value to resolve.  Under
#' `method = "pool"` this is the pool's suboptimality depth (defaults to `10`);
#' clades not broken within it are censored.  Under `method = "constraint"` it
#' is ignored (exact values are computed).
#' @param optimalScore Numeric: the optimal tree length \eqn{L^*}{L*}.  If
#' `NULL` (default) it is taken from `attr(tree, "score")` when available, or
#' computed by search.  Supplying it (or a scored `multiPhylo`) avoids a
#' redundant search.
#' @param format Character specifying return format, as in [`JackLabels()`]:
#' `"numeric"` (default) returns named numeric values for further analysis;
#' `"character"` returns a vector shaped for `phylo$node.label`.
#' @inheritParams MaximizeParsimony
#' @param \dots Further arguments passed to [`MaximizeParsimony()`] /
#' [`SuboptimalTrees()`], e.g. `maxReplicates`, `maxSeconds`, `strategy`,
#' `nThreads`, `verbosity`.
#'
#' @return A numeric vector (or, if `format = "character"`, a
#' `phylo$node.label`-shaped character vector) giving the Bremer support of each
#' resolved clade in `tree`, named by node number (the row names of
#' [`TreeTools::as.Splits()`]).  Annotate a plot with
#' [`TreeTools::LabelSplits()`], or assign to `tree$node.label`.  Under
#' `method = "pool"` the numeric result carries a logical `censored` attribute
#' marking clades whose support exceeds `maxBremer`.
#'
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' \donttest{
#' trees <- MaximizeParsimony(dataset, maxReplicates = 8, verbosity = 0)
#' decay <- Bremer(trees, dataset, maxReplicates = 8, verbosity = 0)
#' decay
#'
#' # Annotate a tree
#' TreeTools::LabelSplits(trees[[1]], decay)
#' }
#' @references \insertAllCited{}
#' @template MRS
#' @seealso
#' Other clade support measures: [`JackLabels()`], [`SiteConcordance`];
#' [`SuboptimalTrees()`] collects the pool used by `method = "pool"`.
#' @family split support functions
#' @importFrom TreeTools as.Splits NTip SplitFrequency TipLabels
#' @export
Bremer <- function(tree, dataset,
                   concavity = Inf, extended_iw = TRUE, xpiwe_r = 0.5,
                   xpiwe_max_f = 5, hierarchy = NULL, inapplicable = "bgs",
                   hsj_alpha = 1.0,
                   method = c("constraint", "pool"),
                   maxBremer = Inf, optimalScore = NULL,
                   format = "numeric", ...) {
  method <- match.arg(method)

  scoringArgs <- list(concavity = concavity, extended_iw = extended_iw,
                      xpiwe_r = xpiwe_r, xpiwe_max_f = xpiwe_max_f,
                      hierarchy = hierarchy, inapplicable = inapplicable,
                      hsj_alpha = hsj_alpha)

  # --- Resolve the reference tree, its clades, and L* ---
  ref <- .BremerReference(tree, dataset, optimalScore)

  if (length(ref$splitNames) == 0L) {
    warning("Reference tree has no resolved internal clades; ",
            "nothing to calculate.")
    return(.BremerFormat(setNames(numeric(0), character(0)),
                         logical(0), ref$reference, format))
  }

  engine <- if (method == "pool") .BremerPool else .BremerConstraint
  res <- engine(ref, dataset, scoringArgs, maxBremer, list(...))

  .BremerFormat(res$bremer, res$censored, ref$reference, format)
}

# Resolve reference tree + clades (Splits) + optimal score L*.
#' @importFrom TreeTools as.Splits NTip
.BremerReference <- function(tree, dataset, optimalScore) {
  Lstar <- optimalScore
  if (inherits(tree, "multiPhylo")) {
    if (is.null(Lstar)) {
      s <- attr(tree, "score")
      if (!is.null(s) && is.finite(s)) Lstar <- s
    }
    reference <- if (length(tree) == 1L) tree[[1L]] else ape::consensus(tree, p = 1)
  } else if (inherits(tree, "phylo")) {
    reference <- tree
  } else {
    stop("`tree` must be a `phylo` or `multiPhylo` object.")
  }

  splits <- as.Splits(reference, tipLabels = names(dataset))
  splitNames <- rownames(as.matrix(splits))
  if (is.null(splitNames)) splitNames <- character(0)

  list(reference = reference, splits = splits, splitNames = splitNames,
       Lstar = Lstar)
}

# Format a node-indexed numeric result, mirroring JackLabels().
#' @importFrom TreeTools NTip
.BremerFormat <- function(values, censored, reference, format) {
  numericFmt <- c("numeric", "number", "double")
  characterFmt <- c("character", "text")
  returnMode <- c(rep("numeric", length(numericFmt)),
                  rep("character", length(characterFmt)))[
                    pmatch(tolower(format), c(numericFmt, characterFmt))]
  if (is.na(returnMode)) returnMode <- "numeric"

  switch(returnMode,
    "character" = {
      ret <- character(reference[["Nnode"]])
      if (length(values)) {
        idx <- as.integer(names(values)) - NTip(reference)
        ret[idx] <- as.character(values)
      }
      ret
    },
    {
      if (length(censored)) attr(values, "censored") <- censored
      values
    })
}

# Approximate Bremer from a pool of suboptimal trees.
# Returns list(bremer = named numeric, censored = logical).
#' @importFrom TreeTools SplitFrequency
.BremerPool <- function(ref, dataset, scoringArgs, maxBremer, dots) {
  poolK <- if (is.finite(maxBremer)) maxBremer else 10
  if (!is.finite(maxBremer)) {
    message("method = \"pool\": using maxBremer = ", poolK,
            " for the pool depth (set `maxBremer` to control it).")
  }

  pool <- do.call(SuboptimalTrees,
                  c(list(dataset = dataset, maxSuboptimal = poolK),
                    scoringArgs, dots))
  poolScores <- attr(pool, "scores")
  if (is.null(poolScores)) {
    poolScores <- do.call(TreeLength,
                          c(list(tree = pool, dataset = dataset), scoringArgs))
  }
  poolBest <- min(poolScores)

  Lstar <- ref$Lstar
  tol <- 1e-8
  if (is.null(Lstar)) {
    Lstar <- poolBest
  } else if (poolBest < Lstar - tol) {
    warning("Pool search found a tree (length ", signif(poolBest, 7),
            ") shorter than the supplied optimalScore (", signif(Lstar, 7),
            "); adopting the shorter length as L*.")
    Lstar <- poolBest
  }

  # Per-clade displaying matrix: refFreq gives the node names; then test each
  # pool tree individually (SplitFrequency of a length-1 forest is 0/1).
  refFreq <- SplitFrequency(ref$reference, pool)
  nSplits <- length(refFreq)
  displayMat <- vapply(seq_along(pool), function(j) {
    SplitFrequency(ref$reference, pool[j])
  }, double(nSplits))
  dim(displayMat) <- c(nSplits, length(pool))

  bremer <- numeric(nSplits)
  censored <- logical(nSplits)
  for (i in seq_len(nSplits)) {
    nonDisplaying <- displayMat[i, ] < 0.5
    if (any(nonDisplaying)) {
      bremer[i] <- min(poolScores[nonDisplaying]) - Lstar
    } else {
      bremer[i] <- Inf
      censored[i] <- TRUE
    }
  }
  names(bremer) <- names(refFreq)
  names(censored) <- names(refFreq)
  list(bremer = bremer, censored = censored)
}

# Rigorous Bremer by converse-constraint search: for each clade, the shortest
# tree forced to LACK it.  Uses the negative-constraint engine wired into
# MaximizeParsimony() via the internal `.negativeConstraint` argument.
#
# The converse search keeps the full search machinery -- Wagner starts, TBR,
# ratchet, sectorial search and NNI perturbation -- which all re-optimize
# through the negative-constraint-guarded TBR and so stay in the space of trees
# lacking the clade; the tree pool additionally rejects any tree that displays
# it.  Only the phases that accept score-worsening moves through an unguarded
# path (drift, in-sector drift, simulated annealing) are disabled -- they would
# otherwise wander onto the clade and, being worse-accepting, report it as
# unsupported.  Runs serially: the pool guard is on the serial search path.
.BremerConstraint <- function(ref, dataset, scoringArgs, maxBremer, dots) {
  disabled <- list(
    driftCycles = 0L, sectorGoDrift = 0L, sectorDriftCycles = 0L,
    annealCycles = 0L
  )
  overrides <- dots[intersect(names(dots), names(disabled))]
  converseFixed <- modifyList(disabled, overrides)
  # nThreads is forced to 1 (the negative-constraint pool guard is on the
  # serial search path only); everything else the user passes flows through.
  passthrough <- dots[setdiff(names(dots), c(names(disabled), "nThreads"))]

  # L* (the unconstrained optimum) is found by a full-strength search -- the
  # worse-accepting phases are safe here because there is no forbidden clade.
  Lstar <- ref$Lstar
  if (is.null(Lstar)) {
    res0 <- do.call(MaximizeParsimony,
                    c(list(dataset = dataset, collapse = TRUE, nThreads = 1L),
                      scoringArgs, passthrough))
    Lstar <- attr(res0, "score")
  }

  runConverse <- function(negSplit) {
    args <- c(list(dataset = dataset, collapse = TRUE, nThreads = 1L,
                   .negativeConstraint = negSplit),
              scoringArgs, converseFixed, passthrough)
    attr(do.call(MaximizeParsimony, args), "score")
  }

  splitNames <- ref$splitNames
  scores <- vapply(seq_along(splitNames), function(i) {
    s <- runConverse(ref$splits[[i]])
    # A negative score is the engine's "no tree found" sentinel (empty pool):
    # the converse search failed to locate any tree lacking the clade.
    if (!is.finite(s) || s < 0) NA_real_ else s
  }, double(1))

  if (anyNA(scores)) {
    warning(sum(is.na(scores)), " converse-constraint search(es) found no tree ",
            "lacking the clade; reported as NA. Increase `maxReplicates`.")
  }

  trueBest <- suppressWarnings(min(c(Lstar, scores), na.rm = TRUE))
  if (is.finite(trueBest) && trueBest < Lstar - 1e-8) {
    warning("A converse-constraint search found a tree (length ",
            signif(trueBest, 7), ") shorter than the optimal score (",
            signif(Lstar, 7), "); the original search was suboptimal. ",
            "Adopting the shorter length as L*.")
    Lstar <- trueBest
  }

  bremer <- setNames(scores - Lstar, splitNames)
  list(bremer = bremer, censored = logical(length(bremer)))
}
