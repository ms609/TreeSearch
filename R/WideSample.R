#' Select a topologically diverse subset of trees
#'
#' Selects `n` trees from a `multiPhylo` object that are as topologically
#' distinct from one another as possible, by solving the Max-Min Diversity
#' Problem (MMDP): maximise the *minimum* pairwise distance within the chosen
#' subset.  This is useful when a search returns many most-parsimonious trees
#' and downstream analyses (consensus, tree-space visualization) need a
#' manageable but diverse subset.
#'
#' Uniform random subsampling of MPTs is misleading: the number of trees in a
#' region of tree space reflects the density of the parsimony landscape, not
#' the likelihood or support for that topology.  A random draw over-represents
#' topologies that sit on broad plateaux and under-represents isolated optima.
#' `WideSample()` instead selects for topological *spread*, density-blind, by
#' dispatching to the appropriate MMDP solver from the \pkg{MaxMin}
#' package:
#'
#' \describe{
#'   \item{anchored Gonzalez (`quality = 1`)}{Greedy furthest-point selection
#'     from a deterministic peripheral seed.  Fast; the only feasible option
#'     for very large tree sets, where it runs from an on-demand distance
#'     oracle and never materialises the full distance matrix.}
#'   \item{DropAdd-TS (`quality = 2`)}{Drop-add tabu search
#'     (Porumbel et al. 2011): a ~99%-optimal heuristic.  Requires the full
#'     distance matrix.}
#'   \item{exact (`quality = 3`)}{Node-packing integer program
#'     (Sayyady & Fathi 2016): the proven optimum.  Practical only for small
#'     `length(trees)`; needs the \pkg{highs} package.}
#' }
#'
#' With `quality = NULL` (default) the tier is chosen automatically from
#' `length(trees)`: the exact solver for very small sets (when \pkg{highs} is
#' available), DropAdd-TS while the distance matrix is affordable to build, and
#' anchored Gonzalez beyond that.  A dense distance matrix is roughly
#' `8 * length(trees)^2` bytes (about 1.1 GB at 12,000 trees, 12.8 GB at
#' 40,000), so for the largest sets only the matrix-free anchored Gonzalez is
#' reachable.
#'
#' @param trees A `multiPhylo` object, or a single `phylo` (coerced silently).
#' @param n Integer: number of trees to retain. If `n >= length(trees)`, all
#'   trees are returned unchanged.
#' @param dist Either:
#'   \itemize{
#'     \item A function giving pairwise distances (default:
#'       [TreeDist::ClusteringInfoDistance()]).  It must support the form
#'       `dist(trees)` returning a `dist` object; for the largest tree sets it
#'       is additionally called as `dist(trees[[i]], trees)` and must then
#'       return a numeric vector of length `length(trees)` (the distances from
#'       tree `i` to every tree).  `ClusteringInfoDistance()` satisfies both.
#'     \item A pre-computed `dist` object or square numeric matrix whose size
#'       matches `length(trees)`.
#'   }
#' @param quality Integer solver tier, or `NULL` (default) to choose
#'   automatically by `length(trees)`.  `1` = anchored Gonzalez (fast),
#'   `2` = DropAdd-TS (~99%-optimal), `3` = exact optimum.  Forcing `quality`
#'   2 or 3 with a *distance function* on a tree set too large to build the
#'   matrix is an error rather than a silent downgrade; pass a pre-computed
#'   `dist` or use `quality = 1` for such sets.
#' @param time_budget_s Numeric: wall-clock budget, in seconds, for the
#'   heuristic (`quality = 2`) and exact (`quality = 3`) tiers.  Ignored by the
#'   anchored-Gonzalez tier, which is effectively instantaneous.  Default `2`.
#'
#' @return A `multiPhylo` object of length `min(n, length(trees))`.
#'   Attributes of the input (e.g. `score`, `hits_to_best`) are preserved.
#'
#' @examples
#' library("TreeTools")
#' trees <- as.phylo(0:99, nTip = 8)
#'
#' # Fast anchored-Gonzalez subsample (deterministic, matrix-free)
#' sub10 <- WideSample(trees, 10, quality = 1)
#' length(sub10)  # 10
#'
#' \donttest{
#' # Automatic tier selection (DropAdd-TS heuristic at this size)
#' auto10 <- WideSample(trees, 10)
#'
#' # Pre-computed distances
#' dists <- TreeDist::ClusteringInfoDistance(trees)
#' sub5 <- WideSample(trees, 5, dist = dists)
#'
#' # Force the exact optimum on a small set (needs the 'highs' package)
#' if (requireNamespace("highs", quietly = TRUE)) {
#'   sub4 <- WideSample(trees[1:20], 4, quality = 3)
#' }
#' }
#'
#' @references
#' \insertRef{Porumbel2011}{TreeSearch}
#'
#' \insertRef{Sayyady2016}{TreeSearch}
#'
#' @template MRS
#' @family tree scoring
#' @importFrom TreeDist ClusteringInfoDistance
#' @export
WideSample <- function(
    trees,
    n,
    dist = TreeDist::ClusteringInfoDistance,
    quality = NULL,
    time_budget_s = 2
) {
  # Build ceiling: largest N for which we materialise a dense N x N matrix from
  # a distance function. ~1.1 GB at 12,000; as.matrix.dist overflows near
  # 46,340 (the dist half-vector exceeds .Machine$integer.max).
  buildCeiling <- 12000L

  if (inherits(trees, "phylo")) {
    trees <- c(trees)
  } else if (!inherits(trees, "multiPhylo")) {
    stop("`trees` must be a multiPhylo object")
  }
  nTrees <- length(trees)

  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 0L) {
    stop("`n` must be a single non-negative integer")
  }
  if (n >= nTrees) {
    return(trees)
  }
  if (n == 0L) {
    return(.SubsetMultiPhylo(trees, integer(0)))
  }
  # Tiers 2 and 3 require m >= 2; a single tree has no pairwise distance to
  # maximise, so return one deterministically before tier selection.
  if (n == 1L) {
    return(.SubsetMultiPhylo(trees, 1L))
  }

  if (!is.null(quality)) {
    quality <- as.integer(quality)
    if (length(quality) != 1L || is.na(quality) || !quality %in% 1:3) {
      stop("`quality` must be NULL, 1, 2, or 3")
    }
  }

  # --- classify `dist`: pre-computed matrix vs distance function -------------
  distIsMatrix <- inherits(dist, "dist") ||
    (is.matrix(dist) && is.numeric(dist))
  distIsFun <- is.function(dist)
  if (!distIsMatrix && !distIsFun) {
    stop("`dist` must be a function, a `dist` object, or a numeric matrix")
  }

  dmat <- NULL
  if (distIsMatrix) {
    dmat <- as.matrix(dist)
    if (nrow(dmat) != ncol(dmat)) {
      stop("`dist` matrix must be square")
    }
    if (nrow(dmat) != nTrees) {
      stop("`dist` has ", nrow(dmat), " rows but `trees` has ", nTrees,
           " trees")
    }
  }
  matrixAvailable <- !is.null(dmat)

  # --- select the solver tier on (matrix-available, N) ----------------------
  tier <- .SelectWideSampleTier(quality, matrixAvailable, nTrees, buildCeiling)

  # --- build the matrix when the chosen tier needs one (tiers 2, 3) ---------
  # Tier 1 stays matrix-free: it reads distances through a column oracle, so it
  # never builds an N x N matrix it was not already handed.
  if (tier %in% c(2L, 3L) && !matrixAvailable) {
    # .SelectWideSampleTier guarantees nTrees <= buildCeiling here.
    dmat <- as.matrix(dist(trees))
    matrixAvailable <- TRUE
  }

  # --- dispatch -------------------------------------------------------------
  idx <- switch(
    tier,
    # Tier 1: single anchored Gonzalez, fed by a column oracle. Reading the
    # oracle from a supplied matrix or from the on-demand tree callback feeds
    # Gonzalez() the identical distances, so the selection does not depend
    # on whether distances were pre-computed.
    `1` = {
      colFn <- if (matrixAvailable) {
        function(i) dmat[, i]
      } else {
        .WideSampleColumnOracle(dist, trees, nTrees)
      }
      MaxMin::Gonzalez(colFn, n, N = nTrees)
    },
    `2` = MaxMin::DropAddTS(dmat, n, time_budget_s = time_budget_s)$indices,
    `3` = {
      if (nTrees > 40L) {
        warning("Exact MMDP (quality = 3) on ", nTrees,
                " trees may be very slow; consider quality = 2 (DropAdd-TS) ",
                "or a larger `time_budget_s`.", immediate. = TRUE)
      }
      MaxMin::ExactMaxMin(dmat, n, time_budget_s = time_budget_s)$indices
    }
  )

  .SubsetMultiPhylo(trees, sort(as.integer(idx)))
}

#' Choose the WideSample solver tier
#'
#' Keyed on whether a distance matrix is already available and on
#' `length(trees)`, never on N alone: a supplied matrix keeps the higher tiers
#' reachable past the build ceiling, whereas a distance function past the
#' ceiling cannot reach them (building the matrix would exhaust memory).
#' @return Integer tier (1, 2 or 3); errors when a forced quality is
#'   unreachable.
#' @keywords internal
.SelectWideSampleTier <- function(quality, matrixAvailable, nTrees, ceiling) {
  if (is.null(quality)) {
    if (nTrees <= 40L && requireNamespace("highs", quietly = TRUE)) {
      return(3L)                       # exact, when affordable and available
    }
    if (matrixAvailable || nTrees <= ceiling) {
      return(2L)                       # DropAdd-TS (build matrix if needed)
    }
    return(1L)                         # anchored Gonzalez, matrix-free
  }
  if (quality == 1L) {
    return(1L)
  }
  # quality 2 or 3 both need the full matrix.
  if (matrixAvailable || nTrees <= ceiling) {
    return(quality)
  }
  stop("`quality = ", quality, "` needs a distance matrix, but ", nTrees,
       " trees exceeds the build ceiling (", ceiling, ") and no pre-computed ",
       "`dist` was supplied. Use `quality = 1` (anchored Gonzalez) for sets ",
       "this large, or pass a pre-computed distance matrix.")
}

#' Build a column-oracle closure for the matrix-free Gonzalez path
#'
#' Returns a function of one 1-based index `i` giving the distances from tree
#' `i` to every tree, as required by the distance-column oracle path of
#' [MaxMin::Gonzalez()]. Probes
#' the `(tree, trees)` calling form once up front and fails clearly if the
#' supplied `dist` function does not support it.
#' @keywords internal
.WideSampleColumnOracle <- function(dist, trees, nTrees) {
  probe <- tryCatch(
    dist(trees[[1L]], trees),
    error = function(e) {
      stop("`dist` must accept `dist(trees[[i]], trees)` for tree sets too ",
           "large to build a full distance matrix; calling it raised: ",
           conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.numeric(probe) || length(probe) != nTrees) {
    stop("`dist(trees[[i]], trees)` must return a numeric vector of length ",
         nTrees, "; got ",
         if (is.numeric(probe)) paste0("length ", length(probe))
         else class(probe)[[1L]], ".")
  }
  function(i) as.numeric(dist(trees[[i]], trees))
}

#' Subset a multiPhylo preserving attributes
#' @keywords internal
.SubsetMultiPhylo <- function(trees, idx) {
  saved <- attributes(trees)
  result <- trees[idx]
  # Restore non-standard attributes (e.g. score, hits_to_best)
  standard <- c("names", "class")
  for (nm in setdiff(names(saved), standard)) {
    attr(result, nm) <- saved[[nm]]
  }
  result
}
