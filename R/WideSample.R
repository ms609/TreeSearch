#' Select a topologically diverse subset of trees
#'
#' Selects `n` trees from a `multiPhylo` object that are as topologically
#' distinct from one another as possible, by solving the Max-Min Diversity
#' Problem (MMDP): maximize the *minimum* pairwise distance within the chosen
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
#'   \item{FarFirst (`effort = 1`)}{Greedy farthest-first selection
#'     (Gonzalez 1985) from a deterministic peripheral seed.  Fast and
#'     matrix-free: it reads distances through an on-demand column oracle and
#'     never materializes the full distance matrix, so it is the only feasible
#'     option for very large tree sets.  RNG-free, hence deterministic.}
#'   \item{DropAdd (`effort = 2`)}{Drop-add tabu search
#'     (Porumbel et al. 2011): a ~99%-optimal heuristic that terminates at a
#'     deterministic plateau.  RNG-free, hence deterministic.  Requires the full
#'     distance matrix.}
#'   \item{Grasp (`effort = 3`)}{GRASP with path relinking
#'     (Resende et al. 2010): attains the highest \eqn{T_k} of the package's
#'     heuristics, at a cost that grows steeply with `n`.  Requires the full
#'     distance matrix.  Draws on the session RNG, so the particular trees it
#'     returns vary between runs unless you call [set.seed()] first (the
#'     achieved diversity is essentially unaffected).}
#'   \item{exact (`effort = 4`)}{Node-packing integer program
#'     (Sayyady & Fathi 2016): the proven optimum.  The solver is now
#'     sparse-matrix and heuristic warm-started, so it is practical up to a few
#'     hundred trees; it needs the \pkg{highs} package.  The optimal
#'     *diversity* is deterministic, but when several subsets are tied-optimal
#'     the particular trees returned can vary with the session RNG.}
#' }
#'
#' With `effort = NULL` (default) the tier is chosen automatically from
#' `length(trees)`: the exact solver for small sets (up to ~200 trees, when
#' \pkg{highs} is available), DropAdd while the distance matrix is affordable to
#' build, and FarFirst beyond that.  Grasp (`effort = 3`) is never selected
#' automatically -- its cost grows steeply with `n`, so it is opt-in.  A dense
#' distance matrix is roughly `8 * length(trees)^2` bytes (about 1.1 GB at
#' 12,000 trees, 12.8 GB at 40,000), so for the largest sets only the
#' matrix-free FarFirst tier is reachable.
#'
#' Two size thresholds govern automatic selection; tune them for the host
#' machine with [options()] rather than per call:
#' \describe{
#'   \item{`WideSample.buildCeiling`}{Largest `length(trees)` for
#'     which a dense distance matrix is built from a distance function (default
#'     `12000`; ~1.1 GB).  Beyond it only the matrix-free FarFirst tier is
#'     reachable from a function (a pre-computed matrix is always honoured).}
#'   \item{`WideSample.exactCeiling`}{Largest `length(trees)` at
#'     which automatic selection reaches the exact tier (default `200`).}
#' }
#'
#' @param trees A `multiPhylo` object, or a single `phylo` (coerced silently).
#' @param n Integer: number of trees to retain. If `n >= length(trees)`, all
#'   trees are returned unchanged.  If `n == 1`, the single most central tree
#'   (the medoid) is returned.
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
#' @param effort Integer solver tier, or `NULL` (default) to choose
#'   automatically by `length(trees)`.  `1` = FarFirst (fast, matrix-free),
#'   `2` = DropAdd (~99%-optimal, deterministic), `3` = Grasp (highest-quality
#'   heuristic, higher cost), `4` = exact optimum.  Forcing `effort` 2, 3 or 4
#'   with a *distance function* on a tree set too large to build the matrix is
#'   an error rather than a silent downgrade; pass a pre-computed `dist` or use
#'   `effort = 1` for such sets.
#' @param maxSeconds Numeric: wall-clock budget, in seconds, for the
#'   refinement (`effort = 2`, `3`) and exact (`effort = 4`) tiers; ignored by
#'   the matrix-free FarFirst tier, which is effectively instantaneous.  The
#'   heuristic tiers terminate at an internal deterministic plateau, usually
#'   well within this budget, so it acts as a safety cap; a value small enough
#'   to interrupt the plateau makes the result machine-dependent.  Default `60`.
#'
#' @return A `multiPhylo` object of length `min(n, length(trees))`.
#'   Attributes of the input (e.g. `score`, `hits_to_best`) are preserved.
#'
#' @examples
#' library("TreeTools")
#' trees <- as.phylo(0:99, nTip = 8)
#'
#' # Fast FarFirst subsample (deterministic, matrix-free)
#' sub10 <- WideSample(trees, 10, effort = 1)
#' length(sub10)  # 10
#'
#' \donttest{
#' # Automatic tier selection (exact at this size when 'highs' is installed,
#' # otherwise the DropAdd heuristic)
#' auto10 <- WideSample(trees, 10)
#'
#' # Pre-computed distances
#' dists <- TreeDist::ClusteringInfoDistance(trees)
#' sub5 <- WideSample(trees, 5, dist = dists)
#'
#' # Highest-quality heuristic (Grasp); set a seed for a reproducible selection
#' set.seed(1)
#' sub8 <- WideSample(trees, 8, effort = 3)
#'
#' # Force the exact optimum on a small set (needs the 'highs' package)
#' if (requireNamespace("highs", quietly = TRUE)) {
#'   sub4 <- WideSample(trees[1:20], 4, effort = 4)
#' }
#' }
#'
#' @references
#' \insertRef{Gonzalez1985}{TreeSearch}
#'
#' \insertRef{Porumbel2011}{TreeSearch}
#'
#' \insertRef{Resende2010}{TreeSearch}
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
    effort = NULL,
    maxSeconds = 60
) {
  # Build ceiling: largest N for which we materialize a dense N x N matrix from
  # a distance function. ~1.1 GB at 12,000; as.matrix.dist overflows near
  # 46,340 (the dist half-vector exceeds .Machine$integer.max).
  buildCeiling <- getOption("WideSample.buildCeiling", 12000L)
  # Exact ceiling: largest N at which auto-selection reaches the exact tier.
  # MaxMin::ExactMaxMin() is now a sparse-matrix, heuristic-warm-started solver
  # (~20x faster than the dense form), practical to a few hundred trees at the
  # small `n` of interest; beyond that the node-packing IP wall bites (the
  # MaxMin optimum sits near the diameter, where the threshold graph is
  # near-complete). Kept conservative because the IP cost turns on `n` and
  # instance structure, not on `length(trees)` alone.
  exactCeiling <- getOption("WideSample.exactCeiling", 200L)

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
    # Return:
    return(trees)
  }
  if (n == 0L) {
    # Return:
    return(.SubsetMultiPhylo(trees, integer(0)))
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

  # A single tree has no pairwise distance to maximize; return the medoid (the
  # most central tree) as the most representative single choice. Independent of
  # `effort`/`maxSeconds`, so handled before they are validated.
  if (n == 1L) {
    # Return:
    return(.SubsetMultiPhylo(
      trees, .WideSampleMedoid(dist, trees, nTrees, dmat, buildCeiling)
    ))
  }

  if (!is.null(effort)) {
    effort <- as.integer(effort)
    if (length(effort) != 1L || is.na(effort) || !effort %in% 1:4) {
      stop("`effort` must be NULL, 1, 2, 3, or 4")
    }
  }

  if (!is.numeric(maxSeconds) || length(maxSeconds) != 1L ||
      is.na(maxSeconds) || maxSeconds <= 0) {
    stop("`maxSeconds` must be a single positive number (or Inf)")
  }

  # --- select the solver tier on (matrix-available, N) ----------------------
  tier <- .SelectWideSampleTier(effort, matrixAvailable, nTrees,
                                buildCeiling, exactCeiling)

  # --- build the matrix when the chosen tier needs one (tiers 2-4) ----------
  # Tier 1 (FarFirst) stays matrix-free: it reads distances through a column
  # oracle, so it never builds an N x N matrix it was not already handed. The
  # matrix-bound tiers (DropAdd, Grasp, exact) have no oracle path.
  if (tier > 1L && !matrixAvailable) {
    # .SelectWideSampleTier guarantees nTrees <= buildCeiling here.
    dmat <- as.matrix(dist(trees))
    matrixAvailable <- TRUE
  }

  # --- dispatch -------------------------------------------------------------
  # switch() on an integer selects by position, so the cases MUST stay in tier
  # order (1, 2, 3, 4); the backtick labels are cosmetic.
  idx <- switch(
    tier,
    # Tier 1: FarFirst fed by a column oracle. Reading the oracle from a
    # supplied matrix or from the on-demand tree callback feeds FarFirst the
    # identical distances, so the selection does not depend on whether distances
    # were pre-computed; the deterministic peripheral seed keeps it RNG-free.
    `1` = {
      colFn <- if (matrixAvailable) {
        function(i) dmat[, i]
      } else {
        .WideSampleColumnOracle(dist, trees, nTrees)
      }
      MaxMin::FarFirst(n, colFn, N = nTrees)
    },
    # Tier 2: DropAdd returns the bare (sorted) index vector; it runs to its
    # deterministic plateau, with `maxSeconds` as a safety cap.
    `2` = MaxMin::DropAdd(n, dmat, maxSeconds = maxSeconds,
                          progress = FALSE),
    # Tier 3: Grasp likewise returns the bare index vector (RNG-dependent).
    `3` = MaxMin::Grasp(n, dmat, maxSeconds = maxSeconds),
    # Tier 4: exact solver returns a list; take its `$indices`.
    `4` = {
      if (nTrees > exactCeiling) {
        warning("Exact MMDP (effort = 4) on ", nTrees,
                " trees may be very slow; consider effort = 2 (DropAdd) ",
                "or 3 (Grasp), or a larger `maxSeconds`.",
                immediate. = TRUE)
      }
      MaxMin::ExactMaxMin(n, dmat, maxSeconds = maxSeconds,
                          progress = FALSE)$indices
    }
  )

  # FarFirst returns farthest-first (selection) order; sort to ascending tree
  # order so the subset preserves the input ordering. A no-op for tiers 2-4,
  # which already return ascending indices.
  .SubsetMultiPhylo(trees, as.integer(idx))
}

#' Choose the WideSample solver tier
#'
#' Keyed on whether a distance matrix is already available and on
#' `length(trees)`, never on N alone: a supplied matrix keeps the higher tiers
#' reachable past the build ceiling, whereas a distance function past the
#' ceiling cannot reach them (building the matrix would exhaust memory). The
#' exact tier is additionally gated on a (smaller) exact ceiling and on the
#' \pkg{highs} package being installed; Grasp (`effort = 3`) is reachable only
#' when forced, never auto-selected.
#' @return Integer tier (1, 2, 3 or 4); errors when a forced effort is
#'   unreachable.
#' @keywords internal
.SelectWideSampleTier <- function(effort, matrixAvailable, nTrees, ceiling,
                                  exactCeiling = 200L,
                                  highsAvailable =
                                    requireNamespace("highs", quietly = TRUE)) {
  if (is.null(effort)) {
    if (nTrees <= exactCeiling && highsAvailable) {
      # Return:
      return(4L)                       # exact, when affordable and available
    }
    if (matrixAvailable || nTrees <= ceiling) {
      # Return:
      return(2L)                       # DropAdd (build matrix if needed)
    }
    # Return:
    return(1L)                         # FarFirst, matrix-free
  }
  if (effort == 1L) {
    # Return:
    return(1L)
  }
  # effort 2 (DropAdd), 3 (Grasp) and 4 (exact) all need the full matrix.
  if (matrixAvailable || nTrees <= ceiling) {
    # Return:
    return(effort)
  }
  stop("`effort = ", effort, "` needs a distance matrix, but ", nTrees,
       " trees exceeds the build ceiling (", ceiling, ") and no pre-computed ",
       "`dist` was supplied. Use `effort = 1` (FarFirst) for sets ",
       "this large, or pass a pre-computed distance matrix.")
}

#' The medoid tree, for the single-tree (`n == 1`) case
#'
#' Returns the index of the most central tree -- the medoid, minimizing summed
#' distance to all others. Uses the distance matrix when one is available or
#' affordable to build; when only a distance function is supplied for a set too
#' large to build a matrix, the central medoid is not affordable, so the
#' deterministic peripheral seed ([MaxMin::FarFirst()] with `m = 1`) is returned
#' as a matrix-free fallback.
#' @return Integer index (1-based) of the selected tree.
#' @keywords internal
.WideSampleMedoid <- function(dist, trees, nTrees, dmat, buildCeiling) {
  if (is.null(dmat) && nTrees <= buildCeiling) {
    dmat <- as.matrix(dist(trees))
  }
  if (!is.null(dmat)) {
    # Medoid: smallest summed distance to the rest (diagonal is zero, so it
    # does not bias the sum). which.min breaks ties on the smallest index.
    # Return:
    which.min(rowSums(dmat))
  } else {
    colFn <- .WideSampleColumnOracle(dist, trees, nTrees)
    # Return:
    as.integer(MaxMin::FarFirst(colFn, m = 1L, N = nTrees, progress = FALSE))
  }
}

#' Build a column-oracle closure for the matrix-free FarFirst path
#'
#' Returns a function of one 1-based index `i` giving the distances from tree
#' `i` to every tree, as required by the distance-column oracle path of
#' [MaxMin::FarFirst()]. Probes
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
  # Return:
  result
}
