#' Collect suboptimal trees for landscape analysis
#'
#' `SuboptimalTrees()` performs a parsimony search with [`MaximizeParsimony()`]
#' and returns *every* tree retained within a specified number of steps of the
#' optimum, each annotated with its parsimony score.  This exposes the shape of
#' the parsimony landscape near the optimum -- for example to visualise the
#' distribution of near-optimal scores, to measure tree-to-tree distances among
#' competing resolutions, or as input to a fast approximate
#' [`Bremer()`][Bremer] support calculation (`method = "pool"`).
#'
#' The retained trees are the search engine's internal tree pool.  Its size is
#' bounded by `maxPool`, so memory use is capped even when many trees fall
#' within `maxSuboptimal` steps; once the pool is full the engine evicts trees
#' to preserve topological diversity.  The pool therefore reflects the islands
#' the search actually visited, and is not an exhaustive enumeration of every
#' tree within `maxSuboptimal` steps: raise `maxReplicates` / `maxSeconds`
#' (passed via `...`) for a denser sample.
#'
#' @param dataset A phylogenetic data matrix of class `phyDat`, as accepted by
#' [`MaximizeParsimony()`].
#' @param tree Optional starting tree (of class `phylo`) or `multiPhylo`; if
#' `NULL` (default) the search begins from random addition sequence trees.
#' @param maxSuboptimal Numeric: retain trees scoring up to this many steps
#' worse than the best tree found (in the search's optimality units -- integer
#' steps under equal weights, fractional under implied weights or profile
#' parsimony).  Sets `poolSuboptimal` in [`SearchControl()`].
#' @param maxPool Integer: maximum number of trees to retain in the pool.  Sets
#' `poolMaxSize` in [`SearchControl()`]; raised from the search default (100)
#' so that a suboptimal sample is not prematurely truncated.
#' @param \dots Further arguments passed to [`MaximizeParsimony()`], including
#' scoring options (`concavity`, `inapplicable`, ...) and search effort
#' (`maxReplicates`, `maxSeconds`, `strategy`, `nThreads`, `verbosity`).
#' Named [`SearchControl()`] fields may also be passed here to override the
#' constructed control.
#'
#' @return A `multiPhylo` object listing the retained trees, best tree(s)
#' first.  Each tree carries a `score` attribute giving its parsimony length,
#' and the object carries a `scores` attribute: a numeric vector of those
#' lengths aligned with the returned trees.  [`Suboptimality()`] reports each
#' tree's excess over the optimum.
#'
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' \donttest{
#' trees <- SuboptimalTrees(dataset, maxSuboptimal = 3, maxReplicates = 8,
#'                          verbosity = 0)
#' attr(trees, "scores")        # parsimony length of each retained tree
#' Suboptimality(trees)         # excess over the optimum
#' }
#' @seealso
#' [`MaximizeParsimony()`] performs the underlying search;
#' [`Bremer()`][Bremer] uses this pool for approximate decay indices;
#' [`Suboptimality()`] summarises the excess length of each tree.
#' @template MRS
#' @family split support functions
#' @export
SuboptimalTrees <- function(dataset, tree = NULL,
                            maxSuboptimal = 5, maxPool = 1000L, ...) {
  if (!is.numeric(maxSuboptimal) || length(maxSuboptimal) != 1L ||
      maxSuboptimal < 0) {
    stop("`maxSuboptimal` must be a single non-negative number.")
  }
  maxPool <- as.integer(maxPool)
  if (is.na(maxPool) || maxPool < 1L) {
    stop("`maxPool` must be a single positive integer.")
  }

  # SuboptimalTrees() manages `control` (to size the pool) and `collapse` (which
  # must stay FALSE to retain the full suboptimal sample).  A whole `control=`
  # or `collapse=` argument arriving through `...` would otherwise collide with
  # these explicit arguments ("matched by multiple actual arguments").  Strip
  # them; individual SearchControl() fields (e.g. `ratchetCycles = `) still pass
  # through `...` and are merged by MaximizeParsimony().
  dots <- list(...)
  if ("control" %in% names(dots)) {
    warning("`control` is managed by SuboptimalTrees(); pass individual ",
            "SearchControl() fields (e.g. `ratchetCycles = `) via `...` instead.")
    dots[["control"]] <- NULL
  }
  if ("collapse" %in% names(dots)) {
    if (!identical(dots[["collapse"]], FALSE)) {
      warning("`collapse` is forced to FALSE by SuboptimalTrees() so the full ",
              "suboptimal pool is returned.")
    }
    dots[["collapse"]] <- NULL
  }
  # `poolSuboptimal`/`poolMaxSize` are the raw SearchControl() fields that
  # `maxSuboptimal`/`maxPool` set below.  Passed via `...` they would flow into
  # MaximizeParsimony()'s dots-override-control merge and SILENTLY win over the
  # constructed control (last-writer-wins), defeating `maxSuboptimal`/`maxPool`.
  # Strip them with a warning, directing the caller to the blessed arguments.
  managed <- intersect(c("poolSuboptimal", "poolMaxSize"), names(dots))
  if (length(managed)) {
    warning("Set the suboptimal-pool size via `maxSuboptimal` / `maxPool`, not ",
            "`", paste(managed, collapse = "` / `"),
            "`; the latter would silently override the former and ",
            if (length(managed) > 1L) "are" else "is", " ignored.")
    dots[managed] <- NULL
  }

  control <- SearchControl(poolSuboptimal = maxSuboptimal, poolMaxSize = maxPool)

  do.call(MaximizeParsimony,
          c(list(dataset = dataset, tree = tree, control = control,
                 collapse = FALSE), dots))
}
