#' Collect suboptimal trees for landscape analysis
#'
#' `SuboptimalTrees()` performs a parsimony search with [`MaximizeParsimony()`]
#' and returns every retained tree within a specified number of steps of the
#' optimum, each annotated with its parsimony score.
#'
#' The retained trees represent the search engine's internal tree pool, whose
#' size is bounded by `maxPool`. Once the pool is filled, trees are chosen for
#' eviction so as to preserve topological diversity.
#'
#' @inheritParams MaximizeParsimony
#' @param tree Optional starting tree (of class `phylo`) or `multiPhylo`; if
#' `NULL`, the search begins from random addition sequence trees.
#' @param maxSuboptimal Numeric: only trees whose score is within
#' `maxSuboptimal` of the best score found will be retained.
#' Corresponds to `poolSuboptimal` in [`SearchControl()`].
#' @param maxPool Integer: maximum number of trees to retain in the pool.
#' Corresponds to `poolMaxSize` in [`SearchControl()`].
#' @param \dots Further arguments passed to [`MaximizeParsimony()`], including
#' scoring options (`concavity`, `inapplicable`, ...), search effort
#' (`maxReplicates`, `maxSeconds`, `effort`, `nThreads`, `verbosity`),
#' and [`SearchControl()`] fields.
#'
#' @return `SuboptimalTrees()` returns a `multiPhylo` object listing the
#' retained trees, each supplied with a `score` attribute of each tree that
#' records the score it attained. The global `scores` attribute compiles these
#' scores in a numeric vector.
#' [`Suboptimality()`] reports each tree's excess over the optimum.
#'
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' \donttest{
#' set.seed(0)
#' trees <- SuboptimalTrees(dataset, maxSuboptimal = 3, maxReplicates = 5,
#'                          nThreads = 1, maxPool = 100L, verbosity = 0)
#' table(attr(trees, "scores"))  # parsimony lengths of retained trees
#' table(Suboptimality(trees))   # excess over the optimum
#' }
#' @seealso
#' [`MaximizeParsimony()`] performs the underlying search;
#' [`Bremer()`][Bremer] uses this pool for approximate decay indices;
#' [`Suboptimality()`] summarises the excess length of each tree.
#' @template MRS
#' @family split support functions
#' @family tree scoring
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

#' Compare scores of suboptimal trees
#'
#' @param trees List of trees, perhaps generated through `SuboptimalTrees()`,
#' each of which must bear a numeric attribute `score` giving its score.
#' @param normalize Logical stating whether to normalize results to lowest
#' score.
#' @return `Suboptimality()` returns a numeric vector listing, for each tree,
#' its difference in score from the optimal (lowest) within `trees`.
#' @export
Suboptimality <- function (trees, normalize = FALSE) {
  scores <- vapply(trees, attr, double(1), "score")
  
  # Return:
  if (normalize) {
    (scores - min(scores)) / min(scores)
  } else {
    scores - min(scores)
  }
}
