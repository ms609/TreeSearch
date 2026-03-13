#' Fast driven parsimony search using the C++ tree search engine
#'
#' Performs a multi-replicate driven search for most-parsimonious trees,
#' combining random addition sequence (Wagner) starting trees, TBR
#' rearrangement, exclusive sectorial search (XSS), ratchet perturbation,
#' and tree fusing -- all in compiled C++.
#'
#' Each replicate builds a random Wagner tree, optimizes it with TBR,
#' applies sectorial search and ratchet escape, then adds the result to a
#' pool of unique topologies.
#' Periodically, tree fusing recombines the best trees in the pool.
#' The search stops when the best score has been independently discovered
#' `targetHits` times, or `maxReplicates` replicates have been completed.
#'
#' @inheritParams MaximizeParsimony
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates (default: 100).
#' @param targetHits Integer: stop when the best score has been found
#'   independently this many times (default: `max(10, NTip / 5)`).
#' @param tbrMaxHits Integer: within each TBR phase, the number of
#'   equal-score trees to explore on the plateau before stopping
#'   (default: `max(20, NTip)`).
#' @param ratchetCycles Integer: number of ratchet perturbation cycles
#'   per replicate.
#' @param ratchetPerturbProb Numeric: probability of zeroing each
#'   character's weight during ratchet perturbation.
#' @param xssRounds Integer: number of exclusive sectorial search rounds.
#' @param xssPartitions Integer: number of partitions for XSS.
#' @param sectorMinSize,sectorMaxSize Integer: minimum and maximum clade
#'   sizes for sectorial search.
#' @param fuseInterval Integer: fuse pool trees every this many replicates.
#' @param fuseAcceptEqual Logical: accept equally-scoring fused trees?
#' @param poolMaxSize Integer: maximum number of trees retained in pool.
#' @param poolSuboptimal Numeric: retain trees within this many steps of
#'   the best score (0 = optimal only).
#'
#' @return A `multiPhylo` object containing the best tree(s) found, with
#'   attributes:
#'   \describe{
#'     \item{`score`}{Best parsimony score.}
#'     \item{`replicates`}{Number of replicates completed.}
#'     \item{`hits_to_best`}{Number of independent discoveries of the best
#'       score.}
#'   }
#'
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' result <- MaximizeParsimony2(dataset, maxReplicates = 3L, targetHits = 2L)
#' result
#' attr(result, "score")
#'
#' @family tree scoring
#' @seealso [MaximizeParsimony()] for the original R-level search.
#' @importFrom TreeTools NTip RenumberTips RootTree
#' @importFrom cli cli_alert_success cli_alert_info cli_alert_warning
#' @export
MaximizeParsimony2 <- function(
    dataset,
    tree,
    maxReplicates = 100L,
    targetHits = max(10L, as.integer(NTip(dataset) / 5)),
    tbrMaxHits = max(20L, NTip(dataset)),
    ratchetCycles = 10L,
    ratchetPerturbProb = 0.04,
    xssRounds = 3L,
    xssPartitions = 4L,
    sectorMinSize = 6L,
    sectorMaxSize = 50L,
    fuseInterval = 3L,
    fuseAcceptEqual = FALSE,
    poolMaxSize = 100L,
    poolSuboptimal = 0.0,
    verbosity = 1L
) {

  # --- Input validation ---
  if (!inherits(dataset, "phyDat")) {
    stop("`dataset` must be a phyDat object.")
  }

  nTip <- length(dataset)
  if (nTip < 4L) {
    stop("Need at least 4 taxa for tree search.")
  }

  # --- Starting tree ---
  if (missing(tree)) {
    tree <- AdditionTree(dataset)
  } else if (inherits(tree, "multiPhylo")) {
    tree <- tree[[1L]]
  }
  if (!inherits(tree, "phylo")) {
    stop("`tree` must be of class 'phylo'.")
  }

  # Make bifurcating if needed
  if (dim(tree[["edge"]])[1] != 2L * tree[["Nnode"]]) {
    tree <- MakeTreeBinary(tree)
    if (dim(tree[["edge"]])[1] != 2L * tree[["Nnode"]]) {
      tree <- RootTree(tree, 1L)
    }
    if (dim(tree[["edge"]])[1] != 2L * tree[["Nnode"]]) {
      stop("Could not make `tree` binary.")
    }
  }

  # --- Match tree tips to dataset ---
  leaves <- tree[["tip.label"]]
  taxa <- names(dataset)
  treeOnly <- setdiff(leaves, taxa)
  datOnly <- setdiff(taxa, leaves)
  if (length(treeOnly)) {
    warning("Dropping taxa on tree but not in dataset: ",
            paste0(treeOnly, collapse = ", "))
    tree <- TreeTools::DropTip(tree, treeOnly)
  }
  if (length(datOnly)) {
    warning("Dropping taxa in dataset but not on tree: ",
            paste0(datOnly, collapse = ", "))
    dataset <- dataset[-match(datOnly, taxa)]
  }

  # Reorder tips to match dataset, put in preorder
  tree <- Preorder(RenumberTips(tree, names(dataset)))

  # Ensure root's first child is a tip (for C++ engine compatibility)
  if (tree[["edge"]][1L, 2L] > NTip(tree)) {
    tree <- RootTree(tree, 1L)
  }

  # --- Extract data matrices ---
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  weight <- at$weight
  levels <- at$levels

  # --- Run C++ driven search ---
  result <- ts_driven_search(
    contrast = contrast,
    tip_data = tip_data,
    weight = weight,
    levels = levels,
    maxReplicates = as.integer(maxReplicates),
    targetHits = as.integer(targetHits),
    tbrMaxHits = as.integer(tbrMaxHits),
    ratchetCycles = as.integer(ratchetCycles),
    ratchetPerturbProb = as.double(ratchetPerturbProb),
    xssRounds = as.integer(xssRounds),
    xssPartitions = as.integer(xssPartitions),
    sectorMinSize = as.integer(sectorMinSize),
    sectorMaxSize = as.integer(sectorMaxSize),
    fuseInterval = as.integer(fuseInterval),
    fuseAcceptEqual = as.logical(fuseAcceptEqual),
    poolMaxSize = as.integer(poolMaxSize),
    poolSuboptimal = as.double(poolSuboptimal)
  )

  # --- Reconstruct phylo from edge matrix ---
  bestTree <- tree
  bestTree[["edge"]] <- result$edge

  # --- Output ---
  if (verbosity > 0L) {
    cli_alert_success(paste0(
      "Search complete: score {.strong {signif(result$score, 7)}}, ",
      "{result$replicates} replicate{?s}, ",
      "{result$hits_to_best} hit{?s} to best"
    ))
  }

  structure(
    list(bestTree),
    score = result$score,
    replicates = result$replicates,
    hits_to_best = result$hits_to_best,
    class = "multiPhylo"
  )
}
