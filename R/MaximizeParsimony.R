# Internal helper: prepare constraint data for C++ engine.
# Returns a named list of constraint arguments (empty list if no constraint).
# @param constraint A phyDat, phylo, or NULL.
# @param dataset A phyDat whose names define the tip ordering.
# @keywords internal
.PrepareConstraint <- function(constraint, dataset) {
  if (is.null(constraint)) return(list())

  if (inherits(constraint, "phylo")) {
    constraint <- MatrixToPhyDat(t(as.matrix(constraint)))
  }
  if (!inherits(constraint, "phyDat")) {
    constraint <- MatrixToPhyDat(constraint)
  }

  # Match constraint taxa to dataset
  consTaxa <- names(constraint)
  treeTaxa <- names(dataset)
  treeOnly <- setdiff(treeTaxa, consTaxa)
  if (length(treeOnly)) {
    constraint <- AddUnconstrained(constraint, treeOnly)
  }
  consOnly <- setdiff(consTaxa, treeTaxa)
  if (length(consOnly)) {
    warning("Ignoring taxa in constraint missing on tree: ",
            paste0(consOnly, collapse = ", "))
    constraint <- constraint[-match(consOnly, consTaxa)]
  }
  constraint <- constraint[names(dataset)]

  consContrast <- attr(constraint, "contrast")
  nConsStates <- ncol(consContrast)
  if (nConsStates < 2L) return(list())

  consMat <- matrix(unlist(constraint, use.names = FALSE),
                    nrow = length(constraint), byrow = TRUE)
  consSplits <- matrix(0L, nrow = ncol(consMat), ncol = length(constraint))
  for (ch in seq_len(ncol(consMat))) {
    for (tip in seq_len(length(constraint))) {
      token <- consMat[tip, ch]
      if (consContrast[token, nConsStates] == 1 &&
          consContrast[token, 1] == 0) {
        consSplits[ch, tip] <- 1L
      }
    }
  }

  keep <- apply(consSplits, 1, function(row) {
    s <- sum(row)
    s >= 1 && s < length(constraint) - 1
  })
  consSplits <- consSplits[keep, , drop = FALSE]
  if (nrow(consSplits) == 0L) return(list())

  consWeight <- attr(constraint, "weight")
  consExpectedScore <- sum(
    MinimumLength(constraint, compress = TRUE) * consWeight
  )

  consTipData <- matrix(unlist(constraint, use.names = FALSE),
                        nrow = length(constraint), byrow = TRUE)

  list(
    consSplitMatrix = consSplits,
    consContrast = consContrast,
    consTipData = consTipData,
    consWeight = as.integer(consWeight),
    consLevels = attr(constraint, "levels"),
    consExpectedScore = as.integer(consExpectedScore)
  )
}

#' Find most parsimonious trees
#'
#' Performs a multi-replicate driven search for most-parsimonious trees,
#' combining random addition sequence (Wagner) starting trees, TBR
#' rearrangement, exclusive sectorial search (XSS), ratchet perturbation,
#' drift, and tree fusing -- all in compiled C++.
#'
#' Each replicate builds a random Wagner tree, optimizes it with TBR,
#' applies sectorial search and ratchet escape, then adds the result to a
#' pool of unique topologies.
#' Periodically, tree fusing recombines the best trees in the pool.
#' The search stops when the best score has been independently discovered
#' `targetHits` times, or `maxReplicates` replicates have been completed.
#'
#' Implied weighting is supported natively: set `concavity` to a numeric
#' value (e.g.\sspace{}10).
#' Profile parsimony (`concavity = "profile"`) is supported natively:
#' characters are simplified to binary (max 2 informative states),
#' inapplicable tokens are treated as ambiguous, and per-character
#' information profiles are used for scoring
#' \insertCite{Faith2001}{TreeSearch}.
#'
#' @param dataset A phylogenetic data matrix of \pkg{phangorn} class
#' \code{phyDat}, whose names correspond to the labels of any accompanying tree.
#' @param tree (optional) A bifurcating tree of class \code{\link[ape]{phylo}},
#'   containing only the tips listed in `dataset`, from which the search
#'   should begin.
#'   If unspecified, an [addition tree][AdditionTree()] will be generated.
#'   Edge lengths are not supported and will be deleted.
#' @param concavity Determines the degree to which extra steps beyond the first
#' are penalized.  Specify a numeric value to use implied weighting
#' \insertCite{Goloboff1993}{TreeSearch}; `concavity` specifies _k_ in
#'  _k_ / _e_ + _k_. A value of 10 is recommended;
#' TNT sets a default of 3, but this is too low in some circumstances
#' \insertCite{Goloboff2018,Smith2019}{TreeSearch}.
#' Better still explore the sensitivity of results under a range of
#' concavity values, e.g. `k = 2 ^ (1:7)`.
#' Specify `Inf` to weight each additional step equally.
#' Specify `"profile"` to employ profile parsimony
#' \insertCite{Faith2001}{TreeSearch}.
#' @param constraint Either an object of class `phyDat`, in which case
#' returned trees will be perfectly compatible with each character in
#' `constraint`; or a tree of class `phylo`, all of whose nodes will occur
#' in any output tree.
#' Constraint searches are supported natively: all tree rearrangements
#' are filtered to respect the constraint topology.
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates (default: 100).
#' @param targetHits Integer: stop when the best score has been found
#'   independently this many times (default: `max(10, NTip / 5)`).
#' @param tbrMaxHits Integer: within ratchet/drift TBR phases, the number
#'   of equal-score trees to accept before stopping.
#' @param ratchetCycles Integer: number of ratchet perturbation cycles
#'   per replicate.
#' @param ratchetPerturbProb Numeric: probability of zeroing each
#'   character's weight during ratchet perturbation.
#' @param ratchetPerturbMode Integer: perturbation mode. `0` = zero only
#'   (default), `1` = upweight only (double selected characters), `2` = mixed
#'   (zero some, double others).
#' @param ratchetPerturbMaxMoves Integer: maximum TBR moves during the
#'   perturbed search phase. `0` = automatic (default).
#' @param ratchetAdaptive Logical: if `TRUE`, automatically adjust
#'   perturbation probability based on escape rate.
#' @param driftCycles Integer: number of drift cycles per replicate
#'   (alternating suboptimal and equal-score phases).
#' @param driftAfdLimit Integer: maximum absolute fit difference (steps)
#'   for accepting suboptimal drift moves.
#' @param driftRfdLimit Numeric: maximum relative fit difference for
#'   accepting suboptimal drift moves.
#' @param xssRounds Integer: number of exclusive sectorial search rounds.
#' @param xssPartitions Integer: number of partitions for XSS.
#' @param rssRounds Integer: number of random sectorial search rounds after
#'   XSS. Set to `0` to disable RSS.
#' @param cssRounds Integer: number of constrained sectorial search (CSS)
#'   rounds. CSS runs sector-restricted TBR on the full tree (no HTU
#'   approximation), polishing improvements from XSS/RSS.
#'   Set to `0` to disable CSS.
#' @param cssPartitions Integer: number of partitions for CSS.
#' @param sectorMinSize,sectorMaxSize Integer: minimum and maximum clade
#'   sizes for sectorial search.
#' @param fuseInterval Integer: fuse pool trees every this many replicates.
#' @param fuseAcceptEqual Logical: accept equally-scoring fused trees?
#' @param poolMaxSize Integer: maximum number of trees retained in pool.
#' @param poolSuboptimal Numeric: retain trees within this many steps of
#'   the best score (0 = optimal only).
#' @param tabuSize Integer: size of the tabu list for TBR plateau exploration.
#'   Prevents cycling through recently visited topologies when accepting
#'   equal-score moves. Set to `0` to disable. Default 100.
#' @param wagnerStarts Integer: number of random Wagner starting trees per
#'   replicate. The best-scoring tree is kept before proceeding to TBR.
#'   Default 1.
#' @param nThreads Integer: number of parallel threads for search replicates.
#'   \describe{
#'     \item{`1` (default)}{Serial execution -- identical to previous behavior.}
#'     \item{`0`}{Auto-detect: use one fewer thread than the number of CPU
#'       cores.}
#'     \item{`> 1`}{Use the specified number of worker threads.}
#'   }
#'   In parallel mode, each replicate runs independently with a shared tree
#'   pool. Results may vary across runs with the same `set.seed()` due to
#'   thread scheduling nondeterminism. Use `nThreads = 1` for reproducible
#'   results.
#' @param verbosity Integer specifying level of messaging; higher values give
#' more detail. Set to `0` to run silently.
#' @param progressCallback Optional function called with a single list
#'   argument containing search progress information.
#'   The list includes elements: `replicate`, `max_replicates`,
#'   `best_score`, `hits_to_best`, `target_hits`, `pool_size`,
#'   `phase` (character), `elapsed` (seconds), and `phase_score`.
#'   When `NULL` (default) and `verbosity >= 1` in an interactive session,
#'   a `cli` progress bar is created automatically.
#'   Supply a custom function (e.g. using [shiny::setProgress()])
#'   to control progress display.
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
#' result <- MaximizeParsimony(dataset, maxReplicates = 3L, targetHits = 2L)
#' result
#' attr(result, "score")
#'
#' @template MRS
#' @family tree scoring
#' @seealso [`Morphy()`] for fine-grained control over the R-level search loop.
#' @references
#' \insertAllCited{}
#' @importFrom TreeTools NTip RenumberTips RootTree MakeTreeBinary Preorder
#' @importFrom cli cli_alert_success cli_alert_info cli_alert_warning
#' @encoding UTF-8
#' @export
MaximizeParsimony <- function(
    dataset,
    tree,
    concavity = Inf,
    constraint,
    maxReplicates = 100L,
    targetHits = max(10L, as.integer(NTip(dataset) / 5)),
    tbrMaxHits = 1L,
    ratchetCycles = 10L,
    ratchetPerturbProb = 0.04,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 0L,
    ratchetAdaptive = FALSE,
    driftCycles = 6L,
    driftAfdLimit = 3L,
    driftRfdLimit = 0.1,
    xssRounds = 3L,
    xssPartitions = 4L,
    rssRounds = 1L,
    cssRounds = 1L,
    cssPartitions = 4L,
    sectorMinSize = 6L,
    sectorMaxSize = 50L,
    fuseInterval = 3L,
    fuseAcceptEqual = FALSE,
    poolMaxSize = 100L,
    poolSuboptimal = 0.0,
    tabuSize = 100L,
    wagnerStarts = 1L,
    nThreads = 1L,
    verbosity = 1L,
    progressCallback = NULL
) {

  # --- Progress callback: build default cli bar if needed ---
  if (is.null(progressCallback) && verbosity >= 1L && interactive()) {
    pb_env <- new.env(parent = emptyenv())
    pb_env$id <- cli::cli_progress_bar(
      total = as.integer(maxReplicates),
      format = paste0(
        "Rep {cli::pb_current}/{cli::pb_total}",
        " | Best: {pb_env$best}",
        " | Hits: {pb_env$hits}/{pb_env$target}"
      ),
      .auto_close = FALSE,
      .envir = pb_env
    )
    pb_env$best <- "?"
    pb_env$hits <- 0L
    pb_env$target <- as.integer(targetHits)
    progressCallback <- function(info) {
      pb_env$best <- signif(info$best_score, 6)
      pb_env$hits <- info$hits_to_best
      pb_env$target <- info$target_hits
      if (identical(info$phase, "done")) {
        cli::cli_progress_done(id = pb_env$id, .envir = pb_env)
      } else if (identical(info$phase, "replicate")) {
        cli::cli_progress_update(
          id = pb_env$id, set = info$replicate, .envir = pb_env
        )
      }
    }
    on.exit(
      tryCatch(
        cli::cli_progress_done(id = pb_env$id, .envir = pb_env),
        error = function(e) NULL
      ),
      add = TRUE
    )
  }

  # --- Profile parsimony: prepare data ---
  useProfile <- !missing(concavity) && identical(concavity, "profile")
  if (useProfile) {
    dataset <- PrepareDataProfile(dataset)
    concavity <- Inf  # EW on the simplified binary data; profile scores via lookup
  }

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

  # --- Prepare constraint for C++ engine ---
  consArgs <- .PrepareConstraint(
    constraint = if (!missing(constraint)) constraint,
    dataset = dataset
  )
  if (length(consArgs) > 0L && verbosity > 0L) {
    cli_alert_info("Constraint: {nrow(consArgs$consSplitMatrix)} split{?s}")
  }

  # --- Profile parsimony: extract info_amounts ---
  profileArgs <- list()
  if (useProfile) {
    infoAmounts <- attr(dataset, "info.amounts")
    if (!is.null(infoAmounts) && length(infoAmounts) > 0L) {
      profileArgs$infoAmounts <- infoAmounts
    }
  }

  # --- Run C++ driven search ---
  searchArgs <- list(
    contrast = contrast,
    tip_data = tip_data,
    weight = weight,
    levels = levels,
    maxReplicates = as.integer(maxReplicates),
    targetHits = as.integer(targetHits),
    tbrMaxHits = as.integer(tbrMaxHits),
    ratchetCycles = as.integer(ratchetCycles),
    ratchetPerturbProb = as.double(ratchetPerturbProb),
    ratchetPerturbMode = as.integer(ratchetPerturbMode),
    ratchetPerturbMaxMoves = as.integer(ratchetPerturbMaxMoves),
    ratchetAdaptive = as.logical(ratchetAdaptive),
    driftCycles = as.integer(driftCycles),
    driftAfdLimit = as.integer(driftAfdLimit),
    driftRfdLimit = as.double(driftRfdLimit),
    xssRounds = as.integer(xssRounds),
    xssPartitions = as.integer(xssPartitions),
    rssRounds = as.integer(rssRounds),
    cssRounds = as.integer(cssRounds),
    cssPartitions = as.integer(cssPartitions),
    sectorMinSize = as.integer(sectorMinSize),
    sectorMaxSize = as.integer(sectorMaxSize),
    fuseInterval = as.integer(fuseInterval),
    fuseAcceptEqual = as.logical(fuseAcceptEqual),
    poolMaxSize = as.integer(poolMaxSize),
    poolSuboptimal = as.double(poolSuboptimal),
    tabuSize = as.integer(tabuSize),
    wagnerStarts = as.integer(wagnerStarts),
    verbosity = as.integer(verbosity),
    concavity = as.double(concavity),
    progressCallback = progressCallback,
    nThreads = as.integer(nThreads)
  )
  result <- do.call(ts_driven_search, c(searchArgs, consArgs, profileArgs))

  # --- Reconstruct phylo from edge matrices ---
  treeTpl <- tree
  treeTpl[["edge.length"]] <- NULL
  resultTrees <- result$trees
  if (length(resultTrees) == 0L) {
    resultTrees <- list()
  }
  outTrees <- lapply(resultTrees, function(edgeMat) {
    tr <- treeTpl
    tr[["edge"]] <- edgeMat
    tr
  })
  if (length(outTrees) == 0L) {
    outTrees <- list(treeTpl)
  }

  # --- Output ---
  if (verbosity > 0L) {
    cli_alert_success(paste0(
      "Search complete: score {.strong {signif(result$best_score, 7)}}, ",
      "{result$replicates} replicate{?s}, ",
      "{result$hits_to_best} hit{?s} to best"
    ))
  }

  structure(
    outTrees,
    score = result$best_score,
    replicates = result$replicates,
    hits_to_best = result$hits_to_best,
    class = "multiPhylo"
  )
}


#' @rdname MaximizeParsimony
#' @usage MaximizeParsimony2(...)
#' @section Deprecated:
#' `MaximizeParsimony2()` is a deprecated alias for `MaximizeParsimony()`.
#' @export
MaximizeParsimony2 <- function(...) {
  .Deprecated("MaximizeParsimony")
  MaximizeParsimony(...)
}
