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

# Strategy presets for adaptive search (Phase 6E).
# Wrapped in a function to avoid load-order dependency on SearchControl().
.StrategyPresets <- function() list(
  sprint = SearchControl(
    tbrMaxHits = 1L, ratchetCycles = 3L, ratchetPerturbProb = 0.04,
    ratchetPerturbMode = 0L, ratchetAdaptive = FALSE,
    driftCycles = 0L, xssRounds = 1L, xssPartitions = 4L,
    rssRounds = 0L, cssRounds = 0L, cssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 5L, fuseAcceptEqual = FALSE,
    tabuSize = 0L, wagnerStarts = 1L
  ),
  default = SearchControl(
    tbrMaxHits = 1L, ratchetCycles = 5L, ratchetPerturbProb = 0.04,
    ratchetPerturbMode = 0L, ratchetAdaptive = FALSE,
    driftCycles = 2L, xssRounds = 3L, xssPartitions = 4L,
    rssRounds = 1L, cssRounds = 0L, cssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 3L, fuseAcceptEqual = FALSE,
    tabuSize = 100L, wagnerStarts = 1L
  ),
  thorough = SearchControl(
    tbrMaxHits = 3L, ratchetCycles = 20L, ratchetPerturbProb = 0.04,
    ratchetPerturbMode = 2L, ratchetAdaptive = TRUE,
    driftCycles = 12L, driftAfdLimit = 5L, driftRfdLimit = 0.15,
    xssRounds = 5L, xssPartitions = 6L,
    rssRounds = 3L, cssRounds = 2L, cssPartitions = 6L,
    sectorMinSize = 6L, sectorMaxSize = 80L,
    fuseInterval = 2L, fuseAcceptEqual = TRUE,
    tabuSize = 200L, wagnerStarts = 3L
  )
)

# Select strategy preset based on dataset size and character count.
# @param nTip Integer number of taxa
# @param nChar Integer number of character patterns (unique columns)
# @return Character name of the strategy preset
# @details
# Empirically calibrated on 15 neotrans matrices (61-86 tips) + 4
# inapplicable.phyData datasets.  Key findings:
#   - Datasets with few characters (< 100 patterns) have flat parsimony
#     landscapes where extra search adds zero score improvement (0/6 benefited).
#   - Datasets with >= 100 patterns and >= 65 taxa have structured landscapes
#     where thorough search finds substantially better trees (7/9 benefited,
#     median +14 steps, max +74 steps at 86 tips / 528 chars).
#   - At 62 tips (Agnarsson2004, 242 patterns) thorough adds 0 steps; at 65
#     tips (project3617, 361 patterns) it adds 14 steps.
.AutoStrategy <- function(nTip, nChar) {
  if (nTip <= 30L) return("sprint")
  # Few characters -> flat landscape; thorough search is pointless
  if (nChar < 100L) return("default")
  # Enough characters to have a structured landscape;
  # moderate-to-large datasets benefit from intensive search
  if (nTip >= 65L) return("thorough")
  "default"
}

#' Find most parsimonious trees
#'
#' Performs a multi-replicate driven search for most-parsimonious trees,
#' combining random addition sequence (Wagner) starting trees, TBR
#' rearrangement, exclusive sectorial search (XSS), ratchet perturbation,
#' drift, and tree fusing -- all in compiled C++.
#'
#' The search pipeline follows the "new technology search" approach of
#' \insertCite{Goloboff1999;textual}{TreeSearch}, as implemented in TNT
#' \insertCite{Goloboff2016}{TreeSearch}.
#' Parsimony scoring uses the Fitch
#' \insertCite{Fitch1971}{TreeSearch} algorithm; inapplicable characters
#' are handled with the algorithm of
#' \insertCite{Brazeau2019;textual}{TreeSearch}.
#' Each replicate builds a random addition sequence (Wagner) tree
#' \insertCite{Kluge1969}{TreeSearch}, optimizes it with TBR,
#' applies sectorial search and the parsimony ratchet
#' \insertCite{Nixon1999}{TreeSearch} to escape local optima, then adds
#' the result to a pool of unique topologies.
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
#'   or a `multiPhylo` (first tree used).
#'   When supplied, the first replicate uses this topology as its starting
#'   point (warm-start), skipping the random Wagner tree construction.
#'   Subsequent replicates still begin from random Wagner trees.
#'   This is useful for continuing a search from a previously found optimum.
#'   If unspecified, all replicates start from random Wagner trees.
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
#' @param strategy Character: named strategy preset controlling the search
#'   heuristic parameters. Presets:
#'   \describe{
#'     \item{`"auto"` (default)}{Selects automatically based on dataset size
#'       and signal density (characters per taxon):
#'       `"sprint"` for <=30 taxa, `"thorough"` for >=75 taxa with
#'       fewer than 5 characters per taxon, `"default"` otherwise.}
#'     \item{`"sprint"`}{Fast search: 3 ratchet cycles, no drift, minimal
#'       sectorial. Good for small datasets or quick surveys.}
#'     \item{`"default"`}{Balanced: 5 ratchet + 2 drift + sectorial + fusing.}
#'     \item{`"thorough"`}{Intensive: 20 ratchet cycles, 12 drift, adaptive
#'       perturbation, extra sectorial rounds. Best for large (75+ tips)
#'       datasets with limited phylogenetic signal, where sprint/default
#'       may miss the global optimum.}
#'     \item{`"none"`}{Use only the explicitly supplied parameter values.}
#'   }
#'   Explicit `control` fields always override the preset; for example,
#'   `strategy = "sprint", control = SearchControl(ratchetCycles = 10L)` uses
#'   sprint defaults for everything except `ratchetCycles`.
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates (default: 100).
#'   For large or complex datasets a higher value improves the chance of
#'   finding all MPTs.  A rough minimum is
#'   `max(10, ceiling(NTip * NChar / 5000))`, where `NChar = sum(weight)`.
#'   A warning is issued when an explicit value falls below this threshold
#'   for datasets with 30 or more taxa.
#' @param targetHits Integer: stop when the best score has been found
#'   independently this many times (default: `max(10, NTip / 5)`).
#' @param maxSeconds Numeric: maximum wall-clock time in seconds for the
#'   search. When reached, the current replicate finishes and the search
#'   stops. `0` (default) means no time limit.
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
#' @param control A [`SearchControl`] object (or a named list) of low-level
#'   search parameters.  Most users can rely on the `strategy` presets and
#'   ignore this argument; see [`SearchControl()`] for full documentation
#'   of individual fields.
#' @param ... Backward compatibility: individual control parameters (e.g.
#'   `ratchetCycles = 10L`) may still be passed as named arguments.
#'   These override the corresponding `control` fields and the strategy
#'   preset.
#'   Legacy `Morphy()`-style parameters (e.g. `ratchIter`, `tbrIter`) are
#'   detected and forwarded to [`Morphy()`] with a deprecation warning.
#'
#' @return A `multiPhylo` object containing the best tree(s) found, with
#'   attributes:
#'   \describe{
#'     \item{`score`}{Best parsimony score.}
#'     \item{`replicates`}{Number of replicates completed.}
#'     \item{`hits_to_best`}{Number of independent discoveries of the best
#'       score.}
#'     \item{`timed_out`}{Logical: `TRUE` if the search stopped because
#'       `maxSeconds` was exceeded.}
#'     \item{`timings`}{Named numeric vector of cumulative wall-clock time
#'       (in milliseconds) spent in each search phase across all replicates:
#'       `wagner_ms`, `tbr_ms`, `xss_ms`, `rss_ms`, `css_ms`, `ratchet_ms`,
#'       `drift_ms`, `final_tbr_ms`, `fuse_ms`.}
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
#' [`Resample()`] for jackknife and bootstrap resampling.
#' [`SearchControl()`] for expert-level tuning of the search heuristics.
#' @references
#' \insertAllCited{}
#' @importFrom TreeTools NTip RandomTree Renumber RenumberTips RootTree MakeTreeBinary
#'   Preorder
#' @importFrom cli cli_alert_success cli_alert_info cli_alert_warning
#' @encoding UTF-8
#' @export
MaximizeParsimony <- function(
    dataset,
    tree,
    concavity = Inf,
    constraint,
    strategy = "auto",
    maxReplicates = 100L,
    targetHits = max(10L, as.integer(NTip(dataset) / 5)),
    maxSeconds = 0,
    nThreads = 1L,
    verbosity = 1L,
    progressCallback = NULL,
    control = SearchControl(),
    ...
) {

  # --- Backward compatibility: detect Morphy()-style parameters ---
  dots <- list(...)
  .morphyParams <- c("ratchIter", "tbrIter", "startIter", "finalIter",
                      "maxHits", "maxTime", "quickHits", "ratchEW",
                      "tolerance")
  legacyHits <- intersect(names(dots), .morphyParams)
  if (length(legacyHits)) {
    .Deprecated(
      "Morphy",
      msg = paste0(
        "Parameter", if (length(legacyHits) > 1L) "s", " ",
        paste0(sQuote(legacyHits), collapse = ", "),
        " belong", if (length(legacyHits) == 1L) "s", " to `Morphy()`,",
        " not the new `MaximizeParsimony()`.\n",
        "  Delegating to `Morphy()`. ",
        "Please update your code to call `Morphy()` directly ",
        "or use the new MaximizeParsimony() parameters.\n",
        "  See ?Morphy and ?MaximizeParsimony for details."
      )
    )
    morphyArgs <- dots
    morphyArgs$dataset <- dataset
    if (!missing(tree)) morphyArgs$tree <- tree
    if (!missing(concavity)) morphyArgs$concavity <- concavity
    if (!missing(constraint)) morphyArgs$constraint <- constraint
    if (!missing(verbosity)) morphyArgs$verbosity <- verbosity
    return(do.call(Morphy, morphyArgs))
  }

  # --- Resolve control: merge control + ... overrides ---
  # Coerce a plain list to SearchControl
  if (!inherits(control, "SearchControl")) {
    control <- do.call(SearchControl, control)
  }

  # Named ... args that match SearchControl fields override `control`
  controlFields <- names(SearchControl())
  controlDots <- dots[intersect(names(dots), controlFields)]
  otherDots <- dots[setdiff(names(dots), controlFields)]
  if (length(controlDots)) {
    for (nm in names(controlDots)) {
      control[[nm]] <- controlDots[[nm]]
    }
  }
  if (length(otherDots)) {
    warning("Unknown arguments ignored: ",
            paste0(sQuote(names(otherDots)), collapse = ", "))
  }

  # --- Apply strategy preset ---
  if (!is.null(strategy) && !identical(strategy, "none")) {
    if (identical(strategy, "auto")) {
      strategy <- .AutoStrategy(NTip(dataset),
                                sum(attr(dataset, "weight")))
    }
    preset <- .StrategyPresets()[[strategy]]
    if (!is.null(preset)) {
      # Determine which control fields the user explicitly set.
      # Fields are "explicit" if:
      #   (a) passed via ... (already merged into control above), OR
      #   (b) control was explicitly supplied and differs from SearchControl()
      defaults <- SearchControl()
      explicit_via_dots <- names(controlDots)
      explicit_via_control <- if ("control" %in% names(match.call())) {
        # User passed control = SearchControl(...) — honour all fields in it
        names(control)
      } else {
        character(0)
      }
      explicit <- union(explicit_via_dots, explicit_via_control)

      # Apply preset values for any field the user didn't explicitly set
      for (nm in names(preset)) {
        if (!(nm %in% explicit)) {
          control[[nm]] <- preset[[nm]]
        }
      }
      if (verbosity >= 1L) {
        cli::cli_alert_info("Strategy: {.strong {strategy}}")
      }
    } else if (!identical(strategy, "auto")) {
      warning("Unknown strategy '", strategy, "'; using default parameters.")
    }
  }

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
  if (is.null(attr(dataset, "levels")) || ncol(attr(dataset, "contrast")) == 0L) {
    stop("Dataset contains no informative character states.")
  }

  # --- Starting tree ---
  userTree <- !missing(tree)
  if (!userTree) {
    tree <- TreeTools::RandomTree(nTip, root = TRUE)
    tree[["tip.label"]] <- names(dataset)
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

  # --- Replicate count adequacy check ---
  # Warn only when the user explicitly passed maxReplicates.
  # Formula: max(10, ceiling(nTip * nChar / 5000)) where nChar = sum(weight).
  # Derived from T-069 benchmarks: at 225 taxa / 748 chars a single rep takes
  # ~40s and at least ~34 reps are needed to fill the tree pool reliably.
  if (!missing(maxReplicates) && nTip >= 30L) {
    nChars <- sum(weight)
    minReps <- pmax(10L, ceiling(nTip * nChars / 5000L))
    if (maxReplicates < minReps) {
      warning(
        "With ", nTip, " taxa and ", nChars, " characters, at least ",
        minReps, " replicates are recommended for reliable results ",
        "(you specified ", maxReplicates, "). ",
        "Consider increasing `maxReplicates` or setting `maxSeconds` ",
        "to allow more search time.",
        call. = FALSE
      )
    }
  }

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

  # --- IW: compute minimum step counts per character ---
  if (is.finite(concavity)) {
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
  }

  # --- Run C++ driven search ---
  # Read all control fields from the resolved control object
  ctrl <- control
  searchArgs <- list(
    contrast = contrast,
    tip_data = tip_data,
    weight = weight,
    levels = levels,
    maxReplicates = as.integer(maxReplicates),
    targetHits = as.integer(targetHits),
    tbrMaxHits = as.integer(ctrl$tbrMaxHits),
    ratchetCycles = as.integer(ctrl$ratchetCycles),
    ratchetPerturbProb = as.double(ctrl$ratchetPerturbProb),
    ratchetPerturbMode = as.integer(ctrl$ratchetPerturbMode),
    ratchetPerturbMaxMoves = as.integer(ctrl$ratchetPerturbMaxMoves),
    ratchetAdaptive = as.logical(ctrl$ratchetAdaptive),
    driftCycles = as.integer(ctrl$driftCycles),
    driftAfdLimit = as.integer(ctrl$driftAfdLimit),
    driftRfdLimit = as.double(ctrl$driftRfdLimit),
    xssRounds = as.integer(ctrl$xssRounds),
    xssPartitions = as.integer(ctrl$xssPartitions),
    rssRounds = as.integer(ctrl$rssRounds),
    cssRounds = as.integer(ctrl$cssRounds),
    cssPartitions = as.integer(ctrl$cssPartitions),
    sectorMinSize = as.integer(ctrl$sectorMinSize),
    sectorMaxSize = as.integer(ctrl$sectorMaxSize),
    fuseInterval = as.integer(ctrl$fuseInterval),
    fuseAcceptEqual = as.logical(ctrl$fuseAcceptEqual),
    poolMaxSize = as.integer(ctrl$poolMaxSize),
    poolSuboptimal = as.double(ctrl$poolSuboptimal),
    maxSeconds = as.double(maxSeconds),
    tabuSize = as.integer(ctrl$tabuSize),
    wagnerStarts = as.integer(ctrl$wagnerStarts),
    verbosity = as.integer(verbosity),
    min_steps = if (is.finite(concavity)) minSteps else integer(0),
    concavity = as.double(concavity),
    progressCallback = progressCallback,
    nThreads = as.integer(nThreads),
    startEdge = if (userTree) tree[["edge"]] else NULL,
    sprFirst = as.logical(ctrl$sprFirst)
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
    # C++ edge order may differ from template; renumber to valid preorder
    Renumber(tr)
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
    timed_out = isTRUE(result$timed_out),
    timings = unlist(result$timings),
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
