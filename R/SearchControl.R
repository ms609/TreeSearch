#' Expert search heuristic parameters
#'
#' Construct a list of low-level search parameters for
#' [`MaximizeParsimony()`].  Most users can ignore these and rely on the
#' `strategy` presets (`"sprint"`, `"default"`, `"thorough"`); `SearchControl`
#' is provided for expert tuning.
#'
#' The parameters correspond to heuristics described by
#' \insertCite{Goloboff1999;textual}{TreeSearch}
#' (sectorial search, tree drifting, tree fusing) and
#' \insertCite{Nixon1999;textual}{TreeSearch}
#' (parsimony ratchet), as implemented in TNT
#' \insertCite{Goloboff2016}{TreeSearch}.
#'
#' @param tbrMaxHits Integer; number of equally-scoring trees to accept
#'   before stopping a TBR pass.
#' @param sprFirst Logical; run an SPR pass before TBR in each replicate?
#' @param tabuSize Integer; tabu list size for TBR plateau exploration.
#' @param wagnerStarts Integer; random Wagner starting trees per replicate.
#' @param ratchetCycles Integer; number of ratchet perturbation cycles.
#' @param ratchetPerturbProb Numeric (0--1); probability of perturbing each
#'   character.
#' @param ratchetPerturbMode Integer; 0 = zero-weight only, 1 = up-weight only,
#'   2 = mixed.
#' @param ratchetPerturbMaxMoves Integer; maximum TBR moves per perturbation
#'   cycle (0 = automatic).
#' @param ratchetAdaptive Logical; adjust perturbation probability based on
#'   escape rate?
#' @param driftCycles Integer; number of drift search cycles.
#' @param driftAfdLimit Integer; maximum absolute fit difference (steps) for
#'   accepting a suboptimal drift move.
#' @param driftRfdLimit Numeric; maximum relative fit difference for
#'   accepting a suboptimal drift move.
#' @param xssRounds Integer; rounds of exclusive sectorial search.
#' @param xssPartitions Integer; number of partitions in XSS.
#' @param rssRounds Integer; rounds of random sectorial search.
#' @param cssRounds Integer; rounds of constrained (sector-restricted TBR)
#'   sectorial search.
#' @param cssPartitions Integer; number of partitions in CSS.
#' @param sectorMinSize,sectorMaxSize Integer; minimum and maximum clade
#'   sizes for sectorial search.
#' @param fuseInterval Integer; fuse pool trees every _n_ replicates.
#' @param fuseAcceptEqual Logical; accept equally-scoring fused trees?
#' @param poolMaxSize Integer; maximum trees retained in the pool.
#' @param poolSuboptimal Numeric; retain trees that are this many steps
#'   worse than the best tree.  0 (default) keeps only optimal trees.
#' @param consensusStableReps Integer; stop when the strict consensus of
#'   best-score pool trees has been unchanged for this many consecutive
#'   replicates.
#'   0 (default) disables this criterion; a typical value is 3--5.
#'   When both `consensusStableReps` and `targetHits` are active, the search
#'   stops when either criterion is met first.
#' @param adaptiveLevel Logical; dynamically scale ratchet and drift effort
#'   based on the observed hit rate?  When `TRUE`, easy landscapes
#'   (high hit rate) trigger reduced effort per replicate, while hard
#'   landscapes trigger increased effort.  Default `FALSE`.
#' @param consensusConstrain Logical; lock the strict consensus of pool
#'   trees as topological constraints for subsequent replicates?  When
#'   `TRUE`, after enough replicates (≥5), splits present in ALL
#'   best-score pool trees are enforced as constraints, focusing search on
#'   uncertain regions.  Constraints are cleared whenever a new best score
#'   is found.  Only active when no user-supplied `constraint` is
#'   present.  Default `FALSE`.
#'
#' @return A named list of class `"SearchControl"`.
#'
#' @examples
#' # Use defaults
#' SearchControl()
#'
#' # Light ratchet, no drift
#' SearchControl(ratchetCycles = 5L, ratchetPerturbProb = 0.04,
#'               driftCycles = 0L)
#'
#' @family tree search functions
#' @seealso [`MaximizeParsimony()`]
#' @references
#' \insertAllCited{}
#' @export
SearchControl <- function(
    # TBR
    tbrMaxHits = 1L,
    sprFirst = FALSE,
    tabuSize = 100L,
    wagnerStarts = 1L,
    # Ratchet
    ratchetCycles = 12L,
    ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = FALSE,
    # Drift
    driftCycles = 2L,
    driftAfdLimit = 5L,
    driftRfdLimit = 0.15,
    # Sectorial
    xssRounds = 3L,
    xssPartitions = 4L,
    rssRounds = 1L,
    cssRounds = 0L,
    cssPartitions = 4L,
    sectorMinSize = 6L,
    sectorMaxSize = 50L,
    # Fuse / pool
    fuseInterval = 3L,
    fuseAcceptEqual = FALSE,
    poolMaxSize = 100L,
    poolSuboptimal = 0,
    # Stopping criteria
    consensusStableReps = 0L,
    adaptiveLevel = FALSE,
    consensusConstrain = FALSE
) {
  structure(
    list(
      tbrMaxHits = as.integer(tbrMaxHits),
      sprFirst = as.logical(sprFirst),
      tabuSize = as.integer(tabuSize),
      wagnerStarts = as.integer(wagnerStarts),
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
      consensusStableReps = as.integer(consensusStableReps),
      adaptiveLevel = as.logical(adaptiveLevel),
      consensusConstrain = as.logical(consensusConstrain)
    ),
    class = "SearchControl"
  )
}

#' @export
print.SearchControl <- function(x, ...) {
  groups <- list(
    "TBR" = c("tbrMaxHits", "sprFirst", "tabuSize", "wagnerStarts"),
    "Ratchet" = c("ratchetCycles", "ratchetPerturbProb", "ratchetPerturbMode",
                   "ratchetPerturbMaxMoves", "ratchetAdaptive"),
    "Drift" = c("driftCycles", "driftAfdLimit", "driftRfdLimit"),
    "Sectorial" = c("xssRounds", "xssPartitions", "rssRounds",
                     "cssRounds", "cssPartitions",
                     "sectorMinSize", "sectorMaxSize"),
    "Fuse/Pool" = c("fuseInterval", "fuseAcceptEqual",
                     "poolMaxSize", "poolSuboptimal"),
    "Stopping" = c("consensusStableReps", "adaptiveLevel",
                    "consensusConstrain")
  )
  cat("SearchControl object\n")
  for (gname in names(groups)) {
    cat(sprintf("  %s:\n", gname))
    for (pname in groups[[gname]]) {
      cat(sprintf("    %-25s = %s\n", pname, format(x[[pname]])))
    }
  }
  invisible(x)
}
