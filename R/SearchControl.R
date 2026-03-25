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
#' @param nniFirst Logical; run an NNI pass before SPR/TBR in each replicate?
#'   At small tree sizes (\eqn{\le}88 tips) overhead is negligible; at \eqn{\ge}100 tips
#'   this significantly accelerates the initial descent from the Wagner tree.
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
#'   within-replicate escape rate?
#' @param ratchetTaper Logical; taper ratchet perturbation probability across
#'   replicates as the pool stabilizes?  When `TRUE`, early replicates use
#'   the full `ratchetPerturbProb`; later replicates (with high hit rates)
#'   use a reduced probability for finer local exploration.  The effective
#'   probability is `ratchetPerturbProb * max(floor, 1 - strength * hitRate)`
#'   where `hitRate` is the fraction of replicates that found the current
#'   best score.  Default `FALSE`.
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
#' @param perturbStopFactor Integer; stop after
#'   `nTip * perturbStopFactor` consecutive replicates that fail to improve
#'   the best score.  0 (default) disables this criterion.
#'   Inspired by IQ-TREE's unsuccessful-perturbation stopping rule
#'   \insertCite{Nguyen2015}{TreeSearch}; adapted from per-perturbation to
#'   per-replicate granularity.  Small values (1--3) are typical.
#' @param adaptiveLevel Logical; dynamically scale ratchet and drift effort
#'   based on the observed hit rate?  When `TRUE`, easy landscapes
#'   (high hit rate) trigger reduced effort per replicate, while hard
#'   landscapes trigger increased effort.  Default `FALSE`.
#' @param nniPerturbCycles Integer; number of stochastic NNI-perturbation
#'   cycles per replicate.  Each cycle randomly applies NNI swaps to a
#'   fraction of internal branches, then runs TBR to find a new local
#'   optimum.  Complementary to the weight-perturbation ratchet: the ratchet
#'   perturbs the objective function, while NNI-perturbation perturbs the
#'   topology directly.
#'   0 (default) disables NNI perturbation.
#'   Inspired by `doRandomNNIs()` in IQ-TREE
#'   \insertCite{Nguyen2015}{TreeSearch}.
#' @param nniPerturbFraction Numeric (0--1); fraction of internal branches
#'   to swap during each NNI-perturbation cycle.  Default 0.5.
#' @param consensusConstrain Logical; lock the strict consensus of pool
#'   trees as topological constraints for subsequent replicates?  When
#'   `TRUE`, after enough replicates (\eqn{\ge}5), splits present in ALL
#'   best-score pool trees are enforced as constraints, focusing search on
#'   uncertain regions.  Constraints are cleared whenever a new best score
#'   is found.  Only active when no user-supplied `constraint` is
#'   present.  Default `FALSE`.
#' @param wagnerBias Integer; criterion for biasing taxon addition order
#'   during Wagner tree construction.  0 = random (default),
#'   1 = Goloboff (2014) non-ambiguous-character priority,
#'   2 = entropy-based state-specificity priority.  Biased orders use
#'   softmax-weighted sampling for diversity across replicates.
#' @param wagnerBiasTemp Numeric; softmax temperature controlling
#'   selectivity of biased Wagner addition (default 0.3).  Lower values
#'   concentrate sampling on the highest-scoring taxa; higher values
#'   approach uniform random.
#' @param outerCycles Integer; number of outer search cycles per replicate
#'   (default 1).  Each outer cycle runs the full
#'   \[XSS/RSS/CSS → ratchet → NNI-perturbation → drift → TBR\] sequence,
#'   with perturbation cycles divided evenly among outer iterations.
#'   Matches the interleaved sectorial + ratchet pattern of TNT's `xmult`
#'   \insertCite{Goloboff1999}{TreeSearch}.
#' @param maxOuterResets Integer; maximum number of improvement-triggered
#'   resets of the outer cycle counter (default 0 = no resets, so
#'   `outerCycles` is exact).  When the search finds a new best score during
#'   an outer cycle, the counter resets up to this many times, allowing
#'   productive re-exploration.  Set to \eqn{-1} for unlimited resets.
#'   Strategy presets (`"default"`, `"thorough"`) set 2–3.
#' @param annealCycles Integer; number of simulated annealing perturbation
#'   cycles (PCSA) per replicate.  Each cycle perturbs the current best tree
#'   via scheduled SA cooling, then reconverges with TBR.  If the result
#'   improves on the best, it becomes the new starting point.  Effective at
#'   escaping deep basins under equal-weights parsimony at \eqn{\ge}100 tips.
#'   0 (default) disables SA perturbation.
#' @param annealPhases Integer; number of temperature steps in the linear
#'   cooling schedule per SA cycle (default 5).
#' @param annealTStart Numeric; initial Boltzmann temperature for SA cooling
#'   schedule (default 20).  Higher temperatures accept more suboptimal moves.
#' @param annealTEnd Numeric; final Boltzmann temperature (default 0 =
#'   strict hill-climbing at end of each cycle).
#' @param annealMovesPerPhase Integer; stochastic TBR moves per temperature
#'   step (default 0 = number of tips).
#' @param enumTimeFraction Numeric between 0 and 0.5; fraction of `maxSeconds`
#'   reserved for MPT enumeration (TBR plateau walk to discover additional
#'   equal-score topologies).  The main search loop exits at
#'   `maxSeconds * (1 - enumTimeFraction)`.  Set to 0 to disable the reserve
#'   (pre-v1.6 behaviour: enumeration skipped if the main loop times out).
#'   Default: `0.1` (10%).
#' @param adaptiveStart Logical; use Thompson-sampling (bandit) strategy
#'   selection for starting trees?  When `TRUE`, each replicate draws its
#'   starting strategy from a pool of options (random Wagner, biased Wagner,
#'   random tree, pool ratchet, pool NNI-perturb), adapting to which
#'   strategies yield the best scores.  Default `FALSE`.
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
    nniFirst = TRUE,
    sprFirst = FALSE,
    tabuSize = 100L,
    wagnerStarts = 1L,
    # Wagner biased addition (Goloboff 2014 §3.3)
    # 0L = random (default), 1L = Goloboff non-ambiguous score, 2L = entropy
    wagnerBias = 0L,
    wagnerBiasTemp = 0.3,
    # Outer search cycle count (Goloboff 1999 §2.3)
    # Repeat [XSS → Ratchet → NNI-perturb → Drift → TBR] this many times.
    # Cycles are divided evenly; default 1 = single pipeline pass.
    outerCycles = 1L,
    # Max improvement-triggered resets of the outer cycle counter.
    # 0 = no resets (outerCycles is exact); -1 = unlimited.
    # Strategy presets set 2-3 for productive re-exploration.
    maxOuterResets = 0L,
    # Ratchet
    ratchetCycles = 12L,
    ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = FALSE,
    ratchetTaper = FALSE,
    # NNI perturbation
    nniPerturbCycles = 0L,
    nniPerturbFraction = 0.5,
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
    perturbStopFactor = 0L,
    adaptiveLevel = FALSE,
    consensusConstrain = FALSE,
    # Simulated annealing perturbation (PCSA, T-207)
    annealCycles = 0L,
    annealPhases = 5L,
    annealTStart = 20,
    annealTEnd = 0,
    annealMovesPerPhase = 0L,
    # Adaptive starting-tree strategy (T-190)
    # When TRUE, each replicate draws its starting strategy via Thompson
    # sampling from {Wagner-random, Wagner-Goloboff, Wagner-entropy,
    # random-tree, pool-ratchet, pool-NNI-perturb}. Overrides wagnerBias.
    adaptiveStart = FALSE,
    enumTimeFraction = 0.1
) {
  structure(
    list(
      tbrMaxHits = as.integer(tbrMaxHits),
      nniFirst = as.logical(nniFirst),
      sprFirst = as.logical(sprFirst),
      tabuSize = as.integer(tabuSize),
      wagnerStarts = as.integer(wagnerStarts),
      wagnerBias = as.integer(wagnerBias),
      wagnerBiasTemp = as.double(wagnerBiasTemp),
      outerCycles = as.integer(outerCycles),
      maxOuterResets = as.integer(maxOuterResets),
      ratchetCycles = as.integer(ratchetCycles),
      ratchetPerturbProb = as.double(ratchetPerturbProb),
      ratchetPerturbMode = as.integer(ratchetPerturbMode),
      ratchetPerturbMaxMoves = as.integer(ratchetPerturbMaxMoves),
      ratchetAdaptive = as.logical(ratchetAdaptive),
      ratchetTaper = as.logical(ratchetTaper),
      nniPerturbCycles = as.integer(nniPerturbCycles),
      nniPerturbFraction = as.double(nniPerturbFraction),
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
      perturbStopFactor = as.integer(perturbStopFactor),
      adaptiveLevel = as.logical(adaptiveLevel),
      consensusConstrain = as.logical(consensusConstrain),
      annealCycles = as.integer(annealCycles),
      annealPhases = as.integer(annealPhases),
      annealTStart = as.double(annealTStart),
      annealTEnd = as.double(annealTEnd),
      annealMovesPerPhase = as.integer(annealMovesPerPhase),
      adaptiveStart = as.logical(adaptiveStart),
      enumTimeFraction = as.double(enumTimeFraction)
    ),
    class = "SearchControl"
  )
}

#' @export
print.SearchControl <- function(x, ...) {
  groups <- list(
    "TBR" = c("tbrMaxHits", "nniFirst", "sprFirst", "tabuSize",
              "wagnerStarts", "wagnerBias", "wagnerBiasTemp", "outerCycles"),
    "Ratchet" = c("ratchetCycles", "ratchetPerturbProb", "ratchetPerturbMode",
                   "ratchetPerturbMaxMoves", "ratchetAdaptive",
                   "ratchetTaper"),
    "NNI Perturbation" = c("nniPerturbCycles", "nniPerturbFraction"),
    "Drift" = c("driftCycles", "driftAfdLimit", "driftRfdLimit"),
    "Annealing" = c("annealCycles", "annealPhases", "annealTStart",
                     "annealTEnd", "annealMovesPerPhase"),
    "Sectorial" = c("xssRounds", "xssPartitions", "rssRounds",
                     "cssRounds", "cssPartitions",
                     "sectorMinSize", "sectorMaxSize"),
    "Fuse/Pool" = c("fuseInterval", "fuseAcceptEqual",
                     "poolMaxSize", "poolSuboptimal"),
    "Stopping" = c("consensusStableReps", "perturbStopFactor",
                    "adaptiveLevel",
                    "consensusConstrain", "adaptiveStart")
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
