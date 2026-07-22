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
#' @param clipOrder Integer (experimental); clip-ordering strategy for TBR
#'   search.  Determines the order in which edges are tried as clip points.
#'   0 = random (default); 1 = inverse-weight (fewest descendant taxa first);
#'   2 = tips-first (terminal edges before internal); 3 = bucket ordering;
#'   4 = anti-tip (internal before terminal); 5 = large-first (most descendant
#'   taxa first).  On datasets with \eqn{\ge}65 tips, \code{clipOrder = 2L}
#'   (tips-first) typically increases replicate throughput by 5--15\% by
#'   evaluating higher-probability improvement candidates earlier.
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
#' @param stallEscalateFactor Numeric (>= 1); cross-replicate stall escalation.
#'   When a driven search stalls -- no improvement for `ceiling(nTip / 10)`
#'   consecutive replicates -- the ratchet perturbation probability is
#'   multiplied by this factor for each further `ceiling(nTip / 10)` replicates
#'   without improvement (capped at 0.5), and adaptive perturbation
#'   (`ratchetAdaptive`) is engaged, until an improvement resets the strength to
#'   its base value.  This lets a search discover at runtime the perturbation
#'   strength a difficult dataset needs, rather than relying on a fixed value.
#'   The default `1` disables escalation, leaving search behaviour unchanged.
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
#' @param rasStarts Integer; random-addition sequence (\acronym{RAS}) + TBR
#'   restarts per sector in XSS/RSS.  `1` (default) polishes the current sector
#'   with a single TBR pass;
#'   `n > 1` rebuilds the sector from scratch `n` times and keeps the best,
#'   following \insertCite{Goloboff1999;textual}{TreeSearch} RSS (TNT uses 3).
#'   Lets the search escape sector-local optima that a single TBR cannot leave.
#' @param sectorAcceptEqual Logical; accept equal-score sector resolutions in
#'   XSS/RSS (default `FALSE`).  On flat (e.g. missing-data) landscapes this lets
#'   the search traverse equally-parsimonious plateaus laterally rather than
#'   reverting every non-improving sector move, following Goloboff (2014).
#' @param sectorMaxHits Integer; equal-length trees the internal sector TBR holds
#'   while swapping a sector (default `1`).  TNT holds many; higher values let the
#'   sector search traverse equally-parsimonious plateaus (pairs with
#'   `sectorAcceptEqual`).
#' @param sectorCollapseTarget Integer; when `> 0`, a selected sector clade larger
#'   than this is **collapsed** into approximately this many composite terminals
#'   (deep sub-clades replaced by their first-pass state sets), so the sector
#'   search rearranges major sub-clades as a coarse skeleton rather than shuffling
#'   tips within a contiguous clade -- the reduced-dataset construction of
#'   \insertCite{Goloboff1999;textual}{TreeSearch}.  `0` (default) keeps the full
#'   fully-resolved clade.
#' @param rssPicks Integer; number of sector picks made **between** successive
#'   rounds of global \acronym{TBR} in random sectorial search.  `0` (default)
#'   auto-sizes to `2 * nTip / meanSectorSize` (~5).  Higher values chain more
#'   sequential sector replacements before the next global cleanup, matching the
#'   ~20-25 of \insertCite{Goloboff1999;textual}{TreeSearch} / \acronym{TNT}.
#' @param sectorGoDrift,sectorGoComb Integer; size thresholds (in real sector
#'   tips) above which a sector is solved by tree-drifting (`sectorGoDrift`) or by
#'   combined analysis (`sectorGoComb`) instead of plain
#'   \acronym{RAS}+\acronym{TBR}, following \acronym{TNT}'s `sectsch`
#'   `godrift`/`gocomb`.  Small sectors are cheap to optimise exhaustively; large
#'   sectors have more reach but need drift to escape their own local optima.
#'   `sectorGoComb` solves the sector with `sectorCombStarts` \acronym{RAS}+drift
#'   starts and then *fuses* them (recombines shared clades across the starts over
#'   `sectorFuseRounds` rounds), keeping the fused tree only if it beats the best
#'   start.  `sectorGoComb` takes precedence when both trigger.  `0` (default)
#'   disables each, leaving all sectors on \acronym{RAS}+\acronym{TBR}.  To
#'   exercise these, also raise `sectorMaxSize` so large sectors are selected.
#' @param sectorDriftCycles Integer; drift cycles used when a sector is solved by
#'   drift or combined analysis (\acronym{TNT} `drift`).
#' @param sectorDriftAfd Integer; absolute-fit-difference limit (steps) for
#'   in-sector drift.
#' @param sectorDriftRfd Numeric; relative-fit-difference limit for in-sector
#'   drift.
#' @param sectorCombStarts Integer; \acronym{RAS}+drift starts per combined
#'   sector (\acronym{TNT} `combstarts`).
#' @param sectorFuseRounds Integer; max fuse rounds for the combined-analysis
#'   recombination sub-step (\acronym{TNT} `fuse`).
#' @param postRatchetSectorial Logical; when `TRUE`, run XSS+RSS+CSS again
#'   after ratchet perturbation using the same round counts.  Approximates
#'   TNT's interleaved sectorial pattern.  Default: `FALSE`.
#' @param fuseInterval Integer; fuse pool trees every _n_ replicates.
#' @param fuseAcceptEqual Logical; accept equally-scoring fused trees?
#' @param intraFuse Logical; fuse the current tree against pool donors
#'   within each replicate, after TBR polish.  This approximates TNT's
#'   within-replicate fusing pattern. Default: `FALSE`.
#' @param poolMaxSize Integer; maximum trees retained in the pool.
#' @param poolSuboptimal Numeric; retain trees that are this many steps
#'   worse than the best tree.  0 (default) keeps only optimal trees.
#' @param consensusStableReps Integer; stop when the strict consensus of
#'   best-score pool trees has been unchanged for this many consecutive
#'   replicates.
#'   0 (default) disables this criterion; a typical value is 3--5.
#'   When both `consensusStableReps` and `targetHits` are active, the search
#'   stops when either criterion is met first.
#' @param perturbStopFactor Integer; stop when the number of consecutive
#'   replicates that fail to improve the best score exceeds
#'   `(targetHits / hits) * nTip * perturbStopFactor`, where `hits` is
#'   the number of replicates that have independently found the best score
#'   so far.  This scales patience inversely with progress toward
#'   `targetHits`: with few hits the threshold is large (more persistence);
#'   as hits approach `targetHits` the threshold converges to the flat
#'   `nTip * perturbStopFactor` limit.  Before any hit has been found
#'   (`hits == 0`) the criterion does not fire.
#'   When `targetHits` is disabled (0), falls back to the flat
#'   `nTip * perturbStopFactor` limit.
#'   0 disables this criterion entirely.
#'   Default 2.
#'   Inspired by IQ-TREE's unsuccessful-perturbation stopping rule
#'   \insertCite{Nguyen2015}{TreeSearch}; adapted from per-perturbation to
#'   per-replicate granularity.
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
#' @param pruneReinsertCycles Integer; number of taxon pruning-reinsertion
#'   perturbation cycles per replicate.  Each cycle drops a fraction of leaves,
#'   runs TBR on the reduced tree to let the backbone restructure, then
#'   greedily reinserts the dropped taxa via Wagner addition and TBR-polishes
#'   the full tree.  Complementary to the ratchet (which perturbs character
#'   weights) and NNI-perturbation (which perturbs the topology directly).
#'   0 (default) disables this perturbation.
#' @param pruneReinsertDrop Numeric (0--1); fraction of tips to drop per
#'   cycle.  Default 0.10 (10%).  Always drops at least 3 tips and keeps
#'   at least 4.
#' @param pruneReinsertSelection Integer; tip selection strategy for choosing
#'   which tips to drop:
#'   - `0` = random (default).
#'   - `1` = instability-weighted: tips whose parent-edge split is rare across
#'     pool trees are preferentially dropped.  Requires \eqn{\ge}2 pool trees;
#'     falls back to random otherwise.
#'   - `2` = missing-data-weighted: tips with more ambiguous or inapplicable
#'     characters are preferentially dropped.  High-missingness taxa are
#'     hardest to score correctly and most likely to be trapped in suboptimal
#'     positions.
#'   - `3` = combined: weight = instability × (1 + normalised missingness).
#'     Targets taxa that are both unstably placed and data-poor.
#' @param pruneReinsertTbrMoves Integer; maximum number of TBR moves accepted
#'   during the reduced-tree backbone optimisation phase of each
#'   prune-reinsert cycle.  0 means run to convergence; the default of 5
#'   mirrors the ratchet design (short perturbation, many diverse cycles)
#'   and substantially reduces per-cycle cost on datasets with inapplicable
#'   characters (where Brazeau scoring dominates).  Increase towards 0 if
#'   you prefer thorough backbone optimisation over replicate throughput.
#' @param pruneReinsertFullMoves Integer; maximum TBR moves during the
#'   full-tree polish after each prune-reinsert cycle.  0 (default) runs
#'   to convergence.  Has no effect when `pruneReinsertNni = TRUE`.
#' @param pruneReinsertNni Logical; if `TRUE`, use NNI (nearest-neighbour
#'   interchange) instead of TBR for the full-tree polish step.  NNI
#'   converges roughly 5x faster than TBR at large tip counts (\eqn{\ge}120),
#'   substantially reducing per-cycle cost while still reaching a local
#'   optimum before the outer-loop TBR polish.  Default `FALSE`.
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
    # TBR clip ordering strategy (experimental).
    # 0L=RANDOM (default), 1L=INV_WEIGHT (w=1/(1+s)), 2L=TIPS_FIRST,
    # 3L=BUCKET (tips/small/large), 4L=ANTI_TIP (non-tips first),
    # 5L=LARGE_FIRST (large then small then tips)
    clipOrder = 0L,
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
    # Default 12->6 (T-P5d, 2026-06-19): the ratchet was over-provisioned;
    # halving cycles saved 20-38% wall at zero quality loss on the mid-size EW
    # benchmarks.  The `thorough` preset (which `large` now rebases on) sets 20.
    ratchetCycles = 6L,
    ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = FALSE,
    ratchetTaper = FALSE,
    stallEscalateFactor = 1.0,
    # NNI perturbation
    nniPerturbCycles = 0L,
    nniPerturbFraction = 0.5,
    # Drift
    driftCycles = 0L,
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
    rasStarts = 1L,
    sectorAcceptEqual = FALSE,
    sectorMaxHits = 1L,
    sectorCollapseTarget = 0L,
    rssPicks = 0L,
    sectorGoDrift = 0L,
    sectorGoComb = 0L,
    sectorDriftCycles = 5L,
    sectorDriftAfd = 3L,
    sectorDriftRfd = 0.1,
    sectorCombStarts = 3L,
    sectorFuseRounds = 3L,
    postRatchetSectorial = FALSE,
    # Fuse / pool
    fuseInterval = 3L,
    fuseAcceptEqual = FALSE,
    intraFuse = FALSE,
    poolMaxSize = 100L,
    poolSuboptimal = 0,
    # Stopping criteria
    consensusStableReps = 0L,
    perturbStopFactor = 2L,
    adaptiveLevel = FALSE,
    consensusConstrain = FALSE,
    # Taxon pruning-reinsertion (T-266)
    pruneReinsertCycles = 0L,
    pruneReinsertDrop = 0.10,
    pruneReinsertSelection = 0L,
    pruneReinsertTbrMoves = 5L,
    pruneReinsertFullMoves = 0L,
    pruneReinsertNni = FALSE,
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
  # Record which fields the caller set explicitly (by name or position;
  # `match.call()` normalises positional args to their names).  This lets
  # `MaximizeParsimony()` distinguish a user-supplied `control` field from a
  # default, so a partial `control = SearchControl(...)` merges correctly with
  # a `strategy` preset instead of silently overriding every preset field.
  # A field set explicitly to its default value is still recorded as explicit
  # here, which a value-comparison against `SearchControl()` could not detect.
  .explicit <- names(as.list(match.call())[-1L])
  # Guard the count parameters whose non-positive values crash the C++ kernel:
  # `xssPartitions`/`cssPartitions` divide the tip count in `xss_partition()`
  # (integer division by zero -> SIGFPE), and `poolMaxSize` sizes the tree pool
  # whose eviction branch reads `entries_[0]` once `size >= max_size`
  # (an out-of-bounds read on an empty pool -> segfault). Each must be >= 1.
  for (.p in c("xssPartitions", "cssPartitions", "poolMaxSize")) {
    .v <- as.integer(get(.p))
    if (length(.v) != 1L || is.na(.v) || .v < 1L) {
      stop("`", .p, "` must be a single positive integer")
    }
  }
  # `stallEscalateFactor` multiplies the ratchet perturbation probability when a
  # run stalls; a value < 1 would *shrink* perturbation on stalling (the wrong
  # direction), and the C++ escalator treats exactly 1 as "off".
  .se <- as.double(stallEscalateFactor)
  if (length(.se) != 1L || is.na(.se) || .se < 1) {
    stop("`stallEscalateFactor` must be a single number >= 1")
  }
  structure(
    list(
      tbrMaxHits = as.integer(tbrMaxHits),
      clipOrder = as.integer(clipOrder),
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
      stallEscalateFactor = as.double(stallEscalateFactor),
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
      rasStarts = as.integer(rasStarts),
      sectorAcceptEqual = as.logical(sectorAcceptEqual),
      sectorMaxHits = as.integer(sectorMaxHits),
      sectorCollapseTarget = as.integer(sectorCollapseTarget),
      rssPicks = as.integer(rssPicks),
      sectorGoDrift = as.integer(sectorGoDrift),
      sectorGoComb = as.integer(sectorGoComb),
      sectorDriftCycles = as.integer(sectorDriftCycles),
      sectorDriftAfd = as.integer(sectorDriftAfd),
      sectorDriftRfd = as.double(sectorDriftRfd),
      sectorCombStarts = as.integer(sectorCombStarts),
      sectorFuseRounds = as.integer(sectorFuseRounds),
      postRatchetSectorial = as.logical(postRatchetSectorial),
      fuseInterval = as.integer(fuseInterval),
      fuseAcceptEqual = as.logical(fuseAcceptEqual),
      intraFuse = as.logical(intraFuse),
      poolMaxSize = as.integer(poolMaxSize),
      poolSuboptimal = as.double(poolSuboptimal),
      consensusStableReps = as.integer(consensusStableReps),
      perturbStopFactor = as.integer(perturbStopFactor),
      adaptiveLevel = as.logical(adaptiveLevel),
      consensusConstrain = as.logical(consensusConstrain),
      pruneReinsertCycles = as.integer(pruneReinsertCycles),
      pruneReinsertDrop = as.double(pruneReinsertDrop),
      pruneReinsertSelection = as.integer(pruneReinsertSelection),
      pruneReinsertTbrMoves = as.integer(pruneReinsertTbrMoves),
      pruneReinsertFullMoves = as.integer(pruneReinsertFullMoves),
      pruneReinsertNni = as.logical(pruneReinsertNni),
      annealCycles = as.integer(annealCycles),
      annealPhases = as.integer(annealPhases),
      annealTStart = as.double(annealTStart),
      annealTEnd = as.double(annealTEnd),
      annealMovesPerPhase = as.integer(annealMovesPerPhase),
      adaptiveStart = as.logical(adaptiveStart),
      enumTimeFraction = as.double(enumTimeFraction)
    ),
    class = "SearchControl",
    explicit = .explicit
  )
}

#' @export
print.SearchControl <- function(x, ...) {
  groups <- list(
    "TBR" = c("tbrMaxHits", "clipOrder", "nniFirst", "sprFirst", "tabuSize",
              "wagnerStarts", "wagnerBias", "wagnerBiasTemp", "outerCycles",
              "maxOuterResets"),
    "Ratchet" = c("ratchetCycles", "ratchetPerturbProb", "ratchetPerturbMode",
                   "ratchetPerturbMaxMoves", "ratchetAdaptive",
                   "ratchetTaper", "stallEscalateFactor"),
    "NNI Perturbation" = c("nniPerturbCycles", "nniPerturbFraction"),
    "Drift" = c("driftCycles", "driftAfdLimit", "driftRfdLimit"),
    "Prune-Reinsert" = c("pruneReinsertCycles", "pruneReinsertDrop",
                          "pruneReinsertSelection", "pruneReinsertTbrMoves",
                          "pruneReinsertFullMoves", "pruneReinsertNni"),
    "Annealing" = c("annealCycles", "annealPhases", "annealTStart",
                     "annealTEnd", "annealMovesPerPhase"),
    "Sectorial" = c("xssRounds", "xssPartitions", "rssRounds",
                     "cssRounds", "cssPartitions",
                     "sectorMinSize", "sectorMaxSize", "rasStarts",
                     "sectorAcceptEqual", "sectorMaxHits", "sectorCollapseTarget",
                     "rssPicks", "sectorGoDrift", "sectorGoComb",
                     "sectorDriftCycles", "sectorDriftAfd", "sectorDriftRfd",
                     "sectorCombStarts", "sectorFuseRounds",
                     "postRatchetSectorial"),
    "Fuse/Pool" = c("fuseInterval", "fuseAcceptEqual", "intraFuse",
                     "poolMaxSize", "poolSuboptimal"),
    "Stopping" = c("consensusStableReps", "perturbStopFactor",
                    "adaptiveLevel",
                    "consensusConstrain", "adaptiveStart",
                    "enumTimeFraction")
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
