# Internal helper: count non-missing taxa per character pattern.
# Used by XPIWE (Goloboff 2014) to compute the extrapolation factor.
# @param dataset A phyDat object.
# @return Integer vector of length = number of unique patterns.
# @keywords internal
.ObsCount <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  levels <- at$levels
  # "?" = all-1s contrast row.
  is_missing <- apply(contrast, 1, function(row) all(row == 1))
  # "-" (inapplicable/gap) also counts as missing for XPIWE (Goloboff 2014).
  # TNT counts both ? and - as missing, verified against TNT 1.6.
  inapp_col <- match("-", levels)
  if (!is.na(inapp_col)) {
    is_inapp <- apply(contrast, 1, function(row) {
      row[inapp_col] == 1 && sum(row) == 1
    })
    is_missing <- is_missing | is_inapp
  }
  # dataset is a list of integer vectors (token indices, 1-based) per taxon.
  # tip_data: n_taxa x n_patterns matrix
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  # Count non-missing taxa per pattern
  vapply(seq_len(ncol(tip_data)), function(p) {
    sum(!is_missing[tip_data[, p]])
  }, integer(1))
}

# Internal helper: recode inapplicable ("-") tokens as missing data ("?").
# Backs `inapplicable = "missing"` (pure-Fitch mode).  Every token whose
# contrast includes the gap state is promoted to the fully ambiguous "?"
# token, so the C++ simplification phase sees no genuine inapplicable token
# (`has_genuine_inapp` stays FALSE) and the character is scored with standard
# Fitch parsimony.  Working on the contrast matrix -- rather than
# round-tripping the data through a character matrix -- keeps the pattern
# structure and weights intact, and recodes {state, -} ambiguity tokens
# correctly: "0 or gap" = "0 or anything" = "?" (the round-trip left these as
# genuine inapplicable, which the engine then strips to a pure gap).
# @param dataset A phyDat object.
# @return The phyDat with every gap-bearing token recoded as missing.  If the
#   dataset has no "-" state it is returned unchanged.
# @keywords internal
.GapsAsMissing <- function(dataset) {
  gapCol <- match("-", attr(dataset, "levels"))
  if (is.na(gapCol)) {
    return(dataset)
  }
  contrast <- attr(dataset, "contrast")
  contrast[contrast[, gapCol] == 1, ] <- 1
  attr(dataset, "contrast") <- contrast
  dataset
}

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
    tabuSize = 0L, wagnerStarts = 1L,
    nniFirst = TRUE, sprFirst = FALSE
  ),
  default = SearchControl(
    # ratchetCycles 12->6 (T-P5d, 2026-06-19): profiling found the ratchet
    # over-provisioned -- halving cycles saved 20-38% wall on the mid-size EW
    # benchmarks (Wills/Zanol/Zhu/Giles) at zero quality loss.  Provisional;
    # the planned dataset-property grid will confirm across sizes.
    tbrMaxHits = 1L, ratchetCycles = 6L, ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 0L, ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = FALSE,
    driftCycles = 0L,
    xssRounds = 3L, xssPartitions = 4L,
    rssRounds = 1L, cssRounds = 0L, cssPartitions = 4L,
    sectorMinSize = 6L, sectorMaxSize = 50L,
    fuseInterval = 3L, fuseAcceptEqual = FALSE,
    tabuSize = 100L, wagnerStarts = 3L,
    nniFirst = TRUE, sprFirst = FALSE, adaptiveLevel = TRUE,
    maxOuterResets = 2L
  ),
  thorough = SearchControl(
    tbrMaxHits = 3L, ratchetCycles = 20L, ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = TRUE,
    nniPerturbCycles = 0L,  # T-274: 69% overhead, zero time-adjusted benefit
    driftCycles = 0L,
    xssRounds = 5L, xssPartitions = 6L,
    rssRounds = 3L, cssRounds = 2L, cssPartitions = 6L,
    sectorMinSize = 6L, sectorMaxSize = 80L,
    fuseInterval = 2L, fuseAcceptEqual = TRUE,
    tabuSize = 200L, wagnerStarts = 3L,
    nniFirst = TRUE, sprFirst = FALSE,
    outerCycles = 2L,
    maxOuterResets = 3L,
    adaptiveStart = TRUE
  ),
  # Opt-in "intensive" preset: `thorough` plus extra Wagner starts for more
  # starting-basin diversity.  Never auto-selected (.AutoStrategy returns only
  # sprint/default/thorough/large); the user opts in with strategy = "intensive".
  # Phase-2 sweep (2026-06-16, 5 seeds, EW Fitch): wagnerStarts 3->5 improved the
  # hardest datasets (Wortley2006 -3, Zhu2013 -2 toward the TNT optimum) at
  # neutral-to-lower candidate cost, with a ~+1-step trade-off on a couple of
  # others (Zanol2014, Giles2015) -- hence opt-in rather than a default change.
  # NB rasStarts=3 (TNT-faithful per-sector restarts) was evaluated 2026-06-18:
  # it closes the rss-ONLY gap (+7/+8 -> +1, wins time-matched) but is REDUNDANT
  # in the full thorough pipeline (Zanol/Zhu reach the optimum at rasStarts=1,
  # 60s) -- so NOT adopted.  Revisit for larger datasets / shorter budgets where
  # the full search can't reach the optimum (diag_thorough_rasstarts_tm.R +
  # the Hamilton grid t29_thorough_rasstarts_hamilton.sh).
  intensive = SearchControl(
    tbrMaxHits = 3L, ratchetCycles = 20L, ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = TRUE,
    nniPerturbCycles = 0L,
    driftCycles = 0L,
    xssRounds = 5L, xssPartitions = 6L,
    rssRounds = 3L, cssRounds = 2L, cssPartitions = 6L,
    sectorMinSize = 6L, sectorMaxSize = 80L,
    fuseInterval = 2L, fuseAcceptEqual = TRUE,
    tabuSize = 200L, wagnerStarts = 5L,
    nniFirst = TRUE, sprFirst = FALSE,
    outerCycles = 2L,
    maxOuterResets = 3L,
    adaptiveStart = TRUE
  ),
  # Large-tree preset (>=120 tips): at 180 tips each TBR convergence takes
  # ~5-7s, so phase costs scale sharply. Key design decisions (T-179):
  # - Fewer perturbation cycles: ratchet 12, drift 4 (vs thorough 20/12)
  # - No NNI-perturbation: at ~5.5s/cycle, it dominates the budget; ratchet
  #   provides more diverse escapes per unit time at large-tree scale
  # - Annealing (1 cycle) replaces drift: linear cooling T=20→0 over 5
  #   phases uses stochastic TBR with Boltzmann acceptance — cheaper
  #   per-cycle than drift. 1 cycle (400ms) captures 40% hit rate at
  #   180 tips; 3 cycles (1370ms) showed no significant score gain (T-248)
  # - No outer-cycle interleaving: outerCycles=1 avoids re-running expensive
  #   XSS/RSS/CSS after ratchet (saves ~10s per repeated sectorial pass)
  # - Single biased-Wagner start: saves ~2.6s vs 3 random starts; biased
  #   addition (Goloboff 2014) gives near-optimal Wagner at 180 tips
  # - tbrMaxHits=1: faster TBR passes (fewer equal-score trees explored)
  # - No adaptiveStart: with ~1 replicate per 60s budget, the bandit has
  #   no learning opportunity; adaptiveStart empirically regresses here
  # - Larger sector sizes for proportional tree coverage
  # - Prune-reinsert with NNI polish (T-289f Stage 5, 2026-03-29): 5 cycles,
  #   NNI full-tree polish (pruneReinsertNni=TRUE). TBR polish (Stage 4) was
  #   catastrophic at 206t/60s (0 reps). NNI polish (Stage 5, 5 datasets
  #   131-206t, 10 seeds, 60s+120s) fixes the 0-rep failure and improves
  #   median scores at 131-180t (project3701 146t: -178 steps at 60s;
  #   project804 173t: -9 steps; mbank_X30754 180t: -4 steps at 60s/-7 at
  #   120s). syab07205 (206t) shows +17.5 steps at 60s but neutral at 120s
  #   — acceptable given the gains at smaller sizes in range. See G-006 for
  #   a known limitation (NNI polish ignores ConstraintData; irrelevant here
  #   since the large preset does not use topological constraints).
  # Validated on mbank_X30754 (180t, 418p), 5 seeds at 30/60/120s budgets:
  #   60s:  large median=1255 vs thorough 1259 (+4 steps better)
  #   120s: large median=1250 vs thorough 1250 (tied, 2 reps vs 0-1)
  #   30s:  large median=1276 vs thorough 1283 (+7 steps better)
  large = SearchControl(
    tbrMaxHits = 1L, ratchetCycles = 12L, ratchetPerturbProb = 0.25,
    ratchetPerturbMode = 2L, ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = TRUE,
    nniPerturbCycles = 0L,
    driftCycles = 0L,
    annealCycles = 1L, annealPhases = 5L, annealTStart = 20, annealTEnd = 0,
    xssRounds = 3L, xssPartitions = 6L,
    rssRounds = 2L, cssRounds = 1L, cssPartitions = 6L,
    sectorMinSize = 8L, sectorMaxSize = 100L,
    fuseInterval = 3L, fuseAcceptEqual = TRUE,
    tabuSize = 100L, wagnerStarts = 1L,
    wagnerBias = 1L, wagnerBiasTemp = 0.3,
    nniFirst = TRUE, sprFirst = FALSE,
    outerCycles = 1L,
    pruneReinsertCycles = 5L, pruneReinsertNni = TRUE,
    consensusStableReps = 0L
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
# Merge a strategy preset into a (possibly user-customised) `SearchControl`.
# Fields the user set explicitly are preserved; every other field takes the
# preset's value.  A field counts as explicit if it was either
#   (a) passed as a top-level `...` argument (its name is in `explicitDots`), or
#   (b) supplied inside `control = SearchControl(...)` (its name is recorded in
#       the control's "explicit" attribute by SearchControl()).
# Reading the attribute -- rather than `names(control)` -- is the fix for the
# bug where `SearchControl()` always returns every field, which made the merge
# treat every field as explicit and apply nothing from the preset.
# @param control A SearchControl object (post-`...`-merge).
# @param preset The strategy preset (itself a SearchControl object).
# @param explicitDots Character vector of control-field names passed via `...`.
# @return `control` with preset values applied to non-explicit fields.
.ApplyStrategyPreset <- function(control, preset, explicitDots = character(0)) {
  explicitControl <- attr(control, "explicit")
  if (is.null(explicitControl)) {
    explicitControl <- character(0)
  }
  explicit <- union(explicitDots, explicitControl)
  for (nm in names(preset)) {
    if (!(nm %in% explicit)) {
      control[[nm]] <- preset[[nm]]
    }
  }
  control
}

.AutoStrategy <- function(nTip, nChar) {
  if (nTip <= 30L) return("sprint")
  # Few characters -> flat landscape; thorough search is pointless
  if (nChar < 100L) return("default")
  # Large trees (>=120 tips): per-replicate cost is high; use scaled preset
  # with NNI warmup and biased Wagner (empirically validated on 180-tip data).
  if (nTip >= 120L) return("large")
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
#' @param extended_iw Logical: if `TRUE` (default) and `concavity` is finite,
#'   apply the missing-entries correction of
#'   \insertCite{Goloboff2014;textual}{TreeSearch}.
#'   Characters with missing data receive a reduced effective concavity
#'   _k_c_ = _k_ / _f_c_, making their weights drop off faster.
#'   This compensates for the artificially low homoplasy of poorly sampled
#'   characters.  Set `FALSE` for legacy Goloboff (1993) behaviour.
#'   Ignored when `concavity = Inf` (equal weights) or `"profile"`.
#' @param xpiwe_r Numeric in (0, 1]: proportion of observed homoplasy
#'   expected in unobserved (missing) entries.  Default 0.5 (following TNT).
#'   Only used when `extended_iw = TRUE`.
#' @param xpiwe_max_f Numeric >= 1: maximum extrapolation factor.
#'   Characters with very few observed entries are clamped so that the
#'   extrapolation factor does not exceed this value.  Default 5 (following
#'   TNT).  Only used when `extended_iw = TRUE`.
#' @param hierarchy A [`CharacterHierarchy`] object specifying which
#'   characters are controlling primaries and which are their dependent
#'   secondaries.  Required when `inapplicable` is `"hsj"` or `"xform"`;
#'   ignored when `inapplicable = "bgs"` (the default).
#'   See [`CharacterHierarchy()`] for how to construct one, and
#'   [`HierarchyFromNames()`] for automated construction from
#'   TNT-style character names.
#' @param inapplicable Character: method for handling inapplicable characters.
#'   Case-insensitive.
#'   See `vignette("inapplicable", package = "TreeSearch")` for details.
#'   \describe{
#'     \item{`"bgs"` (default)}{Three-pass algorithm of
#'       \insertCite{Brazeau2019;textual}{TreeSearch}, inferring applicability
#'       regions from the `"-"` token.  No hierarchy required.}
#'     \item{`"missing"`}{Pure Fitch parsimony
#'       \insertCite{Fitch1971}{TreeSearch}: the inapplicable (`"-"`) state is
#'       treated as missing data, so any token that includes a gap is recoded
#'       as fully ambiguous (`"?"`) and contributes no steps -- including
#'       polymorphisms such as `{0,-}`, which become `?`.  Reproduces standard
#'       Fitch analyses (e.g. PAUP*, or TNT with gaps read as missing) that do
#'       not use the Brazeau-Gardner-Smith inapplicable algorithm.  No
#'       hierarchy required.}
#'     \item{`"hsj"`}{Dissimilarity-metric scoring of
#'       \insertCite{Hopkins2021;textual}{TreeSearch}.  Requires a
#'       `hierarchy`; controlled by `hsj_alpha`.}
#'     \item{`"xform"`}{Step-matrix recoding approximating maximum homology
#'       via x-transformations
#'       \insertCite{Goloboff2021;textual}{TreeSearch}.  Requires a
#'       `hierarchy`.}
#'   }
#' @param hsj_alpha Numeric in \[0, 1\]: scaling parameter for secondary-
#'   character contributions under the HSJ method.  0 = secondaries ignored;
#'   1 (default) = secondaries contribute up to 1 per branch per hierarchy
#'   block.  Only used when `inapplicable = "hsj"`.
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
#'       and character count:
#'       `"sprint"` for <=30 taxa; `"large"` for >=120 taxa with >=100
#'       character patterns; `"thorough"` for 65-119 taxa with >=100
#'       character patterns; `"default"` otherwise.}
#'     \item{`"sprint"`}{Fast search: 3 ratchet cycles, no drift, minimal
#'       sectorial. Good for small datasets or quick surveys.}
#'     \item{`"default"`}{Balanced: 12 ratchet + sectorial + fusing.}
#'     \item{`"thorough"`}{Intensive: 20 ratchet cycles, adaptive
#'       perturbation, extra sectorial rounds, NNI perturbation, outer cycle
#'       loop. Best for datasets with 65-119 tips and 100+ character patterns.}
#'     \item{`"large"`}{Large-tree search (>=120 tips): reduced cycle
#'       counts scaled for expensive per-replicate cost, no NNI
#'       perturbation, single biased Wagner start (Goloboff 2014), larger
#'       sector sizes, 1-cycle simulated annealing instead of drift
#'       (linear cooling from T=20 to T=0 over 5 phases).  Empirically matches
#'       or exceeds `"thorough"` at 180 tips across all time budgets.}
#'     \item{`"intensive"`}{Opt-in (never auto-selected): `"thorough"` plus extra
#'       Wagner starts (5) for more starting-basin diversity.  Improves the
#'       hardest datasets by a few steps at neutral-to-lower candidate cost, with
#'       an occasional ~+1-step trade-off elsewhere; choose it explicitly when
#'       pushing for the shortest tree on a difficult matrix.}
#'     \item{`"none"`}{Use only the explicitly supplied parameter values.}
#'   }
#'   Presets stop on `targetHits` and the `perturbStopFactor` no-improvement
#'   rule; `consensusStableReps` (consensus-stability stopping) is off by default
#'   and is not enabled by any preset.
#'   Explicit `control` fields always override the preset; for example,
#'   `strategy = "sprint", control = SearchControl(ratchetCycles = 10L)` uses
#'   sprint defaults for everything except `ratchetCycles`.
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates (default: 96).
#'   The default is a multiple of 48 (= LCM(12, 16)) so that replicates
#'   divide evenly across common 12- or 16-core machines when running in
#'   parallel.
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
#'     \item{`1` (default)}{Serial execution -- identical to previous behaviour.}
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
#'     \item{`n_topologies`}{Number of distinct topologies in the pool at the
#'       best score.}
#'     \item{`last_improved_rep`}{1-based index of the replicate that last
#'       improved the best score (0 if not tracked, e.g. parallel search).}
#'     \item{`timed_out`}{Logical: `TRUE` if the search stopped because
#'       `maxSeconds` was exceeded.}
#'     \item{`consensus_stable`}{Logical: `TRUE` if the search stopped
#'       because the strict consensus was unchanged for
#'       `consensusStableReps` consecutive replicates.}
#'     \item{`perturb_stop`}{Logical: `TRUE` if the search stopped because
#'       `nTip * perturbStopFactor` consecutive replicates failed to improve
#'       the best score (see [`SearchControl()`]).}
#'     \item{`timings`}{Named numeric vector of cumulative wall-clock time
#'       (in milliseconds) spent in each search phase across all replicates:
#'       `wagner_ms`, `tbr_ms`, `xss_ms`, `rss_ms`, `css_ms`, `ratchet_ms`,
#'       `drift_ms`, `final_tbr_ms`, `fuse_ms`.}
#'     \item{`replicate_scores`}{Numeric vector of the best parsimony score
#'       found by each completed replicate.  Passed to [ScoreSpectrum()] for
#'       Chao1-style landscape coverage estimation.}
#'     \item{`candidates_evaluated`}{Number of TBR/SPR-class candidate
#'       rearrangements evaluated across the whole search — the analogue of
#'       TNT's "rearrangements examined", useful for comparing search
#'       efficiency (candidates per unit of score improvement).  Counted only
#'       for single-threaded searches (`0` when `nThreads > 1`); excludes
#'       NNI-warmup and simulated-annealing candidates.}
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
#' @importFrom TreeTools NTip RandomTree Renumber RenumberTips RootTree
#' @importFrom TreeTools MakeTreeBinary Preorder
#' @importFrom cli cli_alert_success cli_alert_info cli_alert_warning
#' @encoding UTF-8
#' @export
MaximizeParsimony <- function(
    dataset,
    tree,
    concavity = Inf,
    extended_iw = TRUE,
    xpiwe_r = 0.5,
    xpiwe_max_f = 5,
    hierarchy = NULL,
    inapplicable = "bgs",
    hsj_alpha = 1.0,
    constraint,
    strategy = "auto",
    maxReplicates = 96L,
    targetHits = NULL,
    maxSeconds = 0,
    nThreads = 1L,
    verbosity = 1L,
    progressCallback = NULL,
    control = SearchControl(),
    ...
) {

  # --- Input validation: check dataset first ---
  if (is.null(dataset)) {
    stop("`dataset` cannot be NULL.")
  }

  # --- Set targetHits default if not provided ---
  if (is.null(targetHits)) {
    targetHits <- max(10L, as.integer(NTip(dataset) / 5))
  }

  # --- Backward compatibility: intercept maxTime → maxSeconds ---
  dots <- list(...)
  if ("maxTime" %in% names(dots)) {
    if (missing(maxSeconds) || maxSeconds == 0) {
      maxSeconds <- as.double(dots[["maxTime"]])
    }
    .Deprecated(msg = paste0(
      "Use `maxSeconds` instead of `maxTime` in MaximizeParsimony().\n",
      "  `maxTime` was a Morphy()-style parameter; `maxSeconds` is the ",
      "equivalent for the new C++ search engine."
    ))
    dots[["maxTime"]] <- NULL
  }

  # --- Backward compatibility: detect Morphy()-style parameters ---
  .morphyParams <- c("ratchIter", "tbrIter", "startIter", "finalIter",
                      "maxHits", "quickHits", "ratchEW",
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
    if (!missing(tree) && !is.null(tree)) morphyArgs$tree <- tree
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
      control <- .ApplyStrategyPreset(control, preset, names(controlDots))
      if (verbosity >= 1L) {
        cli::cli_alert_info("Strategy: {.strong {strategy}}")
      }
    } else if (!identical(strategy, "auto")) {
      warning("Unknown strategy '", strategy, "'; using default parameters.")
    }
  }

  # --- Progress callback: build default cli bar if needed ---
  if (is.null(progressCallback) && verbosity >= 1L && interactive()) {
    pb_env <- new.env(parent = environment())
    pb_env$id <- cli::cli_progress_bar(
      total = as.integer(maxReplicates),
      format = paste0(
        "Rep {cli::pb_current}/{cli::pb_total}",
        " | Best: {best}",
        " | Hits: {hits}/{target}"
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

  # --- Progress file callback (for Shiny background futures) ---
  if (is.null(progressCallback)) {
    progressFile <- Sys.getenv("TREESEARCH_PROGRESS_FILE", "")
    if (nzchar(progressFile)) {
      progressCallback <- function(info) {
        if (identical(info$phase, "replicate")) {
          tryCatch(
            writeLines(paste(info$replicate, info$max_replicates,
                             signif(info$best_score, 8), info$hits_to_best,
                             info$target_hits),
                       progressFile),
            error = function(e) NULL
          )
        }
      }
    }
  }

  # --- Profile parsimony: prepare data ---
  useProfile <- !missing(concavity) && identical(concavity, "profile")
  if (useProfile) {
    profileApprox <- if (!is.null(dots[["profile_approx"]])) {
      dots[["profile_approx"]]
    } else {
      "auto"
    }
    dataset <- PrepareDataProfile(dataset, approx = profileApprox)
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

  # --- Validate inapplicable-handling parameters ---
  inapplicable <- tolower(inapplicable)
  if (inapplicable == "brazeau") inapplicable <- "bgs"
  inapplicable <- match.arg(inapplicable, c("bgs", "hsj", "xform", "missing"))
  # "missing" = pure Fitch: recode every gap-bearing token as missing ("?") so
  # gaps contribute no steps, then score with the standard engine (which on
  # inapplicable-free data reduces to Fitch parsimony).
  if (inapplicable == "missing") {
    dataset <- .GapsAsMissing(dataset)
    inapplicable <- "bgs"
  }
  if (inapplicable != "bgs") {
    if (is.null(hierarchy)) {
      stop("A `hierarchy` is required when inapplicable = \"", inapplicable,
           "\". See ?CharacterHierarchy.")
    }
    if (!inherits(hierarchy, "CharacterHierarchy")) {
      stop("`hierarchy` must be a CharacterHierarchy object.")
    }
    ValidateHierarchy(hierarchy, dataset)
    if (useProfile) {
      stop("Profile parsimony is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
    if (is.finite(concavity)) {
      stop("Implied weighting is not currently supported with inapplicable = \"",
           inapplicable, "\".")
    }
    # xform validation is done; recoding happens below
  }
  if (!is.numeric(hsj_alpha) || length(hsj_alpha) != 1L ||
      hsj_alpha < 0 || hsj_alpha > 1) {
    stop("`hsj_alpha` must be a single number in [0, 1].")
  }
  if (is.finite(concavity) && concavity <= 0) {
    stop("`concavity` must be positive (or Inf for equal weights, ",
         "or \"profile\" for profile parsimony).")
  }

  # --- Starting tree ---
  userTree <- !missing(tree) && !is.null(tree)
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
  weight <- .ScaleWeight(at$weight)
  levels <- at$levels

  # --- Replicate count adequacy check ---
  # Warn only when the user explicitly passed maxReplicates.
  # Formula: max(10, ceiling(nTip * nChar / 5000)) where nChar = sum(weight).
  # Derived from T-069 benchmarks: at 225 taxa / 748 chars a single rep takes
  # ~40s and at least ~34 reps are needed to fill the tree pool reliably.
  if (!missing(maxReplicates) && nTip >= 30L && verbosity > 0L) {
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

  # --- HSJ: prepare hierarchy data for C++ ---
  hsjArgs <- list()
  useHSJ <- !is.null(hierarchy) && identical(inapplicable, "hsj")
  if (useHSJ) {
    hsjArgs$hierarchyBlocks <- .HierarchyToBlocks(hierarchy)
    hsjArgs$hsjTipLabels <- .BuildTipLabels(dataset)
    hsjArgs$hsjAlpha <- as.double(hsj_alpha)
    # 0-based token index of the primary's "absent" state (depends on level
    # ordering, so computed from the data rather than hard-coded).
    hsjArgs$hsjAbsentState <- .HSJAbsentState(dataset)

    # Adjust weights: subtract hierarchy characters so Fitch scores non-hierarchy
    adj_weight <- .NonHierarchyWeights(dataset, hierarchy)
    weight <- as.integer(adj_weight)
  }

  # --- Xform: recode hierarchy into step-matrix characters ---
  xformArgs <- list()
  useXform <- !is.null(hierarchy) && identical(inapplicable, "xform")
  if (useXform) {
    recoded <- RecodeHierarchy(dataset, hierarchy)
    xformArgs$xformChars <- recoded$sankoff_chars

    # Adjust weights: subtract hierarchy characters so Fitch scores non-hierarchy
    adj_weight <- .NonHierarchyWeights(dataset, hierarchy)
    weight <- as.integer(adj_weight)
  }

  # --- IW: compute minimum step counts per character ---
  if (is.finite(concavity)) {
    minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
  }

  # --- XPIWE: compute per-pattern observed-taxa counts ---
  useXpiwe <- isTRUE(extended_iw) && is.finite(concavity) && !useProfile
  if (useXpiwe) {
    obsCount <- .ObsCount(dataset)
  }

  # --- Run C++ driven search ---
  # searchControl: the resolved SearchControl object (already type-coerced)
  # runtimeConfig: session-level params not in SearchControl
  runtimeConfig <- list(
    maxReplicates = as.integer(maxReplicates),
    targetHits = as.integer(targetHits),
    maxSeconds = as.double(maxSeconds),
    verbosity = as.integer(verbosity),
    nThreads = as.integer(nThreads),
    startEdge = if (userTree) tree[["edge"]] else NULL,
    progressCallback = progressCallback
  )

  # scoringConfig: scoring method params
  scoringConfig <- list(
    min_steps = if (is.finite(concavity)) minSteps else integer(0),
    concavity = as.double(concavity),
    xpiwe = useXpiwe,
    xpiwe_r = as.double(xpiwe_r),
    xpiwe_max_f = as.double(xpiwe_max_f),
    obs_count = if (useXpiwe) obsCount else integer(0),
    infoAmounts = profileArgs$infoAmounts
  )

  # constraintConfig / hsjConfig / xformConfig: NULL when empty
  constraintConfig <- if (length(consArgs) > 0L) consArgs
  hsjConfig <- if (length(hsjArgs) > 0L) hsjArgs
  xformConfig <- if (length(xformArgs) > 0L) xformArgs

  result <- ts_driven_search(
    contrast, tip_data, weight, levels,
    control, runtimeConfig, scoringConfig,
    constraintConfig, hsjConfig, xformConfig
  )

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
    total_s <- round(sum(unlist(result$timings), na.rm = TRUE) / 1000, 1)
    stop_reason <- if (isTRUE(result$timed_out)) "timeout"
                   else if (isTRUE(result$consensus_stable)) "consensus stable"
                   else if (isTRUE(result$perturb_stop)) "perturbation limit"
                   else "replicate limit"
    cli_alert_success(paste0(
      "Search complete: score {.strong {signif(result$best_score, 7)}}, ",
      "{result$replicates} replicate{?s} ",
      "(last improved: #{result$last_improved_rep}), ",
      "{result$hits_to_best} hit{?s} to best, ",
      "{result$n_topologies} MPT{?s}, ",
      "stop: {stop_reason}, {total_s}s"
    ))
  }

  structure(
    outTrees,
    score = result$best_score,
    replicates = result$replicates,
    hits_to_best = result$hits_to_best,
    n_topologies = result$n_topologies,
    last_improved_rep = result$last_improved_rep,
    timed_out = isTRUE(result$timed_out),
    consensus_stable = isTRUE(result$consensus_stable),
    perturb_stop = isTRUE(result$perturb_stop),
    timings = unlist(result$timings),
    strategy_diagnostics = result$strategy_diagnostics,
    replicate_scores = result$replicate_scores,
    candidates_evaluated = result$candidates_evaluated,
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
