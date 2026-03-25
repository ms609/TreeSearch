# Backward-compatible wrapper for ts_driven_search.
#
# Accepts the old flat-argument calling convention used by tests and
# packs them into the grouped lists expected by ts_driven_search().
# Production code (MaximizeParsimony, .ResampleHierarchy) calls
# ts_driven_search() directly with pre-built grouped lists.
ts_driven_search <- function(
    contrast,
    tip_data,
    weight,
    levels,
    # --- New grouped-list interface (used when calling with grouped args) ---
    searchControl = NULL,
    runtimeConfig = NULL,
    scoringConfig = NULL,
    constraintConfig = NULL,
    hsjConfig = NULL,
    xformConfig = NULL,
    # --- Old flat-argument interface (used by tests) ---
    maxReplicates = 100L,
    targetHits = 10L,
    tbrMaxHits = 1L,
    ratchetCycles = 10L,
    ratchetPerturbProb = 0.04,
    ratchetPerturbMode = 0L,
    ratchetPerturbMaxMoves = 0L,
    ratchetAdaptive = FALSE,
    ratchetTaper = FALSE,
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
    maxSeconds = 0.0,
    verbosity = 0L,
    min_steps = integer(0),
    concavity = -1.0,
    consSplitMatrix = NULL,
    consContrast = NULL,
    consTipData = NULL,
    consWeight = NULL,
    consLevels = NULL,
    consExpectedScore = 0L,
    infoAmounts = NULL,
    tabuSize = 100L,
    wagnerStarts = 1L,
    progressCallback = NULL,
    nThreads = 1L,
    startEdge = NULL,
    sprFirst = FALSE,
    nniFirst = TRUE,
    hierarchyBlocks = NULL,
    hsjTipLabels = NULL,
    hsjAlpha = 1.0,
    hsjAbsentState = 0L,
    xformChars = NULL,
    xpiwe = FALSE,
    xpiwe_r = 0.5,
    xpiwe_max_f = 5.0,
    obs_count = integer(0),
    consensusStableReps = 0L,
    perturbStopFactor = 0L,
    adaptiveLevel = FALSE,
    consensusConstrain = FALSE,
    nniPerturbCycles = 0L,
    nniPerturbFraction = 0.5,
    wagnerBias = 0L,
    wagnerBiasTemp = 0.3,
    outerCycles = 1L,
    maxOuterResets = 0L,
    adaptiveStart = FALSE,
    enumTimeFraction = 0.1,
    annealConfig = NULL)
{
  # New-style call: grouped lists already provided
  if (!is.null(searchControl)) {
    return(.Call(`_TreeSearch_ts_driven_search`,
      contrast, tip_data, weight, levels,
      searchControl, runtimeConfig, scoringConfig,
      constraintConfig, hsjConfig, xformConfig
    ))
  }

  # Old-style call: pack flat args into grouped lists
  sc <- SearchControl(
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
    adaptiveStart = as.logical(adaptiveStart),
    enumTimeFraction = as.double(enumTimeFraction)
  )

  # Anneal config: fold into SearchControl if provided
  if (!is.null(annealConfig)) {
    sc$annealPhases <- as.integer(annealConfig$phases %||% 0L)
    sc$annealTStart <- as.double(annealConfig$tStart %||% 20)
    sc$annealTEnd <- as.double(annealConfig$tEnd %||% 0)
    sc$annealMovesPerPhase <- as.integer(annealConfig$movesPerPhase %||% 0L)
  }

  rt <- list(
    maxReplicates = as.integer(maxReplicates),
    targetHits = as.integer(targetHits),
    maxSeconds = as.double(maxSeconds),
    verbosity = as.integer(verbosity),
    nThreads = as.integer(nThreads),
    startEdge = startEdge,
    progressCallback = progressCallback
  )

  scoring <- list(
    min_steps = min_steps,
    concavity = as.double(concavity),
    xpiwe = as.logical(xpiwe),
    xpiwe_r = as.double(xpiwe_r),
    xpiwe_max_f = as.double(xpiwe_max_f),
    obs_count = obs_count,
    infoAmounts = infoAmounts
  )

  # Constraint config
  cc <- NULL
  if (!is.null(consSplitMatrix)) {
    cc <- list(
      consSplitMatrix = consSplitMatrix,
      consContrast = consContrast,
      consTipData = consTipData,
      consWeight = consWeight,
      consLevels = consLevels,
      consExpectedScore = as.integer(consExpectedScore)
    )
  }

  # HSJ config
  hc <- NULL
  if (!is.null(hierarchyBlocks)) {
    hc <- list(
      hierarchyBlocks = hierarchyBlocks,
      hsjTipLabels = hsjTipLabels,
      hsjAlpha = as.double(hsjAlpha),
      hsjAbsentState = as.integer(hsjAbsentState)
    )
  }

  # Xform config
  xc <- NULL
  if (!is.null(xformChars)) {
    xc <- list(xformChars = xformChars)
  }

  .Call(`_TreeSearch_ts_driven_search`,
    contrast, tip_data, weight, levels,
    sc, rt, scoring, cc, hc, xc
  )
}
