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
  # Drop the IW minimum-length cache: it is keyed on the (now altered) contrast,
  # and TreeLength() reuses it rather than recomputing when present.
  attr(dataset, "min.length") <- NULL
  dataset
}

# Internal helper: structural sanity check on a user-supplied starting tree.
#
# Deliberately not `ape::checkValidPhylo()`, which prints a report rather than
# signalling a condition.  This checks only the invariants whose violation
# makes TreeTools' C++ rooting and traversal routines index out of bounds --
# a segfault the caller cannot trap, so it has to be pre-empted rather than
# handled.  Reachable in practice: `ape::unroot()` accepts TreeTools' `order =
# "preorder"` attribute and then mishandles it, so unrooting any TreeTools
# tree returns an edge matrix containing NA.
# @param tr A candidate starting tree.
# @param i Index within the supplied pool, or `NA_integer_` for a lone tree.
# @return `tr`, invisibly; called for the error.
# @keywords internal
.CheckStartTree <- function(tr, i) {
  what <- if (is.na(i)) "`tree`" else paste0("`tree[[", i, "]]`")
  edge <- tr[["edge"]]
  if (!is.matrix(edge) || dim(edge)[2L] != 2L || !is.numeric(edge) ||
      anyNA(edge)) {
    stop(what, " has a malformed edge matrix.")
  }
  nTip <- length(tr[["tip.label"]])
  child <- edge[, 2L]
  if (!identical(sort(as.integer(child[child <= nTip])), seq_len(nTip))) {
    stop(what, " is not a valid tree: every leaf must be the child of ",
         "exactly one edge.")
  }
  if (any(edge[, 1L] <= nTip)) {
    stop(what, " is not a valid tree: a leaf cannot be a parent.")
  }
  invisible(tr)
}

# Internal helper: prepare constraint data for C++ engine.
# Returns a named list of constraint arguments (empty list if no constraint).
# @param constraint A phyDat, phylo, or NULL.
# @param dataset A phyDat whose names define the tip ordering.
# @keywords internal
#' @importFrom TreeTools AddUnconstrained
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
  # For each constraint character, record the tips unambiguously in the "1"
  # group (derived state present, ancestral absent) and, separately, those in
  # the "0" group (ancestral present, derived absent).  Tips ambiguous for the
  # character ("?", or unconstrained taxa) belong to neither group and are free
  # to plot on either side of the split.
  consSplits <- matrix(0L, nrow = ncol(consMat), ncol = length(constraint))
  consZero   <- matrix(0L, nrow = ncol(consMat), ncol = length(constraint))
  for (ch in seq_len(ncol(consMat))) {
    for (tip in seq_len(length(constraint))) {
      token <- consMat[tip, ch]
      if (consContrast[token, nConsStates] == 1 &&
          consContrast[token, 1] == 0) {
        consSplits[ch, tip] <- 1L
      } else if (consContrast[token, 1] == 1 &&
                 consContrast[token, nConsStates] == 0) {
        consZero[ch, tip] <- 1L
      }
    }
  }

  keep <- apply(consSplits, 1, function(row) {
    s <- sum(row)
    s >= 1 && s < length(constraint) - 1
  })
  consSplits <- consSplits[keep, , drop = FALSE]
  consZero   <- consZero[keep, , drop = FALSE]
  if (nrow(consSplits) == 0L) return(list())

  # Every returned tree must display all constraint splits simultaneously.
  # Two splits are jointly displayable iff they are compatible in the
  # four-gamete sense: treating each as a bipartition of the tips it
  # constrains (its "1" group vs its "0" group, ambiguous tips excluded), the
  # pair is compatible iff at least one of the four group intersections is
  # empty.  A laminar (nested-or-disjoint) test alone is too strict: it rejects
  # the case where the two "0" groups are disjoint -- i.e. the splits' "1"
  # sides jointly cover the constrained tips -- which is perfectly displayable,
  # e.g. ab | cef and abcd | ef coexist on ((a,b),(d,(c,(e,f)))).
  nSplits <- nrow(consSplits)
  if (nSplits > 1L) {
    for (i in seq_len(nSplits - 1L)) {
      aOne  <- consSplits[i, ] == 1L
      aZero <- consZero[i, ] == 1L
      for (j in seq(i + 1L, nSplits)) {
        bOne  <- consSplits[j, ] == 1L
        bZero <- consZero[j, ] == 1L
        compatible <- !any(aOne & bOne) || !any(aOne & bZero) ||
                      !any(aZero & bOne) || !any(aZero & bZero)
        if (!compatible) {
          stop("Constraint is impossible to satisfy: splits ", i, " and ", j,
               " are incompatible (all four taxon groupings co-occur)")
        }
      }
    }
  }

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

# Ratchet depth for implied weights under `thorough`/`large`, applied after the
# strategy preset (see MaximizeParsimony()). Kept out of `.StrategyPresets()` so
# the preset table stays scorer-agnostic: this depth is calibrated for implied
# weights only, and equal weights measurably does not want it.
.iwRatchetCycles <- 48L
# Largest depth with supporting measurements; user escalation is capped here
# rather than extrapolated. Quoted as a literal in the `targetHits` docs -- keep
# the two in step if this changes.
.iwRatchetMaxCycles <- 115L

# Ratchet depth to impose for this call, or NULL to leave the preset's value.
# `userSet` names the fields the caller set themselves (never overridden).
# See the call site in MaximizeParsimony() for the calibration behind it.
.IwRatchetDepth <- function(strategy, concavity, targetHits, defaultHits,
                            userSet = character(0)) {
  if (!length(strategy) || !strategy %in% c("thorough", "large")) {
    return(NULL)
  }
  # `concavity` may still be the "profile" sentinel here: profile parsimony is a
  # different objective and is left alone, as is equal weights (infinite).
  if (length(concavity) != 1L || !is.numeric(concavity) ||
      !is.finite(concavity)) {
    return(NULL)
  }
  if ("ratchetCycles" %in% userSet) {
    return(NULL)
  }
  escalation <- if (length(defaultHits) == 1L && is.finite(defaultHits) &&
                    defaultHits > 0 && length(targetHits) == 1L &&
                    is.finite(targetHits)) {
    max(1, targetHits / defaultHits)
  } else {
    1
  }
  min(.iwRatchetMaxCycles, as.integer(round(.iwRatchetCycles * escalation)))
}

# Implied-weights operating point for `sprint`/`default`: a deeper ratchet paid
# for by a flat replicate patience.  Same scoping rules as .IwRatchetDepth()
# above (implied weights only, never override the caller), and deliberately
# disjoint from it by strategy so the two can never both fire.
# See the call site in MaximizeParsimony() for the measurements.
.iwStopPackage <- list(
  sprint  = list(ratchetCycles = 12L, ratchetPerturbProb = 0.25,
                 stopPatience = 20L),
  # `ratchetPerturbProb` is already 0.25 in the preset, so it is absent here:
  # this list names only what the implied-weights measurement actually moved.
  default = list(ratchetCycles = 20L, stopPatience = 15L)
)

# Named list of control fields to impose for this call, or NULL for none.
# `userSet` names the fields the caller set themselves; those are dropped from
# the returned list rather than filtered at the call site, keeping the
# never-override-the-user rule in one place.
.IwStopPackage <- function(strategy, concavity, userSet = character(0)) {
  if (!length(strategy) || !strategy %in% names(.iwStopPackage)) {
    return(NULL)
  }
  # As in .IwRatchetDepth(): `concavity` may still be the "profile" sentinel, and
  # equal weights is infinite.  Both are different objectives, and neither was
  # measured here.
  if (length(concavity) != 1L || !is.numeric(concavity) ||
      !is.finite(concavity)) {
    return(NULL)
  }
  out <- .iwStopPackage[[strategy]]
  out <- out[setdiff(names(out), userSet)]
  if (!length(out)) NULL else out
}

# Strategy presets for adaptive search (Phase 6E).
# Wrapped in a function to avoid load-order dependency on SearchControl().
.StrategyPresets <- function() {
  presets <- list(
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
    # driftCycles 0->2 + wagnerStarts 3->5 (two-island sweep 2026-06-25, 30 seeds):
    # drift recovers equal-score trees on TBR-disconnected islands (uphill tunnelling
    # across the barrier; Zhu2013 two-island recovery 0.73 -> 0.95; ws5 alone hurts
    # it, 0.70).  COST (anytime study, 20 training matrices 65-120 tips, 2026-07-02):
    # drift is per-replicate overhead -> at a fixed budget thorough completes fewer
    # reps and reaches the optimum somewhat LESS reliably on the general pool; it is
    # a deliberately higher-effort/slower tier that needs a larger replicate budget
    # to converge. NB thorough is auto-selected for 65-119 tips, so this cost lands on the
    # default path there.
    driftCycles = 2L,
    xssRounds = 5L, xssPartitions = 6L,
    rssRounds = 3L, cssRounds = 2L, cssPartitions = 6L,
    sectorMinSize = 6L, sectorMaxSize = 80L,
    # In-sector drifting for large sectors (TNT `godrift`).  Preset-level matched-
    # wall A/Bs (2026-07-08, arrays 17836031 then 17836637; sector-resolve + general-
    # pool 68-88t + large training 131-205t x 5 seeds): rasStarts = 3 + sectorGoDrift
    # = 25 + sectorDriftCycles = 3 beats stock `thorough` on every class at matched
    # wall with no regression (mid-size sector -1.5, general -0.4; large -2.7, with
    # the rep-starved 205t project3763 -8, hard-floor project4138 reaching optimum).
    # sectorGoDrift = 25 is calibrated to this preset's sector geometry: xss/css
    # sectors are ~12-15 tips (xss/cssPartitions = 6), so drift engages via the RSS
    # large-clade picks; 40 is near-inert here.  rasStarts >= 2 is REQUIRED for the
    # drift retention channel and is coupled -- rasStarts = 3 ALONE (no drift)
    # regresses (triples every sector-solve for fewer reps), but the drift redeems
    # the cost; rasStarts = 2 + drift ties 3 + drift, so 3 (marginally best) is kept.
    rasStarts = 3L,
    sectorGoDrift = 25L, sectorDriftCycles = 3L,
    fuseInterval = 2L, fuseAcceptEqual = TRUE,
    tabuSize = 200L, wagnerStarts = 5L,
    nniFirst = TRUE, sprFirst = FALSE,
    outerCycles = 2L,
    maxOuterResets = 3L,
    adaptiveStart = TRUE
  )
  )

  # `intensive` is retained as a backward-compatible alias of `thorough`.  The
  # 2026-06-25 two-island sweep (30 seeds) folded wagnerStarts = 5 (intensive's
  # sole distinguishing feature) into `thorough` together with driftCycles = 2;
  # ws5 showed no score gain over thorough while costing wall-clock, so the two
  # presets are merged.  Retained as an internal alias only: the effort ladder
  # never names it, and there is no user-facing way to ask for it.
  presets$intensive <- presets$thorough

  # Large-tree preset (>=120 tips).  REBASED on `thorough` (2026-07-07).
  # The former bespoke `large` (T-179) was a cost-cut that predated the
  # thorough/auto overhaul: outerCycles=1 (never re-ran sectorial after
  # ratchet), wagnerStarts=1, driftCycles=0, adaptiveStart=FALSE, xss/rss/css
  # 3/2/1.  Two sweeps (fixed engine, MPT-reach metric, project175 held out):
  #   - Long matched wall (1200s, 6 matrices 125-482t): thorough reach 0.61 vs
  #     large 0.28; `outerCycles=2` is the load-bearing knob (ablation).
  #   - Short-budget gate (30/60/120s, 5 matrices 125-199t x 5 seeds): thorough
  #     dominates on reach at >=60s and ties `large+oc2` at 30s; the apparent
  #     30s dip was one seed on one already-solved matrix (noise, tied gap).
  # thorough ran unstarved at 482t/1200s (rate run 17819704, ~187 reps/seed),
  # confirming its heavier provisioning does not rep-starve at large scale now
  # that maxReplicates=500 is the `large` default (raised this session).
  # DROPPED vs old large (all superseded by thorough's machinery in-sweep):
  #   annealCycles (drift replaces it), biased-Wagner start, prune-reinsert NNI
  #   polish (T-289f), tbrMaxHits=1.  ADOPTED: ratchet 20, drift 2, xss/rss/css
  #   5/3/2, sectorMax 80, wagnerStarts 5, outerCycles 2, adaptiveStart TRUE.
  # RESIDUAL untested corner: >199t at a *tight* (<60s) user budget — the
  # min-replicate / budget-aware fallback (auto-routing arm) is the follow-up,
  # not a reason to withhold the rebase.  The per-replicate reach deficit that
  # survives (project5432/4138: 0 hits over 1870 reps) is an ENGINE limit
  # (cross-set sectorial re-solve), which no preset provisioning can close.
  presets$large <- presets$thorough

  presets
}

# Calibration behind .AutoRung()'s size/character thresholds.
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

# --- The effort ladder -----------------------------------------------------
#
# Rungs 1-3 are the provisioning presets.  Rung 4 is `thorough`'s provisioning
# with a raised replicate cap -- which is exactly what `large` already was
# (`presets$large <- presets$thorough`, plus `maxReplicates = 500`).  So the
# ladder generalises an axis the package was already using; it does not invent
# one.  Above rung 4 only the BUDGET climbs, because provisioning saturates at
# `thorough`: there is nothing further to provision.
#
# These names are internal labels for menu entries, not a user-facing argument.
# Users ask for effort relative to the automatic choice; only the package (and
# its tests) name a rung, via the internal `.rung` argument.
.effortLadder <- c("sprint", "default", "thorough", "large")

# A REPRESENTABILITY limit, not a policy one -- and the distinction matters.
# `.iwRatchetMaxCycles = 115` caps at the largest ratchet depth actually tested,
# because extra ratchet depth is not known to be free.  Extra replicates ARE:
# raising the cap only appends later replicates and never delays an earlier
# improvement, so reach is monotone non-decreasing in `maxReplicates` and the
# only cost is wall -- which is exactly what someone raising `effort` is asking
# to spend.  There is therefore no measured or principled level at which the
# ladder should refuse to go further, and an arbitrary ceiling would just
# obstruct the request.
#
# Rung 26 is where `500 * 2^(rung - 4)` stops fitting in R's integer type
# (500 * 2^22 = 2 097 152 000; one more doubling overflows).  Requests beyond it
# are clamped WITH A MESSAGE, so `effort = 40` announces that it means the same
# as `effort = 26` rather than silently pretending otherwise.
.effortMaxRung <- 26L

# Everything a rung means, in ONE place.  `maxReplicates = NA` means "leave the
# SearchControl default alone".
.RungSpec <- function(rung) {
  rung <- as.integer(rung)
  list(
    preset = .effortLadder[[min(rung, length(.effortLadder))]],
    # 96 (the SearchControl default) through rung 3; 500 at rung 4 -- the value
    # `large` already used -- then doubling.  Raising this cap ANYTIME-DOMINATES
    # (a higher cap only appends later replicates; it never delays an earlier
    # improvement) and easy datasets still stop early on `targetHits`, so the
    # cost falls only on the genuinely hard tail that runs to the cap.  That is
    # what licenses extrapolating this knob past the measured 500 when the
    # ratchet depth may not be extrapolated.
    maxReplicates = if (rung <= 3L) NA_integer_ else
      as.integer(500 * 2^(rung - 4L)),
    # `targetHits` multiplier: 1 through rung 4, then doubling in step with the
    # replicate budget.
    #
    # BOTH knobs double, so that one notch means the same thing -- roughly twice
    # the work -- whichever population a dataset falls in.  They govern disjoint
    # populations (see below), so mixing rates would make a notch 2x the work on
    # hard matrices but only (k+1)/k on easy ones, i.e. notches would shrink as
    # you climb on precisely the population `targetHits` controls.  That is the
    # only argument for the shape; it is an OPERATING POINT, not a fitted
    # constant.  What is measured is that reach was still climbing at 500
    # replicates with no knee (34-matrix 120-180t sweep, reach@96 = 0.68 ->
    # reach@250 = 0.79) -- so more is better, and nothing measures where that
    # stops or what shape the approach has.  A doubling grid over rungs 4-8 on
    # the hard tail is what would replace this guess with a measurement.
    #
    # `maxReplicates` deliberately leads and `targetHits` follows, because the
    # two bite on DISJOINT populations.  `targetHits` ends a run early on easy
    # datasets, so raising it lengthens those; on hard datasets it is never
    # reached and `maxReplicates` binds first.  Measured (array 18096945,
    # equal weights): Zanol2014 ran the full 96 replicates at hits-to-best = 1
    # against a target of 14, and tripling `targetHits` to 42 changed score,
    # replicate count and wall not at all.  A rung that raised `targetHits`
    # alone would therefore do nothing on precisely the datasets someone turns
    # effort up for.
    #
    # It still earns its place from rung 5: it buys MPT completeness on easy
    # data, and under IMPLIED weights it additionally deepens the ratchet
    # through .IwRatchetDepth()'s targetHits/defaultHits escalation (capped at
    # .iwRatchetMaxCycles), which is a genuine reach lever the equal-weights
    # measurement above cannot see.
    hitMultiplier = if (rung <= 4L) 1L else as.integer(2^(rung - 4L))
  )
}

# Automatic rung, from dataset size and character count.  Returns an INDEX into
# .effortLadder, so `effort = 0` reproduces the previous `strategy = "auto"`
# choice exactly and the default stays size-aware.
# @param nTip Integer number of taxa
# @param nChar Integer number of character patterns (unique columns)
.AutoRung <- function(nTip, nChar) {
  if (nTip <= 30L) return(1L)                       # sprint
  # Few characters -> flat landscape; thorough search is pointless
  if (nChar < 100L) return(2L)                      # default
  # Large trees (>=120 tips): `large` is thorough's provisioning with a raised
  # replicate cap (2026-07-07 rebase; see .StrategyPresets).
  if (nTip >= 120L) return(4L)                      # large
  # Enough characters to have a structured landscape;
  # moderate-to-large datasets benefit from intensive search
  if (nTip >= 65L) return(3L)                       # thorough
  2L                                                # default
}

# Resolve the requested rung.  `effort` is an OFFSET from the automatic choice,
# so that the default (0) is exactly what the package chose before this argument
# existed, on every dataset size.  Clamped to [1, .effortMaxRung]; clamping at
# the bottom means a large negative offset reliably selects `sprint` whatever
# the dataset, which is what most callers wanting "just make it quick" mean.
.EffortRung <- function(autoRung, effort, verbosity = 1L) {
  if (length(effort) != 1L || is.na(effort) || !is.finite(effort) ||
      effort != as.integer(effort)) {
    stop("`effort` must be a single whole number (an offset from the ",
         "automatic setting; 0 keeps it).")
  }
  wanted <- autoRung + as.integer(effort)
  rung <- max(1L, min(.effortMaxRung, wanted))
  if (wanted > .effortMaxRung && verbosity >= 1L) {
    message("`effort` clamped to rung ", .effortMaxRung, " (",
            .RungSpec(.effortMaxRung)[["maxReplicates"]],
            " replicates): the largest replicate budget representable as an ",
            "integer. Set `maxReplicates` and `targetHits` directly if you ",
            "need more.")
  }
  rung
}

#' Find most parsimonious trees
#'
#' `MaximizeParsimony()` performs a multi-replicate driven search for
#' most-parsimonious trees, combining random addition sequence (Wagner)
#' starting trees, tree bisection and reconnection  (\acronym{TBR})
#' rearrangement, exclusive sectorial search (\acronym{XSS}),
#' ratchet perturbation, drift, and tree fusing.
#'
#' The search pipeline follows the "new technology search" approach of
#' \insertCite{Goloboff1999;textual}{TreeSearch}, and resembles the
#' implementation in TNT \insertCite{Goloboff2016}{TreeSearch}.
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
#' @section Completeness of the returned tree set:
#' `MaximizeParsimony()` returns the distinct, fully-resolved optimal
#' topologies held in its tree pool; it does not guarantee that every
#' most-parsimonious tree (\acronym{MPT}) is recovered.
#' The size of the returned set is bounded by, in order:
#' \enumerate{
#'   \item **`poolMaxSize`** (default `100`) — a hard ceiling on the number of
#'     trees retained.  Raise it (via [`SearchControl()`]) to keep more MPTs;
#'     with the default you will never see more than 100.
#'   \item **MPT-enumeration time.** After the main search, a TBR plateau walk
#'     enumerates equal-score neighbours of each pool tree, within a time
#'     reserve of `maxSeconds * enumTimeFraction`.  If this phase times out it
#'     returns a *partial* set (the run reports `stop = "timeout"`); allow more
#'     `maxSeconds`, or raise `enumTimeFraction`, for a more complete set.
#'   \item **TBR-island coverage.** The plateau walk only explores islands of
#'     equal-score trees that a main-loop replicate actually landed on.  MPTs
#'     in unvisited islands are never enumerated, however large `poolMaxSize`
#'     is; increase `maxReplicates` to seed more islands.
#' }
#' By default (`collapse = TRUE`), zero-length (unsupported) branches are
#' contracted into polytomies and the returned set is deduplicated on the
#' resulting collapsed topologies, so `n_topologies` counts distinct *collapsed*
#' topologies — the same convention TNT applies under "collapse zero-length
#' branches".
#' This matters because a single soft polytomy (an unsupported clade of 
#' \eqn{k} taxa) has \eqn{(2k-3)!!}
#' equally-parsimonious binary resolutions, so leaving branches resolved can
#' inflate the apparent number of optimal trees by orders of magnitude without
#' adding any biological information.  Set `collapse = FALSE` to return
#' fully-resolved trees instead (one arbitrary resolution per distinct collapsed
#' topology).
#'
#' Implied weighting is supported natively: set `concavity` to a numeric
#' value (e.g.\sspace{}10).
#' Profile parsimony (`concavity = "profile"`) is supported natively.
#' Inapplicable tokens are treated as ambiguous, and each character is scored
#' by its information profile \insertCite{Faith2001}{TreeSearch}; see
#' [`PrepareDataProfile()`] for how multi-state profiles are computed.
#'
#' @param dataset A phylogenetic data matrix of \pkg{phangorn} class
#' \code{phyDat}, whose names correspond to the labels of any accompanying tree.
#' @param tree (optional) A bifurcating tree of class \code{\link[ape]{phylo}},
#'   or a `multiPhylo` containing a pool of such trees, which must all bear the
#'   same tip labels.
#'   Replicate _i_ starts from tree _i_ of the pool (warm-start), skipping the
#'   random Wagner tree construction; any further replicates begin from random
#'   Wagner trees.  Supplying a single tree thus warm-starts the first
#'   replicate only.
#'   This is useful for continuing a search from previously found optima: a
#'   whole `multiPhylo` of most-parsimonious trees seeds the search with the
#'   topological diversity that tree fusing exploits, which a single tree
#'   cannot.
#'   One tree is consumed per replicate actually run, so a search that
#'   converges early — on `targetHits`, `maxSeconds` or the perturbation
#'   limit, whichever fires first — draws on only part of a large pool, and
#'   says so in a warning.  Raise `targetHits` as well as `maxReplicates` to
#'   use more of it.
#'   If unspecified, all replicates start from random Wagner trees.
#'   Edge lengths are not supported and will be deleted.
#'   Rooted and unrooted trees are both accepted; an unrooted tree is rooted
#'   arbitrarily (on its first tip) before the search begins, which may
#'   affect how any polytomies it contains are resolved.
#' @param concavity Determines the degree to which extra steps beyond the first
#' are penalized.  Specify a numeric value to use implied weighting
#' \insertCite{Goloboff1993}{TreeSearch}; `concavity` specifies _k_ in
#'  _k_ / (_e_ + _k_). A value of 10 is recommended;
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
#'       `hierarchy`.  **Scores are rooting-sensitive**: the step matrix of
#'       this recoding is asymmetric -- gaining the controlling character costs one more
#'       than the number of secondaries it brings into existence, against 1 to
#'       lose it -- so a tree's length depends on where its root sits, whereas
#'       parsimony under the other methods does not.  Lengths are therefore
#'       reported at a canonical rooting, on the first taxon of `dataset`, which
#'       is the rooting the returned trees carry; `TreeLength()` canonicalises
#'       identically, so it reproduces the reported score and one topology has
#'       one length.  That value is an upper bound on the rooting-free minimum,
#'       exceeding it by at most the total number of secondary characters across
#'       hierarchy blocks, and attaining it for 87-98% of rootings in
#'       simulation.}
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
#' @param effort Integer: how much search effort to spend, **relative to the
#'   amount chosen automatically** for this dataset.  `0` (the default) accepts
#'   the automatic choice; `1` asks for one notch more, `-1` one notch less.
#'
#'   The automatic choice is made from dataset size and character count, since
#'   those predict how much search a matrix repays: `sprint` for <=30 taxa;
#'   `large` for >=120 taxa with >=100 character patterns; `thorough` for
#'   65-119 taxa with >=100 character patterns; `default` otherwise.  Because
#'   `effort` is an offset rather than an absolute level, `effort = 0` gives a
#'   30-taxon and a 300-taxon matrix quite different searches -- which is the
#'   intent.
#'
#'   The rungs, in order:
#'   \describe{
#'     \item{1, `sprint`}{Fast: 3 ratchet cycles, no drift, minimal sectorial.
#'       Small datasets and quick surveys.}
#'     \item{2, `default`}{Balanced: 6 ratchet cycles, sectorial search and
#'       fusing.}
#'     \item{3, `thorough`}{Intensive: 20 ratchet cycles, adaptive perturbation,
#'       extra sectorial rounds, drift (2 cycles), 5 Wagner starts and an outer
#'       cycle loop.  The drift cycles also recover equal-score trees on
#'       TBR-disconnected islands that random restarts alone miss.}
#'     \item{4, `large`}{`thorough`'s provisioning with `maxReplicates` raised
#'       to 500, to suit the higher per-replicate cost of big trees.}
#'     \item{5 and up}{`thorough`'s provisioning, with both the replicate budget
#'       and the hit target doubling each notch (1000, 2000, 4000 ...
#'       replicates), so that one notch always means roughly twice the work.
#'       There is no policy ceiling: extra replicates cannot cost reach, only
#'       wall, which is what you asked to spend.  The ladder stops only at rung
#'       26, where the replicate budget outgrows R's integer type.}
#'   }
#'
#'   Above rung 4 the **replicate budget** is what buys reach on hard datasets.
#'   Raising `targetHits` alone does not: it ends a run early on easy datasets,
#'   but on hard ones it is never reached and `maxReplicates` binds first.  It
#'   is raised in step all the same, because it governs when *easy* runs stop --
#'   and under implied weights it additionally deepens the ratchet (see
#'   `targetHits`).
#'
#'   The rung-4 budget of 500 is measured: a 34-matrix, 120--180-tip sweep found
#'   the fraction of runs reaching the best score climbing from 0.68 at 96
#'   replicates to 0.79 at 250, with the hard-matrix subset **still climbing at
#'   500 and no knee**.  The doubling *above* that is an operating point rather
#'   than a fitted constant: nothing measures where the reach curve flattens, so
#'   the ladder simply keeps offering more in even steps.  Treat rungs 5+ as
#'   "spend about twice as long again", not as calibrated levels.
#'
#'   Anything you set yourself wins: `maxReplicates` and `targetHits` you supply
#'   are never rescaled by `effort`, and explicit `control` fields always
#'   override the rung's preset -- for example
#'   `effort = -2, control = SearchControl(ratchetCycles = 10L)` uses the lower
#'   rung's settings for everything except `ratchetCycles`.
#'
#'   Every rung stops on `targetHits` and the `perturbStopFactor`
#'   no-improvement rule; `consensusStableReps` (consensus-stability stopping) is
#'   off by default and no rung enables it.  Under implied weights only, rungs 1
#'   and 2 additionally stop on a flat replicate patience (`stopPatience` 20 and
#'   15 respectively) and deepen the ratchet to match (`ratchetCycles` 12 and
#'   20); the pair is a package, since each half fails on its own.  Equal weights
#'   is unaffected.
#' @param maxReplicates Integer: maximum number of independent search
#'   replicates (default: 96).
#'   The default is a multiple of 48 (= LCM(12, 16)) so that replicates
#'   divide evenly across common 12- or 16-core machines when running in
#'   parallel.
#'   When `effort` resolves to rung 4 (`large` -- chosen automatically for
#'   datasets of \eqn{\ge}{>=} 120 tips and \eqn{\ge}{>=} 100 characters) and
#'   `maxReplicates` is left at its default, the
#'   cap is raised to 500: a 120--180-tip sweep showed the fraction of runs
#'   reaching the best-known score climbing from 0.68 at 96 replicates to 0.79
#'   at 250 and still rising at 500, with no plateau.  Raising the cap only
#'   appends later replicates, so it never delays an earlier improvement, and
#'   easy datasets still stop early once `targetHits` is met; only genuinely
#'   hard datasets run the extra replicates.  An explicit `maxReplicates`
#'   is always respected.
#'   For large or complex datasets a higher value improves the chance of
#'   finding all MPTs.  A rough minimum is
#'   `max(10, ceiling(NTip * NChar / 5000))`, where `NChar = sum(weight)`.
#'   A warning is issued when an explicit value falls below this threshold
#'   for datasets with 30 or more taxa.
#' @param targetHits Integer: stop a replicate series once the best score has
#'   been re-found this many times without further improvement
#'   (default: `max(10, NTip / 5)`).  This is the main control over *how hard the
#'   search tries to be sure it is finished*, and rungs 1-4 of `effort` all
#'   share it -- they differ in per-replicate effort, not in when they stop.
#'   (Rung 5 and above raise it, alongside the replicate budget.)  It sets the balance between the two goals a user may bring to a
#'   search:
#'   \describe{
#'     \item{A single tree one can be reasonably confident is
#'       most-parsimonious}{Use a small `targetHits` (e.g. 4--10).  The search
#'       stops soon after the score stops improving: fast, and safe on datasets
#'       whose optimum is reached early.  On hard datasets the score can still
#'       improve after a long unproductive stretch (a better tree may lie many
#'       replicates away), so a small `targetHits` trades a chance at the true
#'       optimum for speed; raise it (or `maxReplicates`) when certainty matters
#'       more than wall-clock.}
#'     \item{A set of trees representing the full range of most-parsimonious
#'       trees}{Use a large `targetHits` with a high `maxReplicates`.  Distinct
#'       equally-parsimonious topologies -- and whole \acronym{TBR}-disconnected
#'       islands of them -- keep being discovered for as long as replicates run,
#'       and the terminal enumeration step can only fill in trees on islands a
#'       replicate has already reached, so a larger budget samples more islands.
#'       No stopping rule can *detect* that every island has been found: a long
#'       run with no new topology is not proof that none remain, so completeness
#'       is bought with search effort, never inferred.}
#'   }
#'   Under implied weights (finite `concavity`) at `effort` rung 3 (`thorough`)
#'   or above, raising `targetHits` above its default also deepens the ratchet
#'   in proportion, up to 115 cycles: no dataset property reliably predicts how
#'   much character reweighting a matrix needs, so a raised `targetHits` is taken
#'   as the user's own signal that this one needs more.  Lowering `targetHits`
#'   does not make the ratchet shallower than its default depth (fewer cycles
#'   were slower to the optimum on every matrix tested), and setting
#'   `ratchetCycles` yourself overrides this entirely.
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
#'   At `1` (default) each replicate reports its score, pool size and hit count;
#'   at `2` and above each search phase reports on completion.
#'
#'   On a large dataset a single phase can run for many minutes, during which
#'   neither level would print anything: on a 182-tip, 420-character matrix with
#'   inapplicable tokens throughout, one TBR phase took 582 s and one ratchet
#'   549 s, together 96% of a 1173 s replicate.  A *heartbeat* therefore reports
#'   from inside the long phases -- overwriting one console line at a terminal,
#'   or emitting discrete lines to a batch log -- so a slow search is
#'   distinguishable from a hung one.  See the environment variables below.
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
#'   Note that supplying a callback *replaces* the per-replicate console line
#'   rather than adding to it, and that the callback fires only when a replicate
#'   or a fuse completes -- so on a dataset whose replicates take many minutes,
#'   nothing arrives until the first one finishes.  The heartbeat described under
#'   `verbosity` is independent of the callback and reports throughout.
#' @section Progress reporting in non-interactive sessions:
#'
#' The automatic `cli` progress bar requires an interactive session.  Under
#' `Rscript` (including a batch or cluster job) two environment variables control
#' reporting instead:
#'
#' \describe{
#'   \item{`TREESEARCH_PROGRESS_FILE`}{Path to a status file.  After each
#'     replicate, a single line is written -- and the file truncated, so it
#'     always holds current state rather than a history -- containing
#'     `replicate`, `max_replicates`, `best_score`, `hits_to_best` and
#'     `target_hits`, space-separated.  Poll it to monitor a long job:
#'     `TREESEARCH_PROGRESS_FILE=progress.txt Rscript analysis.R`.  Only
#'     consulted when `progressCallback` is `NULL`.}
#'   \item{`TS_HEARTBEAT_SECONDS`}{Heartbeat cadence in seconds; fractional
#'     values are allowed.  Defaults to 30 at a terminal (where the line
#'     overwrites itself) and 120 to a batch log (where every heartbeat is a
#'     permanent line).  Set to `0` to disable.  An unparseable value falls back
#'     to the default rather than disabling.}
#' }
#'
#' The heartbeat reports only from searches of the whole tree under the real
#' character weights.  Sectorial searches score a subtree, and the ratchet's
#' perturbation phase scores a reweighted matrix; both legitimately run far below
#' the true optimum, so reporting them would look like erratic progress.  Any
#' score the heartbeat prints is therefore directly comparable with the final
#' tree score.
#' @param .rung Internal.  Pins a named entry of the effort ladder
#'   (`"sprint"`, `"default"`, `"thorough"`, `"large"`), or `"none"` to apply no
#'   preset at all; `NULL` (default) selects the rung from `effort` and the
#'   dataset's size, which is what every ordinary call should do.  Exists for
#'   controlled experiments and the preset smoke tests, which need to name a rung
#'   absolutely rather than relative to the automatic choice.  Not part of the
#'   stable interface: prefer `effort`.
#' @param control A [`SearchControl`] object (or a named list) of low-level
#'   search parameters.  Most users can rely on `effort` and
#'   ignore this argument; see [`SearchControl()`] for full documentation
#'   of individual fields.
#' @param collapse Logical: if `TRUE` (default), contract zero-length
#'   (unsupported) branches in the returned trees into polytomies before
#'   returning, and de-duplicate the result on the resulting collapsed
#'   topologies, akin to TNT's "collapse zero-length branches".  A branch is
#'   treated as zero-length when it has minimum possible length 0 (there exists
#'   a most-parsimonious reconstruction with no change along it), evaluated
#'   under the same scoring method used for the search.  `n_topologies` then
#'   counts distinct collapsed topologies, which is comparable across programs.
#'   This is the recommended behaviour: a fully-resolved tree containing a
#'   zero-length branch asserts a grouping the data do not support, and a single
#'   soft polytomy can otherwise inflate the apparent number of optimal trees by
#'   orders of magnitude.  Set `FALSE` to return fully-resolved trees instead
#'   (one arbitrary resolution per distinct collapsed topology), e.g. when a
#'   downstream step requires binary trees.
#'   Collapsing is applied to the best-score trees (the MPTs); any suboptimal
#'   pool trees retained via `poolSuboptimal` are omitted from the collapsed set.
#'   When a `constraint` is supplied, its enforced splits are protected from
#'   collapse, so an enforced-but-unsupported clade (a zero-length branch) stays
#'   visible (a constraint encodes external evidence the matrix does not
#'   capture); unsupported non-constraint branches still collapse.
#' @param ... Backward compatibility.
#'
#' @return A `multiPhylo` object containing the best tree(s) found, with
#'   attributes:
#'   \describe{
#'     \item{`score`}{Best parsimony score.}
#'     \item{`replicates`}{Number of replicates completed.}
#'     \item{`hits_to_best`}{Number of independent discoveries of the best
#'       score.}
#'     \item{`n_topologies`}{Number of distinct best-score topologies returned.
#'       With `collapse = TRUE` (default) this counts distinct *collapsed*
#'       topologies (equal to `length()` of the result); with `collapse = FALSE`
#'       it is the number of distinct fully-resolved topologies in the pool at
#'       the best score.}
#'     \item{`last_improved_rep`}{1-based index of the replicate that last
#'       improved the best score (0 if not tracked, e.g. parallel search).}
#'     \item{`timed_out`}{Logical: `TRUE` if the search stopped because
#'       `maxSeconds` was exceeded.}
#'     \item{`consensus_stable`}{Logical: `TRUE` if the search stopped
#'       because the strict consensus was unchanged for
#'       `consensusStableReps` consecutive replicates.}
#'     \item{`perturb_stop`}{Logical: `TRUE` if the search stopped because a
#'       run of replicates failed to improve the best score -- either the
#'       `nTip * perturbStopFactor` dry-spell limit or the flat `stopPatience`
#'       count (see [`SearchControl()`]).  The flag does not distinguish which of
#'       the two fired; in a serial search, comparing
#'       `last_improved_rep + stopPatience` against `replicates` will tell you.}
#'     \item{`timings`}{Named numeric vector of cumulative wall-clock time
#'       (in milliseconds) spent in each search phase across all replicates:
#'       `wagner_ms`, `tbr_ms`, `xss_ms`, `rss_ms`, `css_ms`, `ratchet_ms`,
#'       `drift_ms`, `final_tbr_ms`, `fuse_ms`, `nni_ms`, `nni_perturb_ms`,
#'       `anneal_ms`, `prune_reinsert_ms`.}
#'     \item{`replicate_scores`}{Numeric vector of the best parsimony score
#'       found by each completed replicate.  Passed to [ScoreSpectrum()] for
#'       Chao1-style landscape coverage estimation.}
#'     \item{`candidates_evaluated`}{Number of TBR/SPR-class candidate
#'       rearrangements evaluated across the whole search — the analogue of
#'       TNT's "rearrangements examined", useful for comparing search
#'       efficiency (candidates per unit of score improvement).  Counted only
#'       for single-threaded searches (`0` for any parallel search, i.e.
#'       `nThreads != 1`, including `nThreads = 0` auto-detect); excludes
#'       NNI-warmup and simulated-annealing candidates.}
#'   }
#'
#' @examples
#' data("inapplicable.phyData", package = "TreeSearch")
#' dataset <- inapplicable.phyData[["Vinther2008"]]
#' result <- MaximizeParsimony(
#'   dataset,
#'   inapp = "missing",
#'   maxReplicates = 12L,
#'   targetHits = 4L
#' )
#' result
#' attr(result, "score")
#'
#' # Ask for one notch more search than this dataset would get by default,
#' # whatever its size:
#' harder <- MaximizeParsimony(dataset, effort = 1L, maxReplicates = 12L)
#'
#' @template MRS
#' @family tree scoring
#' @seealso [`Resample()`] for jackknife and bootstrap resampling.
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
    effort = 0L,
    maxReplicates = 96L,
    targetHits = NULL,
    maxSeconds = 0,
    nThreads = 1L,
    verbosity = 1L,
    progressCallback = NULL,
    control = SearchControl(),
    collapse = TRUE,
    .rung = NULL,
    ...
) {

  # --- Input validation: check dataset first ---
  if (is.null(dataset)) {
    stop("`dataset` cannot be NULL.")
  }

  # Record whether the user explicitly supplied `maxReplicates` BEFORE any
  # reassignment: assigning to the formal (e.g. the strategy-scaled default
  # below) would immediately flip `missing()`, so this top-of-body capture is
  # the only reliable read.
  userSetReps <- !missing(maxReplicates)

  # `maxReplicates < 1` runs the search loop zero times: the pool stays
  # empty, `best_score` never leaves its C++ sentinel of -1, and the
  # empty-pool fallback below would silently return the random starting
  # tree tagged with that bogus score instead of throwing an error.
  if (length(maxReplicates) != 1L || is.na(maxReplicates) ||
      as.integer(maxReplicates) < 1L) {
    stop("`maxReplicates` must be a single integer of at least 1.")
  }

  # --- Set targetHits default if not provided ---
  # `defaultHits` is retained even when the user supplies `targetHits`: the
  # implied-weights ratchet depth below scales with the user's *escalation*
  # (targetHits / defaultHits), not with the absolute value.
  defaultHits <- max(10L, as.integer(NTip(dataset) / 5))
  # Captured before the assignment below, for the same reason as `userSetReps`:
  # the effort ladder must not scale a hit target the user chose themselves.
  userSetHits <- !is.null(targetHits)
  if (is.null(targetHits)) {
    targetHits <- defaultHits
  }

  # --- Backward compatibility: intercept maxTime → maxSeconds ---
  dots <- list(...)
  if ("maxTime" %in% names(dots)) {
    if (missing(maxSeconds) || maxSeconds == 0) {
      maxSeconds <- as.double(dots[["maxTime"]])
    }
    .Deprecated(msg = paste0(
      "Use `maxSeconds` instead of `maxTime` in MaximizeParsimony().",
    ))
    dots[["maxTime"]] <- NULL
  }

  # --- Reject legacy parameters ---
  .morphyParams <- c("ratchIter", "tbrIter", "startIter", "finalIter",
                     "maxHits", "quickHits", "ratchEW", "tolerance")
  legacyHits <- intersect(names(dots), .morphyParams)
  if (length(legacyHits)) {
    stop("Parameter", if (length(legacyHits) > 1L) "s", " ",
         paste0(sQuote(legacyHits), collapse = ", "),
         if (length(legacyHits) == 1L) "are" else "were",
         " discontinued in v2.0.0.\n",
         "  Use this function's own controls instead ",
         "(see `?SearchControl`, `maxReplicates`, `maxSeconds`).",
         call. = FALSE)
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

  # --- Resolve the effort rung ---
  # `effort` is an OFFSET from the automatic choice, not an absolute level, so
  # `effort = 0` reproduces the size-aware selection exactly on every dataset
  # size -- the previous `strategy = "auto"` behaviour, unchanged.
  #
  # `.rung` is INTERNAL (leading dot): it pins a named menu entry, or "none" to
  # apply no preset at all, which controlled experiments and the preset smoke
  # tests need.  Deliberately not user-facing: users ask for effort relative to
  # the automatic choice, and "no preset at all" must not be reachable by an
  # accidentally-missing variable propagating in.
  autoRung <- .AutoRung(NTip(dataset), sum(attr(dataset, "weight")))
  if (is.null(.rung)) {
    rung <- .EffortRung(autoRung, effort, verbosity)
    rungName <- .RungSpec(rung)[["preset"]]
  } else if (identical(.rung, "none")) {
    rung <- NA_integer_
    rungName <- "none"
  } else {
    rung <- match(.rung, .effortLadder)
    if (is.na(rung)) {
      stop("Internal `.rung` must be one of ",
           paste(sQuote(.effortLadder), collapse = ", "), ", or \"none\".")
    }
    rungName <- .rung
  }

  # --- Apply the rung ---
  if (!identical(rungName, "none")) {
    spec <- .RungSpec(rung)
    strategy <- spec[["preset"]]        # menu label, used by the IW packages below
    preset <- .StrategyPresets()[[strategy]]
    {
      control <- .ApplyStrategyPreset(control, preset, names(controlDots))
      if (verbosity >= 1L) {
        cli::cli_alert_info(
          "Effort {.strong {effort}}: {.emph {strategy}}, rung {rung}"
        )
      }
      # Rung-scaled replicate cap.  The value comes from .RungSpec() -- the ONE
      # place a rung's replicate cap is defined -- rather than a switch on the
      # preset name, so rung 4 cannot end up with two disagreeing sources.
      # Rung 4 (the `large` band, >=120 tips) needs many more independent
      # restarts than the 96 default to reliably reach the optimum: a 34-matrix
      # 120-180t sweep found reach@96 = 0.68 climbing to reach@250 = 0.79, with
      # the hard-matrix subset still climbing at 500 and no knee.  Only override
      # when the user did not set `maxReplicates` themselves.
      if (!userSetReps && !is.na(spec[["maxReplicates"]])) {
        maxReplicates <- spec[["maxReplicates"]]
      }

      # Rung-scaled hit target (rung 5 and up).  Applied HERE, before
      # .IwRatchetDepth() below, because that reads `targetHits / defaultHits`
      # as the user's escalation signal -- so an effort-raised hit target also
      # deepens the implied-weights ratchet, which is the point.  Skipped when
      # the user named `targetHits` themselves: their number is a statement
      # about this dataset and outranks the ladder.
      if (!userSetHits && spec[["hitMultiplier"]] > 1L) {
        targetHits <- as.integer(targetHits * spec[["hitMultiplier"]])
      }

      # Implied-weights ratchet depth. Under implied weights the optimum often
      # sits in a small basin at fine score resolution, separated from an
      # easy-to-find near-optimum by a fraction of a step; character reweighting
      # (the ratchet) is what crosses that gap, and extra *replicates* cannot
      # substitute for it: on one 106-tip matrix 20 000 random-addition restarts
      # all plateau above the optimum that a deeper ratchet reaches.  A 36-matrix
      # grid over
      # `ratchetCycles` in {6, 12, 20, 48, 96} (implied weights, k = 10) found
      # expected wall-clock-to-optimum minimised at 48: on the 4 cycle-sensitive
      # matrices the mean fell 1435 s -> 709 s, while the 32 others paid a median
      # +0.2 s with reach unchanged.  The curve is flat from ~20 to ~96 and rises
      # steeply below 20, so 48 is a broad optimum rather than a knife-edge --
      # hence a constant, not a per-dataset function: dataset size cannot target
      # the need (94-, 106- and 110-tip matrices each appear as both
      # cycle-sensitive and insensitive), and a Wagner-tree consistency gate,
      # though it does correlate with the need, beats the constant by nothing
      # once the constant sits in the flat region.
      #
      # `targetHits` is the user's own statement of how hard this dataset is, so
      # raising it deepens the ratchet in proportion -- the one signal available
      # that no dataset feature supplies.  Escalation only: de-escalating (the
      # documented `targetHits = 4` "one tree, quickly" idiom) must not drop
      # below 48, since fewer cycles were slower for *every* stratum measured.
      # Capped at the largest depth actually tested.
      #
      # Equal weights is excluded deliberately: the same 3-arm test over 68
      # matrices found no reach gain there (0.970 vs 0.965) for a small wall
      # cost, the integer landscape lacking the fractional basins this escapes.
      # Scoped to `thorough`/`large`, whose other knobs match the grid; `default`
      # and `sprint` co-tuned their ratchet with different sectorial settings and
      # are untouched.
      iwCycles <- .IwRatchetDepth(
        strategy, concavity, targetHits, defaultHits,
        userSet = union(names(controlDots), attr(control, "explicit"))
      )
      if (!is.null(iwCycles)) {
        control[["ratchetCycles"]] <- iwCycles
      }

      # Implied-weights operating point for `sprint` and `default`: a deeper
      # ratchet, paid for by a flat replicate patience (`stopPatience`).
      #
      # The two knobs are a package because each fails the other's gate alone.
      # The ratchet is the quality lever: on `default` (44 training matrices,
      # 65-385 tips, 6 seeds, k = 10) `ratchetCycles = 20` alone scored better on
      # 11 matrices and worse on 0 (p = 0.001) and raised reach 0.78 -> 0.84, but
      # cost +25 s of a 151 s mean (36 matrices slower, p = 2.5e-05).  Patience
      # is the wall lever, and alone it degrades score (`sprint` 0/6; `default`
      # 1 better/15 worse, p = 5e-04): stopping early without deepening the
      # replicate simply searches less.  Together, at the values below:
      #   sprint   median matrix -26% wall (19 faster/5 slower), score 4/0,
      #            distinct MPTs unchanged, reach 0.847 -> 0.861
      #   default  median matrix -18% wall (33/11, p = 0.001), score 9/3 --
      #            a favourable direction only, NOT significant (p = 0.15)
      # A 6-arm sweep over patience {10, 15, 20, 25, 30} (2448 cells) found score
      # and wall both MONOTONE in the value with no spike at any of them, so
      # these are operating points chosen on a smooth trade-off, not fitted
      # constants: loosening patience buys score and gives back wall.  Values
      # were selected against `auto`'s regression-averse objective, i.e. on the
      # MEDIAN per-matrix wall change and its sign count, not the mean -- the
      # mean is dominated by the largest matrices and reverses the choice.
      # Residual cost, deliberately accepted and worth stating plainly: 9 of 44
      # `default` matrices are still >10% slower (worst +110%), those where
      # patience does not bite and the deeper ratchet is not paid for.
      #
      # All of it was measured with `nThreads = 1`.  The parallel path implements
      # the same rule over the shared pool but evaluates it on the coordinating
      # thread's poll, so patience bites later there and the wall saving will be
      # smaller; the deeper ratchet applies unchanged either way.
      #
      # `default` sets `adaptiveLevel = TRUE`, so 20 is a BASE that the hit-rate
      # rescaler moves within ~10-30 at runtime; the measured arm had exactly
      # that, so this matches its measurement -- do not "fix" it to a fixed 20.
      #
      # Equal weights and profile parsimony are excluded: neither was measured.
      # `thorough`/`large` are excluded for the same reason, and take their own
      # implied-weights depth from .IwRatchetDepth() above.
      iwStop <- .IwStopPackage(
        strategy, concavity,
        userSet = union(names(controlDots), attr(control, "explicit"))
      )
      for (.f in names(iwStop)) {
        control[[.f]] <- iwStop[[.f]]
      }
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

  # --- Normalize `concavity` ---
  # Route the profile-mode test through the same lenient matcher used at the
  # scoring entry points (`.UseProfile()`, called from tree_length.R and
  # PolEscapa.R) so concavity = "Profile" or "prof" search in profile mode
  # exactly as later re-scoring the result would.  Everything else must
  # resolve to a valid positive number (or Inf) *here*: letting a bad string
  # such as "10" reach as.double() unchecked would silently coerce to 10 while
  # leaving IW's min_steps unpopulated downstream, so the C++ engine would run
  # IW uncorrected for homoplasy with no error or warning (see min_steps in
  # ts_data.cpp / ts_fitch.cpp).
  useProfile <- !missing(concavity) && .UseProfile(concavity)
  if (!useProfile) {
    rawConcavity <- concavity
    concavity <- suppressWarnings(as.numeric(concavity))
    if (length(concavity) != 1L || is.na(concavity)) {
      stop("`concavity` must be a single positive number, Inf (for equal ",
           "weights), or \"profile\" (for profile parsimony); got ",
           deparse(rawConcavity), ".")
    }
  }

  # --- Profile parsimony: prepare data ---
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

  # --- Starting tree(s) ---
  # `tree` may be a single `phylo` or a `multiPhylo` holding a whole pool of
  # warm starts: replicate i then begins from tree i, and replicates beyond
  # the pool build random Wagner trees as usual.  Resuming from a previous
  # run's MPTs is the motivating case -- the pool's topological diversity is
  # exactly what the fusing machinery needs, and one tree cannot supply it.
  userTree <- !missing(tree) && !is.null(tree)
  if (!userTree) {
    tree <- TreeTools::RandomTree(nTip, root = TRUE)
    tree[["tip.label"]] <- names(dataset)
    startTrees <- list(tree)
  } else if (inherits(tree, "multiPhylo")) {
    # `[[` rather than unclass(): a compressed `multiPhylo` stores tip labels
    # once in a shared `TipLabel` attribute, and only `[[` restores them.
    startTrees <- lapply(seq_along(tree), function(i) tree[[i]])
    if (length(startTrees) == 0L) {
      stop("`tree` contains no trees.")
    }
  } else {
    startTrees <- list(tree)
  }
  if (!all(vapply(startTrees, inherits, logical(1), "phylo"))) {
    stop("`tree` must be of class 'phylo'.")
  }
  if (length(startTrees) > 1L) {
    refLabels <- sort(startTrees[[1L]][["tip.label"]])
    sameTips <- vapply(startTrees[-1L], function(x) {
      identical(sort(x[["tip.label"]]), refLabels)
    }, logical(1))
    if (!all(sameTips)) {
      stop("All trees in `tree` must bear the same tip labels.")
    }
  }

  # --- Match tree tips to dataset ---
  # Every starting tree shares a tip set, so resolve the mismatch once.
  leaves <- startTrees[[1L]][["tip.label"]]
  taxa <- names(dataset)
  treeOnly <- setdiff(leaves, taxa)
  datOnly <- setdiff(taxa, leaves)
  if (length(treeOnly)) {
    warning("Dropping taxa on tree but not in dataset: ",
            paste0(treeOnly, collapse = ", "))
  }
  if (length(datOnly)) {
    warning("Dropping taxa in dataset but not on tree: ",
            paste0(datOnly, collapse = ", "))
    dataset <- dataset[-match(datOnly, taxa)]
  }

  # Normalize each start into the form the C++ engine expects.
  startTrees <- lapply(seq_along(startTrees), function(i) {
    tr <- startTrees[[i]]

    # Reject a structurally invalid `phylo` before any traversal code sees it.
    # These objects are not exotic: ape::unroot() accepts TreeTools' `order =
    # "preorder"` attribute and then mishandles it, so unrooting any TreeTools
    # tree yields an edge matrix carrying NA entries.  Rooting or reordering
    # one segfaults inside the dependency, below the level at which R can
    # catch anything, so the guard has to sit ahead of the repair block.
    .CheckStartTree(tr, if (length(startTrees) > 1L) i else NA_integer_)

    # Root before checking for bifurcation: MakeTreeBinary() assumes a rooted
    # tree, where the root's "effective" degree needs +1 for its absent
    # parent edge.  Applied to an already-unrooted tree, that +1 misreads the
    # root's legitimate degree-3 trifurcation as a polytomy and inserts a
    # spurious node.  TreeTools::TreeIsRooted() is used (not ape::is.rooted(),
    # which returns NA for some valid trees here and would break this `if`).
    if (!TreeTools::TreeIsRooted(tr)) {
      tr <- RootTree(tr, 1L)
    }

    # Make bifurcating if needed
    if (dim(tr[["edge"]])[1] != 2L * tr[["Nnode"]]) {
      tr <- MakeTreeBinary(tr)
      # Re-check: MakeTreeBinary() can itself return a malformed object, and
      # the RootTree() below is exactly where such an object kills the session.
      .CheckStartTree(tr, if (length(startTrees) > 1L) i else NA_integer_)
      if (dim(tr[["edge"]])[1] != 2L * tr[["Nnode"]]) {
        tr <- RootTree(tr, 1L)
      }
      if (dim(tr[["edge"]])[1] != 2L * tr[["Nnode"]]) {
        stop("Could not make `tree` binary.")
      }
    }
    if (length(treeOnly)) {
      tr <- TreeTools::DropTip(tr, treeOnly)
    }

    # Reorder tips to match dataset, put in preorder
    tr <- Preorder(RenumberTips(tr, names(dataset)))

    # Ensure root's first child is a tip (for C++ engine compatibility)
    if (tr[["edge"]][1L, 2L] > NTip(tr)) {
      tr <- RootTree(tr, 1L)
    }
    tr
  })
  tree <- startTrees[[1L]]


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
  if (userSetReps && nTip >= 30L && verbosity > 0L) {
    # `weight` here is the .ScaleWeight()-integerised value used by the C++
    # engine (up to ~1260x the original for fractional weights); the
    # recommendation formula is about the number of characters in the
    # dataset, so it must read `at$weight` (pre-scaling) rather than `weight`.
    nChars <- sum(at$weight)
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
    startEdge = if (userTree) lapply(startTrees, `[[`, "edge") else NULL,
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

  # A pool is consumed one tree per replicate *run*, which is bounded by
  # whichever stopping rule fires first -- usually `targetHits`, not
  # `maxReplicates`.  Only the completed count is a truthful bound, so report
  # it after the fact rather than guessing beforehand.  Ungated by `verbosity`,
  # like the taxon-dropping warnings above: silently ignoring supplied data
  # warrants a warning however quiet the search itself is.
  if (length(startTrees) > 1L && result$replicates < length(startTrees)) {
    warning("Used ", result$replicates, " of the ", length(startTrees),
            " trees supplied to `tree`: the search ran ", result$replicates,
            " replicate", if (result$replicates == 1L) "" else "s",
            " and each starts from one tree. Raise `targetHits` or ",
            "`maxReplicates` to draw on more of the pool.", call. = FALSE)
  }

  # --- Reconstruct phylo from edge matrices ---
  treeTpl <- tree
  treeTpl[["edge.length"]] <- NULL
  resultTrees <- result$trees
  if (length(resultTrees) == 0L) {
    resultTrees <- list()
  }
  nTopologies <- result$n_topologies
  if (isTRUE(collapse) && length(resultTrees) > 0L) {
    # Contract zero-length (unsupported) branches into polytomies, à la TNT's
    # "collapse zero-length branches" -- done entirely in C++ (ts_collapse_pool)
    # to avoid a per-tree R surgery quagmire.  The kernel re-roots each tree on
    # tip 0 (so root-adjacent edges are trivial -> rooting-invariant *contraction*;
    # note the LENGTH is not rooting-invariant under HSJ/XFORM, T-374, which is
    # why the XFORM pool is rescored at this rooting below),
    # flags aggressive (min-length-0) internal edges in the *search's* scoring
    # mode, contracts them, and deduplicates on the collapsed topology.
    #
    # Collapse the MPTs (best-score trees) only.  A collapsed topology has a
    # unique min-resolution length, so a suboptimal pool tree (poolSuboptimal > 0)
    # can never share a collapsed shape with a best-score tree; restricting to the
    # best score keeps n_topologies's documented "at the best score" meaning (and
    # matches the collapse = FALSE count when no branch is unsupported).
    # result$scores aligns with result$trees; default poolSuboptimal = 0 keeps
    # every tree.  Edge matrices come straight from the engine, so tip i already
    # maps to tip_data row i -- no R rerooting (which would permute tips and
    # mis-score against tip_data; see na-validation-alignment-gotcha).
    nTip <- length(treeTpl[["tip.label"]])
    bestTrees <- resultTrees[result$scores == result$best_score]
    # Under a constraint, protect the enforced splits from collapse ("show the
    # enforced clade"): a constraint is external evidence for a grouping the
    # matrix doesn't capture, so it stays visible even at zero length, while the
    # unsupported non-constraint branches still collapse.  consSplitMatrix rows
    # are the enforced bipartitions in tip_data order (see .PrepareConstraint).
    consSplits <- if (!is.null(constraintConfig)) {
      constraintConfig[["consSplitMatrix"]]
    }
    collapsed <- ts_collapse_pool(
      bestTrees, contrast, tip_data, weight, levels,
      scoringConfig, hsjConfig, xformConfig, consSplits
    )
    outTrees <- lapply(collapsed$trees, function(edgeMat) {
      tr <- list(
        edge = edgeMat,
        Nnode = max(edgeMat) - nTip,        # contiguous ids: max id = nTip + Nnode
        tip.label = treeTpl[["tip.label"]]
      )
      class(tr) <- "phylo"
      Renumber(tr)
    })
    nTopologies <- collapsed$n_topologies
  } else {
    outTrees <- lapply(resultTrees, function(edgeMat) {
      tr <- treeTpl
      tr[["edge"]] <- edgeMat
      # C++ edge order may differ from template; renumber to valid preorder
      Renumber(tr)
    })
  }
  if (length(outTrees) == 0L) {
    outTrees <- list(treeTpl)
  }

  # --- XFORM: report the score of the tree we are actually returning ---
  # `result$best_score` is recorded mid-search at whatever rooting the replicate
  # held.  XFORM's step matrix is asymmetric, so the score is rooting-dependent,
  # and `ts_collapse_pool()` above hands back every tree re-rooted on tip 0.
  # Reporting `best_score` therefore gives the user a number that `TreeLength()`
  # of the returned tree does not reproduce -- measured at 178 reported against
  # 183 returned (T-385; repro in dev/red-team/heavy-tests/).  Rescore the
  # returned pool at the canonical rooting instead: |pool| evaluations, negligible
  # against a search, and `TreeLength()` canonicalises identically, so the two
  # agree by construction.
  #
  # This deliberately does NOT change what the search optimises.  The reported
  # value stays a rooting-dependent upper bound on the min-over-rootings
  # objective, exceeding it by at most the sum of `nSec` over hierarchy blocks
  # (measured: attained by 87-98% of rootings, mean overstatement 0.02-0.17
  # steps).
  #
  # Min-over-rootings reporting -- the variant the plan calls better -- is NOT used.
  # It reports a quantity the search never compared, but the deciding objection is
  # cost: evaluated naively it is (2 * nTip - 3) x on the Sankoff term, so a
  # 4000-tip pool of 100 trees would need ~800k evaluations at the boundary.
  # Doing it affordably needs an all-rootings up-down DP, which is its own piece
  # of work (Option 4), not a line in a reporting fix.
  # See dev/plans/2026-07-29-t374b-xform-rooting-policy.md (Option 3).
  bestScore <- result$best_score
  if (useXform && length(outTrees) > 0L) {
    canonicalScores <- TreeLength(
      structure(outTrees, class = "multiPhylo"),
      dataset, inapplicable = "xform", hierarchy = hierarchy
    )
    bestScore <- min(canonicalScores)
    if (diff(range(canonicalScores)) > sqrt(.Machine$double.eps)) {
      # Pool membership is chosen on search-time scores taken at differing
      # rootings (`result$scores` above), so trees held to be equally
      # parsimonious can differ once scored at one rooting.  Not silently
      # averaged away: this is the open residue of T-374, and staying quiet about
      # it is what let the reporting gap survive this long.
      warning("Returned trees do not share a length at a common rooting (",
              paste(signif(range(canonicalScores), 8), collapse = " to "),
              "); reporting the smallest.  The x-transformation's score is ",
              "rooting-dependent -- see ?MaximizeParsimony.")
    }
  }

  # --- Output ---
  if (verbosity > 0L) {
    total_s <- round(sum(unlist(result$timings), na.rm = TRUE) / 1000, 1)
    stop_reason <- if (isTRUE(result$timed_out)) "timeout"
                   else if (isTRUE(result$consensus_stable)) "consensus stable"
                   else if (isTRUE(result$perturb_stop)) "perturbation limit"
                   else "replicate limit"
    cli_alert_success(paste0(
      "Search complete: score {.strong {signif(bestScore, 7)}}, ",
      "{result$replicates} replicate{?s} ",
      "(last improved: #{result$last_improved_rep}), ",
      "{result$hits_to_best} hit{?s} to best, ",
      "{nTopologies} MPT{?s}, ",
      "stop: {stop_reason}, {total_s}s"
    ))
  }

  structure(
    outTrees,
    score = bestScore,
    replicates = result$replicates,
    hits_to_best = result$hits_to_best,
    n_topologies = nTopologies,
    last_improved_rep = result$last_improved_rep,
    timed_out = isTRUE(result$timed_out),
    consensus_stable = isTRUE(result$consensus_stable),
    perturb_stop = isTRUE(result$perturb_stop),
    timings = unlist(result$timings),
    strategy_diagnostics = result$strategy_diagnostics,
    replicate_scores = result$replicate_scores,
    candidates_evaluated = result$candidates_evaluated,
    # NA-certification counters (`exact_verify_sweep` calls executed vs skipped
    # by `TBRParams::certify_unrooted`).  Diagnostic: `naDiag$n_evs_skipped` is
    # the only way an A/B can prove the certification gate reached a live call
    # site, since `do_reroot` needs `tabuSize == 0` -- which the shipped presets
    # do not set.  Timing fields are populated only under `TS_NA_TIMING`, and a
    # threaded run reports the main thread's copy only (each worker owns a
    # private dataset), so read them from serial runs.
    naDiag = result$na_diag,
    class = "multiPhylo"
  )
}

#' Launch tree search graphical user interface
#'
#' Opens a "shiny" app for interactive parsimony tree search and results
#' exploration.
#'
#' @return Opens a Shiny application; does not return a value.
#' @seealso [`MaximizeParsimony()`]
#' @importFrom TreeDist ClusteringInfoDistance
#' @export
EasyTrees <- function () {#nocov start
  needed <- c("cluster", "future", "PlotTools", "promises",
              "protoclust", "Rogue", "shiny", "shinyjs")
  missing <- needed[!vapply(needed, requireNamespace,
                            logical(1L), quietly = TRUE)]
  if (length(missing)) {
    stop("EasyTrees() requires additional packages: ",
         paste(missing, collapse = ", "), ".\n",
         "Install with: install.packages(",
         paste0("\"", missing, "\"", collapse = ", "), ")",
         call. = FALSE)
  }
  shiny::runApp(system.file("Parsimony", package = "TreeSearch"))
}

#' @rdname EasyTrees
#' @export
EasyTreesy <- EasyTrees
#nocov end

.UseProfile <- function (concavity) {
  pmatch(tolower(concavity), "profile", -1L) == 1L
}
