# Basin-diversity instrumentation harness  (audit item B1, 2026-06-22)
# =====================================================================
# Tests whether the QUALITY gap (e.g. Zanol2014: TS ~1265 vs TNT 1261) is
# bounded by the single-tree-multistart architecture (audit B1) or by
# something else.  See dev/plans/2026-06-22-architecture-assumptions-audit.md.
#
# The interpretation is GATED on score-reaching, not diversity alone
# (diversity is neither good nor bad on its own):
#
#                         most restarts reach optimum   |  most miss optimum
#   low distinct-count    reliable convergence (GOOD;    |  stuck in one bad
#                         a diverse set adds nothing)    |  basin -> per-restart
#                                                        |  STRENGTH problem
#   high distinct-count   (rare) wide optimal island     |  scatter across
#                                                        |  suboptimal basins ->
#                                                        |  B1 *plausible* iff the
#                                                        |  optimum recombines them
#
# The DECISIVE probe is the fuse on/off A/B (tree-fusing = recombination of
# pool features = exactly B1's proposed mechanism), all pure config:
#   * fuse-OFF misses, fuse-ON reaches  -> recombination is the active
#       ingredient -> B1 real, a stronger diverse-set tour is worth building.
#   * both miss but TNT reaches         -> pool not diverse enough for fusing
#       to have material (diversity bottleneck) OR optimum is not a
#       recombination (per-restart strength).  INDEP diversity disambiguates.
#   * INDEP already reaches optimum     -> B1 is not bounding this config.
#
# Three TS arms (recipe byte-identical except the sharing knobs):
#   INDEP    K x MaximizeParsimony(maxReplicates = 1)  -- no cross-restart
#            sharing; one full-recipe restart per seed.  THE architecture probe.
#   PROD-ON  1 x MaximizeParsimony(maxReplicates = R)  -- default fuseInterval.
#   PROD-OFF 1 x MaximizeParsimony(maxReplicates = R)  -- fusing disabled.
#   (TNT     xmult, hold many: separate OPTIMAL-ISLAND axis, not per-restart.)
#
# Usage (env-driven, mirrors dev/benchmarks/bench_cell.R conventions):
#   BD_MODE=smoke   Rscript dev/benchmarks/basin_diversity.R
#       small end-to-end run on the local machine: plumbing + gap sanity check.
#   BD_MODE=cell    Rscript dev/benchmarks/basin_diversity.R <cellIndex>
#       run ONE grid cell (one INDEP seed or one PROD variant/seed), write a
#       partial .rds.  This is what the Hamilton array fans out over.
#   BD_MODE=analyze Rscript dev/benchmarks/basin_diversity.R
#       merge the partials, compute metrics, write summary + figure + verdict.
#
# Env knobs: TS_LIB, TNT_EXE, BD_DATASET, BD_BEST, BD_INDEP_K, BD_PROD_REPS,
#            BD_PROD_SEEDS, BD_PARTIAL_DIR, BD_OUT_PREFIX, BD_TNT (0/1),
#            BD_TNT_HITS, BD_TNT_REPLIC, BD_NPERM.

suppressMessages({
  library(TreeSearch,
          lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"), winslash = "/",
                                  mustWork = FALSE))
  library(TreeTools)
  library(TreeDist)
})

# --- Best-known Fitch (EW) scores: apples-to-apples TNT targets ---------------
# SOURCE OF TRUTH: dev/benchmarks/headtohead_phase0.csv, column `tnt` (TNT's best
# Fitch score under the SAME "-"->"?" regime this harness uses).  Read at load so
# the targets can never silently drift from the canonical table again.
#
# WHY THE GUARD EXISTS (2026-06-22): the hardcoded fallback below previously held
# Zhu2013=1761, Giles2015=458, Dikow2009=1075 -- all WRONG and from different
# datasets/regimes (1761 is Conrad2008's gap-aware score; 458/1075 are *below*
# the most-permissive Fitch optimum so cannot be these datasets at all).  They
# nearly derailed the B1 read (the reach-gate keys off these targets); only a TNT
# ground-truth run recovered the true optima.  The whole rest of the benchmark
# suite (bench_iterate.R, bench_beam.R, diag_*) already used the correct values.
# Lesson: a literal target table is a silent-drift hazard -- derive it, don't
# retype it.  See dev/benchmarks/README-best-known-targets.md.
.bdBestKnownFallback <- c(Wortley2006 = 479, Eklund2004 = 440, Zanol2014 = 1261,
                          Zhu2013 = 624, Giles2015 = 670, Dikow2009 = 1606)
.BdBestKnown <- function() {
  csv <- file.path("dev", "benchmarks", "headtohead_phase0.csv")
  if (!file.exists(csv)) {
    # Return: hardcoded fallback (kept in sync with the CSV's `tnt` column)
    return(.bdBestKnownFallback)
  }
  h <- utils::read.csv(csv, stringsAsFactors = FALSE)
  tapply(h$tnt, h$dataset, min)   # best TNT Fitch score per dataset
}
.bestKnown <- .BdBestKnown()

# Inapplicable "-" -> missing "?" so the search runs standard Fitch (the
# apples-to-apples regime against TNT).  Mirrors dev/benchmarks/bench_iterate.R.
.FitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

BdLoadDataset <- function(name) {
  e <- new.env()
  utils::data("inapplicable.phyData", package = "TreeSearch", envir = e)
  .FitchPhy(e[["inapplicable.phyData"]][[name]])
}

# The recipe shared by every TS arm.  Only `fuseInterval` is ever varied
# (the sharing treatment); everything else — including poolSuboptimal — is held
# byte-identical, so the A/B isolates the fuse knob alone.
# poolSuboptimal > 0 is REQUIRED: pairwise fuse (TS_FUSE_PAIRWISE=1) needs
# genuinely suboptimal pool trees to use as recipients; with the default
# poolSuboptimal=0 the pool holds only tied-best trees and fuse has no recipient
# that a donor clade could strictly improve.  Held identical in BOTH arms.
.bdPoolSuboptimal <- as.double(Sys.getenv("BD_POOL_SUBOPT", "5"))
BdControl <- function(fuseInterval = 3L) {
  SearchControl(fuseInterval = as.integer(fuseInterval),
                poolSuboptimal = .bdPoolSuboptimal)
}

# --- Arm runners --------------------------------------------------------------

# One independent full-recipe restart.  maxReplicates = 1 => no inter-replicate
# fusing/consensus can fire, so this is a pure single-tree hill-climb.
BdRunIndep <- function(dataset, seed, control = BdControl()) {
  set.seed(seed)
  r <- MaximizeParsimony(dataset, maxReplicates = 1L, targetHits = 99999L,
                         maxSeconds = 0, nThreads = 1L, verbosity = 0L,
                         control = control)
  list(variant = "indep", seed = seed,
       score = as.double(attr(r, "score")),
       replicateScores = as.double(attr(r, "replicate_scores")),
       candidates = as.double(attr(r, "candidates_evaluated")),
       trees = .AsMultiPhylo(r))
}

# One production run, fusing on (PAIRWISE) or off.
#
# B1 EXPERIMENT (the decisive arm): `prod_on` enables PAIRWISE tree-fusing
# (TS_FUSE_PAIRWISE=1 — fuse into suboptimal pool recipients, the Goloboff-1999
# config where recombination can strictly improve a weaker tree), vs `prod_off`
# which never fuses.  Both share the SAME poolSuboptimal (BdControl), so the
# only difference is functional recombination on/off.
#
# TIME-MATCHED on WALL-CLOCK, not replicates.  Pairwise fuse does extra
# tree_fuse+score_tree work per interval, so a fixed-replicate budget would
# silently hand the fuse arm more compute (and candidate-count matching is
# wrong too — tree_fuse scores on a separate uncounted path).  When
# BD_PROD_SECONDS>0 each arm gets exactly that many seconds (reps becomes a high
# cap); fuse-ON wins only if recombination beats the restarts it crowds out.
BdRunProd <- function(dataset, seed, reps, fuse = TRUE) {
  # Toggle the prototype flag for THIS run.  Read per-fuse in the kernel, so
  # setting it here (before the call) scopes it correctly even in-process.
  Sys.setenv(TS_FUSE_PAIRWISE = if (fuse) "1" else "0")
  on.exit(Sys.unsetenv("TS_FUSE_PAIRWISE"), add = TRUE)

  ctrl <- if (fuse) BdControl(fuseInterval = 3L)
          else BdControl(fuseInterval = reps + 1L)   # never fires => fusing off

  budgetSec <- as.double(Sys.getenv("BD_PROD_SECONDS", "0"))
  maxReps <- if (budgetSec > 0) 9999L else as.integer(reps)

  set.seed(seed)
  r <- MaximizeParsimony(dataset, maxReplicates = maxReps,
                         targetHits = 99999L, maxSeconds = budgetSec,
                         nThreads = 1L, verbosity = 0L, control = ctrl)
  list(variant = if (fuse) "prod_on" else "prod_off", seed = seed,
       score = as.double(attr(r, "score")),
       repsDone = length(as.double(attr(r, "replicate_scores"))),
       replicateScores = as.double(attr(r, "replicate_scores")),
       candidates = as.double(attr(r, "candidates_evaluated")),
       trees = .AsMultiPhylo(r))
}

# TNT xmult retaining many MPTs.  Reported as a SEPARATE optimal-island axis
# (best score + count of distinct best-score trees within TNT's own set).
BdRunTnt <- function(dataset, name, seed,
                     tntExe = Sys.getenv("TNT_EXE",
                                "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"),
                     hits = 10L, replic = 10L) {
  if (!nzchar(Sys.which(tntExe)) && !file.exists(tntExe)) {
    message("  [TNT] executable not found; skipping TNT arm.")
    # Return: NULL
    return(NULL)
  }
  wd <- file.path(tempdir(), paste0("bd_tnt_", Sys.getpid(), "_", name, "_", seed))
  unlink(wd, recursive = TRUE)
  dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(wd, recursive = TRUE), add = TRUE)
  WriteTntCharacters(dataset, file.path(wd, "data.tnt"))
  lines <- c("mxram 1024;", "proc data.tnt;",
             sprintf("rseed %d;", seed), "hold 10000;",
             sprintf("xmult = hits %d replic %d;", hits, replic),
             "tsave *finalt.tre;", "save;", "tsave/;", "quit;")
  writeLines(lines, file.path(wd, "runme.run"))
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  invisible(suppressWarnings(
    system2(tntExe, args = "runme.run;", stdout = TRUE, stderr = TRUE)))
  setwd(old)
  tr <- tryCatch(ReadTntTree(file.path(wd, "finalt.tre")),
                 error = function(e) NULL)
  if (is.null(tr)) {
    message("  [TNT] no trees parsed; skipping.")
    # Return: NULL
    return(NULL)
  }
  tr <- .AsMultiPhylo(tr)
  sc <- vapply(tr, function(x) TreeLength(x, dataset), numeric(1))
  list(variant = "tnt", seed = seed, score = min(sc),
       trees = tr, treeScores = sc)
}

# --- Coercion + metrics -------------------------------------------------------

.AsMultiPhylo <- function(x) {
  if (inherits(x, "multiPhylo")) {
    # Return: drop attributes that confuse downstream list ops
    structure(unclass(x)[seq_along(x)], class = "multiPhylo")
  } else if (inherits(x, "phylo")) {
    structure(list(x), class = "multiPhylo")
  } else {
    structure(x, class = "multiPhylo")
  }
}

# Distinct-topology class id per tree (RF == 0 is an equivalence relation, so a
# single greedy pass suffices: every j with RF(i,j)==0 is identical to i).
BdDistinctClasses <- function(trees) {
  n <- length(trees)
  if (n <= 1L) {
    # Return:
    return(rep(1L, n))
  }
  rf <- as.matrix(RobinsonFoulds(trees))
  cls <- integer(n)
  k <- 0L
  for (i in seq_len(n)) {
    if (cls[i] == 0L) {
      k <- k + 1L
      cls[which(rf[i, ] == 0)] <- k
    }
  }
  # Return:
  cls
}

# Mean normalized pairwise clustering-information distance (diversity magnitude).
# ClusteringInfoDist, never Robinson-Foulds, for the magnitude (memory: RF is
# inflated by a single rogue tip; RF==0 above is only used for exact identity).
BdMeanCid <- function(trees) {
  if (length(trees) <= 1L) {
    # Return:
    return(NA_real_)
  }
  # Return:
  mean(as.dist(ClusteringInfoDist(trees, normalize = TRUE)))
}

# Rarefaction: expected number of distinct topologies after m restarts,
# averaged over nPerm random restart orderings.  Lightweight (secondary).
BdDiscoveryCurve <- function(classes, nPerm = 200L) {
  n <- length(classes)
  if (n == 0L) {
    # Return:
    return(numeric(0))
  }
  acc <- numeric(n)
  for (p in seq_len(nPerm)) {
    ord <- sample.int(n)
    seen <- logical(max(classes))
    cum <- 0L
    for (m in seq_len(n)) {
      cl <- classes[ord[m]]
      if (!seen[cl]) { seen[cl] <- TRUE; cum <- cum + 1L }
      acc[m] <- acc[m] + cum
    }
  }
  # Return:
  acc / nPerm
}

# Per-arm summary block (gate-first: score-reaching, then conditional diversity).
BdArmStats <- function(trees, scores, bestKnown, label) {
  scores <- scores[is.finite(scores)]
  best <- if (length(scores)) min(scores) else NA_real_
  classesAll <- BdDistinctClasses(trees)
  list(label = label,
       nRuns = length(scores),
       nTrees = length(trees),
       scoreMin = best,
       scoreMedian = if (length(scores)) stats::median(scores) else NA_real_,
       gap = best - bestKnown,
       fracOptimal = if (length(scores)) mean(scores <= bestKnown + 1e-6) else NA_real_,
       fracWithin1 = if (length(scores)) mean(scores <= bestKnown + 1) else NA_real_,
       nDistinct = length(unique(classesAll)),
       rediscoveryRate = if (length(trees)) 1 - length(unique(classesAll)) / length(trees) else NA_real_,
       meanCid = BdMeanCid(trees))
}

# --- The B1 verdict (encodes the score-gated decision tree above) -------------

BdInterpret <- function(indep, prodOn, prodOff, tnt, bestKnown, pairing = NULL) {
  v <- character(0)
  add <- function(...) v[[length(v) + 1L]] <<- sprintf(...)

  # Reliability-framed, NOT "ever reaches": fusing's payoff is frequency, so the
  # gate and the A/B both key off fracOptimal (per-arm reach rate) and the
  # paired per-seed comparison -- never the across-seed minimum.
  indepReaches  <- isTRUE(indep$fracOptimal > 0)
  indepReliable <- isTRUE(indep$fracOptimal >= 0.8)
  onReach  <- if (is.null(prodOn))  NA_real_ else prodOn$fracOptimal
  offReach <- if (is.null(prodOff)) NA_real_ else prodOff$fracOptimal
  tntReaches <- !is.null(tnt) && isTRUE(tnt$score <= bestKnown + 1e-6)
  diverse <- isTRUE(indep$meanCid > 0.1)   # coarse; report the number, don't over-trust the cut

  add("Best-known target = %g.  INDEP: min=%g, median=%g (gap %+g), %.0f%% of restarts reach optimum, %d distinct topologies / %d trees, mean CID=%.3f.",
      bestKnown, indep$scoreMin, indep$scoreMedian, indep$gap, 100 * indep$fracOptimal,
      indep$nDistinct, indep$nTrees, indep$meanCid)

  if (indepReliable) {
    add("VERDICT: Independent restarts already reach the optimum reliably (%.0f%%). A maintained diverse set would return the same answer -- B1 is NOT bounding quality on this config.",
        100 * indep$fracOptimal)
  } else if (is.null(prodOn) || is.null(prodOff)) {
    # No fuse A/B available: fall back to the INDEP-vs-TNT reach gate, but flag
    # it provisional -- recombination might still rescue the optimum.
    if (!indepReaches && tntReaches) {
      add("VERDICT (provisional -- no PROD arms): TS independent restarts NEVER reach the optimum TNT finds, which LOOKS like a per-restart REACH problem. Run the fuse A/B (prod_on vs prod_off) before concluding: recombination may still rescue the optimum.")
    } else {
      add("VERDICT: PROD arms missing; run the fuse A/B (BD_MODE=cell prod variants) to discriminate B1.")
    }
  } else {
    # Gap exists and the PAIRED fuse A/B is the PRIMARY discriminator -- it is
    # consulted even when no single restart reaches the optimum, because
    # recombining suboptimal restarts to the optimum is exactly B1's mechanism.
    if (!is.null(pairing)) {
      budgetSec <- as.double(Sys.getenv("BD_PROD_SECONDS", "0"))
      budgetStr <- if (budgetSec > 0) sprintf("%gs wall-clock", budgetSec) else "fixed-reps"
      add("Pairwise-fuse A/B (paired over %d seeds; prod_on = TS_FUSE_PAIRWISE recombination, both arms poolSuboptimal=%g, %s budget): fusing strictly helps on %.0f%% of seeds, hurts %.0f%%, median score delta %+.2f (>0 = fusing better). Reach optimum: prod_on %.0f%% of seeds vs prod_off %.0f%%.",
          pairing$nShared, .bdPoolSuboptimal, budgetStr,
          100 * pairing$helpsFrac, 100 * pairing$hurtsFrac,
          pairing$medianDelta, 100 * onReach, 100 * offReach)
    }
    reachGain <- onReach - offReach
    fuseMaterial <- isTRUE(reachGain >= 0.1) ||
      (!is.null(pairing) && isTRUE(pairing$helpsFrac >= 0.3))
    if (fuseMaterial) {
      add("VERDICT: Recombination (fusing) materially improves reach/score over pure restarts (reach +%.0f pts; helps %.0f%% of seeds) -- B1's mechanism is ACTIVE. A stronger diverse-set + fusing tour is worth building.",
          100 * reachGain, if (is.null(pairing)) NA else 100 * pairing$helpsFrac)
    } else if (isTRUE(onReach >= 0.8) && isTRUE(offReach >= 0.8)) {
      add("VERDICT: Both arms reliably reach the optimum at this budget (on %.0f%%, off %.0f%%) and fusing adds little -- a budget/effort effect, not an architecture limit. Re-test at a tighter budget if the production gap persists.",
          100 * onReach, 100 * offReach)
    } else if (isTRUE(onReach < 0.5) && isTRUE(offReach < 0.5) && tntReaches) {
      # Neither TS arm reliably reaches, TNT does, and fusing did NOT help.
      # INDEP diversity disambiguates -- HIGH diversity means fusing HAD varied
      # material and still failed (optimum is not a recombination); LOW
      # diversity means fusing was starved (diversity itself is the bottleneck).
      if (diverse) {
        add("VERDICT: Restarts already produce DIVERSE basins (mean CID=%.3f) so fusing had varied material, yet still cannot reach the optimum TNT finds -- the optimum is NOT a recombination of what restarts reach. This indicts per-restart STRENGTH (reach), not diversity. B1 is NOT the lever.",
            indep$meanCid)
      } else {
        add("VERDICT: Restarts collapse to FEW similar basins (mean CID=%.3f) so fusing is starved of varied material -- diversity itself is the bottleneck. B1 PLAUSIBLE: a maintained diverse-set tour would give fusing material; build it and re-test.",
            indep$meanCid)
      }
    } else {
      add("VERDICT: Inconclusive at this budget (on reach=%.0f%%, off reach=%.0f%%, tnt reaches=%s). Increase seeds/reps.",
          100 * onReach, 100 * offReach, tntReaches)
    }
  }
  if (!is.null(tnt)) {
    tntDistinct <- length(unique(BdDistinctClasses(tnt$trees[tnt$treeScores <= tnt$score + 1e-6])))
    add("TNT optimal island (separate axis): best=%g, %d distinct trees at best score.",
        tnt$score, tntDistinct)
  }
  # Return:
  unlist(v)
}

# --- Grid + cell dispatch (Hamilton array unit) -------------------------------

BdCellGrid <- function(indepK, prodReps, prodSeeds) {
  indep <- data.frame(variant = "indep", seed = seq_len(indepK),
                      reps = 1L, stringsAsFactors = FALSE)
  prod <- expand.grid(variant = c("prod_on", "prod_off"),
                      seed = seq_len(prodSeeds),
                      stringsAsFactors = FALSE)
  prod$reps <- prodReps
  prod <- prod[, c("variant", "seed", "reps")]
  # Return:
  rbind(indep, prod)
}

BdRunCell <- function(cell, dataset) {
  if (cell$variant == "indep") {
    BdRunIndep(dataset, cell$seed)
  } else {
    BdRunProd(dataset, cell$seed, cell$reps, fuse = (cell$variant == "prod_on"))
  }
}

# --- Analysis: merge partials -> stats + verdict + figure ---------------------

BdAnalyze <- function(partialDir, name, bestKnown, outPrefix, tnt = NULL,
                      nPerm = as.integer(Sys.getenv("BD_NPERM", "200"))) {
  files <- list.files(partialDir, pattern = "^cell_.*\\.rds$", full.names = TRUE)
  if (!length(files)) stop("No partial .rds files in ", partialDir)
  parts <- lapply(files, readRDS)
  byVar <- split(parts, vapply(parts, function(p) p$variant, character(1)))

  gather <- function(var) {
    ps <- byVar[[var]]
    if (is.null(ps)) {
      # Return:
      return(NULL)
    }
    trees <- .AsMultiPhylo(do.call(c, lapply(ps, function(p) p$trees)))
    scores <- vapply(ps, function(p) p$score, numeric(1))  # one best score per run
    seeds <- vapply(ps, function(p) as.double(p$seed), numeric(1))
    list(trees = trees, scores = scores, seeds = seeds)
  }
  gi <- gather("indep"); gon <- gather("prod_on"); goff <- gather("prod_off")

  indep <- if (!is.null(gi))
    BdArmStats(gi$trees, gi$scores, bestKnown, "INDEP") else NULL
  prodOn <- if (!is.null(gon))
    BdArmStats(gon$trees, gon$scores, bestKnown, "PROD_ON") else NULL
  prodOff <- if (!is.null(goff))
    BdArmStats(goff$trees, goff$scores, bestKnown, "PROD_OFF") else NULL

  # Paired per-seed fuse comparison (prod_on / prod_off share seeds 1..S).
  # Fusing's benefit is reliability/frequency, NOT whether the optimum is ever
  # touched, so we pair by seed rather than comparing the across-seed minima.
  pairing <- NULL
  if (!is.null(gon) && !is.null(goff)) {
    onBy <- tapply(gon$scores, gon$seeds, min)
    offBy <- tapply(goff$scores, goff$seeds, min)
    shared <- intersect(names(onBy), names(offBy))
    if (length(shared)) {
      o <- onBy[shared]; f <- offBy[shared]
      pairing <- list(nShared = length(shared),
                      helpsFrac = mean(o < f),       # fusing strictly better
                      hurtsFrac = mean(o > f),
                      medianDelta = stats::median(f - o))  # >0 = fusing better
    }
  }

  verdict <- BdInterpret(indep, prodOn, prodOff, tnt, bestKnown, pairing)

  # Durable summary table (one row per arm).
  rows <- Filter(Negate(is.null), list(indep, prodOn, prodOff))
  summary <- do.call(rbind, lapply(rows, function(s) data.frame(
    dataset = name, arm = s$label, nRuns = s$nRuns, nTrees = s$nTrees,
    scoreMin = s$scoreMin, scoreMedian = s$scoreMedian, gap = s$gap,
    fracOptimal = s$fracOptimal, fracWithin1 = s$fracWithin1,
    nDistinct = s$nDistinct, rediscoveryRate = s$rediscoveryRate,
    meanCid = s$meanCid, stringsAsFactors = FALSE)))
  utils::write.csv(summary, paste0(outPrefix, "_summary.csv"), row.names = FALSE)
  writeLines(verdict, paste0(outPrefix, "_verdict.txt"))

  # Figure: (1) INDEP score gate, (2) fuse A/B, (3) discovery curve.
  .BdFigure(paste0(outPrefix, "_figure.pdf"), name, bestKnown,
            gi, gon, goff, tnt, nPerm)

  cat("\n==== Basin-diversity verdict (", name, ") ====\n", sep = "")
  cat(verdict, sep = "\n")
  cat("\n\nSummary written to ", outPrefix, "_summary.csv\n", sep = "")
  # Return:
  invisible(list(summary = summary, verdict = verdict,
                 indep = indep, prodOn = prodOn, prodOff = prodOff, tnt = tnt))
}

.BdFigure <- function(path, name, bestKnown, gi, gon, goff, tnt, nPerm) {
  ok <- tryCatch({ grDevices::pdf(path, width = 11, height = 4); TRUE },
                 error = function(e) FALSE)
  if (!ok) {
    # Return:
    return(invisible(NULL))
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), cex = 0.9)

  # Panel 1: INDEP per-restart converged scores (THE gate).
  if (!is.null(gi)) {
    s <- gi$scores
    graphics::plot(jitter(rep(1, length(s))), s, xlim = c(0.5, 1.5),
                   ylim = range(c(s, bestKnown)), xaxt = "n", xlab = "",
                   ylab = "converged score", pch = 19,
                   col = grDevices::adjustcolor("black", 0.5),
                   main = sprintf("INDEP gate: %.0f%% reach optimum",
                                  100 * mean(s <= bestKnown + 1e-6)))
    graphics::abline(h = bestKnown, col = "firebrick", lwd = 2)
    graphics::axis(1, at = 1, labels = sprintf("n=%d restarts", length(s)))
  }

  # Panel 2: fuse A/B (the decisive probe) -- per-seed, PAIRED (never the min).
  if (!is.null(gon) && !is.null(goff)) {
    onBy <- tapply(gon$scores, gon$seeds, min)
    offBy <- tapply(goff$scores, goff$seeds, min)
    shared <- intersect(names(onBy), names(offBy))
    allv <- c(offBy[shared], onBy[shared], bestKnown,
              if (!is.null(tnt)) tnt$score)
    graphics::plot(NA, xlim = c(0.7, 2.4), ylim = range(allv, na.rm = TRUE),
                   xaxt = "n", xlab = "", ylab = "converged score",
                   main = "Fuse A/B (paired per seed)")
    graphics::axis(1, at = c(1, 2), labels = c("prod_off", "prod_on"))
    for (s in shared)
      graphics::segments(1, offBy[s], 2, onBy[s],
                         col = grDevices::adjustcolor("grey50", 0.6))
    graphics::points(rep(1, length(shared)), offBy[shared], pch = 19, col = "steelblue")
    graphics::points(rep(2, length(shared)), onBy[shared], pch = 19, col = "steelblue")
    graphics::abline(h = bestKnown, col = "firebrick", lwd = 2)
    if (!is.null(tnt))
      graphics::points(2.35, tnt$score, pch = 4, col = "darkgreen", cex = 1.4, lwd = 2)
  }

  # Panel 3: discovery (rarefaction) curve for INDEP (secondary).
  if (!is.null(gi)) {
    cls <- BdDistinctClasses(gi$trees)
    dc <- BdDiscoveryCurve(cls, nPerm)
    graphics::plot(seq_along(dc), dc, type = "l", lwd = 2, col = "black",
                   xlab = "restarts", ylab = "distinct topologies",
                   main = "INDEP discovery curve")
    graphics::abline(h = length(unique(cls)), col = "grey60", lty = 2)
  }
}

# --- Driver -------------------------------------------------------------------

.bdMain <- function() {
  mode <- Sys.getenv("BD_MODE", "smoke")
  name <- Sys.getenv("BD_DATASET", "Zanol2014")
  bestKnown <- as.double(Sys.getenv("BD_BEST",
                  as.character(.bestKnown[[name]] %||% NA)))
  partialDir <- Sys.getenv("BD_PARTIAL_DIR",
                  file.path("dev", "benchmarks", "basin_partials"))
  outPrefix <- Sys.getenv("BD_OUT_PREFIX",
                  file.path("dev", "benchmarks", paste0("basin_", name)))
  wantTnt <- Sys.getenv("BD_TNT", "1") == "1"

  if (mode == "cell") {
    args <- commandArgs(trailingOnly = TRUE)
    idx <- as.integer(if (length(args)) args[[1]]
                      else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))   # 0-based
    indepK <- as.integer(Sys.getenv("BD_INDEP_K", "50"))
    prodReps <- as.integer(Sys.getenv("BD_PROD_REPS", "48"))
    prodSeeds <- as.integer(Sys.getenv("BD_PROD_SEEDS", "8"))
    grid <- BdCellGrid(indepK, prodReps, prodSeeds)
    cell <- grid[idx + 1L, ]
    dataset <- BdLoadDataset(name)
    res <- BdRunCell(cell, dataset)
    dir.create(partialDir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, file.path(partialDir, sprintf("cell_%04d.rds", idx)))
    cat(sprintf("cell %d  %s seed=%d  score=%g\n",
                idx, cell$variant, cell$seed, res$score))

  } else if (mode == "analyze") {
    dataset <- BdLoadDataset(name)
    tnt <- if (wantTnt)
      BdRunTnt(dataset, name, seed = 1L,
               hits = as.integer(Sys.getenv("BD_TNT_HITS", "10")),
               replic = as.integer(Sys.getenv("BD_TNT_REPLIC", "10"))) else NULL
    BdAnalyze(partialDir, name, bestKnown, outPrefix, tnt = tnt)

  } else {  # smoke: small end-to-end run, local; plumbing + gap sanity check.
    indepK <- as.integer(Sys.getenv("BD_INDEP_K", "6"))
    prodReps <- as.integer(Sys.getenv("BD_PROD_REPS", "8"))
    prodSeeds <- as.integer(Sys.getenv("BD_PROD_SEEDS", "2"))
    dataset <- BdLoadDataset(name)
    cat(sprintf("[smoke] %s  K=%d  prodReps=%d  prodSeeds=%d\n",
                name, indepK, prodReps, prodSeeds))
    dir.create(partialDir, recursive = TRUE, showWarnings = FALSE)
    unlink(list.files(partialDir, pattern = "^cell_.*\\.rds$", full.names = TRUE))
    grid <- BdCellGrid(indepK, prodReps, prodSeeds)
    for (i in seq_len(nrow(grid))) {
      res <- BdRunCell(grid[i, ], dataset)
      saveRDS(res, file.path(partialDir, sprintf("cell_%04d.rds", i - 1L)))
      cat(sprintf("  cell %d/%d  %s seed=%d  score=%g\n",
                  i, nrow(grid), grid$variant[i], grid$seed[i], res$score))
    }
    tnt <- if (wantTnt)
      BdRunTnt(dataset, name, seed = 1L, hits = 3L, replic = 3L) else NULL
    BdAnalyze(partialDir, name, bestKnown, outPrefix, tnt = tnt)
  }
}

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

if (sys.nframe() == 0L) .bdMain()
