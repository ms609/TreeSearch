# Fuse EFFICIENCY experiment: does enabling pairwise fuse + pool reach the known
# optimum SOONER (less wall-clock) than plain restarts?  This is the question the
# basin A/B did NOT answer -- that one compared FINAL score at a fixed budget and
# found both arms saturate the optimum.  Here we record the full anytime
# trajectory (elapsed, best-so-far) per seed per arm and compare time-to-optimum.
#
# Design notes (advisor-reviewed 2026-06-22):
# * WALL-CLOCK is the only fair currency (fuse replicates cost more; fuse scoring
#   is uncounted) -- we read `elapsed` (cumulative seconds) from progressCallback,
#   confirmed monotonic-from-start.
# * The anytime CURVE (median best-score at each wall-clock slice, paired by seed)
#   is censoring-free and shows magnitude; time-to-optimum is just where it
#   crosses score==target.
# * Focus is DIKOW2009 (88t): the only roster dataset whose optimum is NOT found
#   in the first restart or two, so the only one where fuse (which can't fire
#   until the pool holds suboptimal recipients) has room to help.  Easy datasets
#   are structurally blind to the effect.
# * Both arms hold poolSuboptimal>0 identical; only TS_FUSE_PAIRWISE differs.  The
#   PRODUCTION baseline is "no fuse, poolSuboptimal=0", so the honest framing is
#   "does enabling the fuse+pool machinery reach optima sooner than plain
#   restarts" -- which is exactly the user's question.
# * nThreads=1 (wall-clock + pool dynamics are nondeterministic in parallel).

suppressWarnings(suppressMessages({
  library(TreeSearch,
          lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-fuse"),
                                  winslash = "/", mustWork = FALSE))
  library(TreeTools)
}))

`%||%` <- function(a, b) if (is.null(a)) b else a

.feFitch <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
FeLoad <- function(name) {
  e <- new.env()
  utils::data("inapplicable.phyData", package = "TreeSearch", envir = e)
  .feFitch(e[["inapplicable.phyData"]][[name]])
}

# One timed run, recording the (elapsed, best) trajectory via the progress
# callback.  fuse = TRUE enables pairwise fuse; both arms share poolSuboptimal.
FeRun <- function(dataset, seed, fuse, capSeconds, poolSubopt) {
  Sys.setenv(TS_FUSE_PAIRWISE = if (fuse) "1" else "0")
  on.exit(Sys.unsetenv("TS_FUSE_PAIRWISE"), add = TRUE)
  ctrl <- SearchControl(
    fuseInterval = if (fuse) 3L else .Machine$integer.max %/% 2L,
    poolSuboptimal = poolSubopt)

  rows <- list()
  cb <- function(info) {
    rows[[length(rows) + 1L]] <<- list(
      elapsed = as.double(info$elapsed %||% NA),
      best    = as.double(info$best_score %||% NA),
      phase   = as.character(info$phase %||% NA))
  }
  set.seed(seed)
  invisible(MaximizeParsimony(
    dataset, maxReplicates = 99999L, targetHits = 99999L,
    maxSeconds = capSeconds, nThreads = 1L, verbosity = 0L,
    progressCallback = cb, control = ctrl))

  tr <- do.call(rbind.data.frame, c(rows, list(stringsAsFactors = FALSE)))
  tr <- tr[!is.na(tr$elapsed) & !is.na(tr$best), , drop = FALSE]
  tr <- tr[order(tr$elapsed), , drop = FALSE]
  tr$bestSoFar <- cummin(tr$best)        # the anytime best-so-far
  list(variant = if (fuse) "on" else "off", seed = seed,
       capSeconds = capSeconds, poolSubopt = poolSubopt,
       elapsed = tr$elapsed, bestSoFar = tr$bestSoFar, phase = tr$phase,
       finalBest = min(tr$bestSoFar))
}

# Time at which best-so-far first reaches `target` (NA = never within cap).
FeTimeToTarget <- function(run, target) {
  hit <- which(run$bestSoFar <= target + 1e-6)
  if (!length(hit)) return(NA_real_)
  run$elapsed[hit[1]]
}

# Best-so-far at a given wall-clock time t (step function; the value carried
# forward from the last callback at or before t; Inf before the first).
FeBestAt <- function(run, t) {
  ix <- which(run$elapsed <= t)
  if (!length(ix)) return(Inf)
  run$bestSoFar[ix[length(ix)]]
}

# --- Grid + cell dispatch (SLURM array unit): 2 arms x seeds -----------------
FeGrid <- function(seeds) {
  expand.grid(variant = c("on", "off"), seed = seq_len(seeds),
              stringsAsFactors = FALSE)
}

.feMain <- function() {
  mode  <- Sys.getenv("FE_MODE", "smoke")
  name  <- Sys.getenv("FE_DATASET", "Dikow2009")
  cap   <- as.double(Sys.getenv("FE_CAP_SECONDS", "60"))
  subo  <- as.double(Sys.getenv("FE_POOL_SUBOPT", "5"))
  seeds <- as.integer(Sys.getenv("FE_SEEDS", "48"))
  pdir  <- Sys.getenv("FE_PARTIAL_DIR",
                      file.path("dev", "benchmarks", paste0("fe_partials_", name)))

  if (mode == "cell") {
    args <- commandArgs(trailingOnly = TRUE)
    idx <- as.integer(if (length(args)) args[[1]]
                      else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
    grid <- FeGrid(seeds)
    cell <- grid[idx + 1L, ]
    d <- FeLoad(name)
    res <- FeRun(d, cell$seed, fuse = (cell$variant == "on"), cap, subo)
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(res, file.path(pdir, sprintf("fe_%04d.rds", idx)))
    cat(sprintf("[cell %d] %s seed=%d final=%g\n", idx, res$variant,
                res$seed, res$finalBest))
    return(invisible())
  }

  if (mode == "smoke") {
    d <- FeLoad(name)
    on1 <- FeRun(d, 1L, TRUE, cap, subo)
    of1 <- FeRun(d, 1L, FALSE, cap, subo)
    cat(sprintf("smoke %s seed1: on final=%g off final=%g\n",
                name, on1$finalBest, of1$finalBest))
    return(invisible())
  }

  stop("FE_MODE analyze is run via fe_analyze() with a target; see comment.")
}

# --- Analysis (run locally after pulling partials) ----------------------------
# target = the known optimum (headtohead_phase0.csv `tnt`).  Produces: time-to-
# optimum paired stats + an anytime survival/score curve.
FeAnalyze <- function(partialDir, target, slices = NULL) {
  fs <- list.files(partialDir, pattern = "fe_", full.names = TRUE)
  parts <- lapply(fs, readRDS)
  vv <- vapply(parts, function(p) p$variant, character(1))
  on  <- parts[vv == "on"];  off <- parts[vv == "off"]
  keyOn  <- vapply(on,  function(p) p$seed, numeric(1))
  keyOff <- vapply(off, function(p) p$seed, numeric(1))
  shared <- intersect(keyOn, keyOff)
  on  <- on[match(shared, keyOn)]
  off <- off[match(shared, keyOff)]

  tOn  <- vapply(on,  FeTimeToTarget, numeric(1), target = target)
  tOff <- vapply(off, FeTimeToTarget, numeric(1), target = target)
  cap  <- max(vapply(parts, function(p) max(p$elapsed), numeric(1)))
  if (is.null(slices)) slices <- seq(cap / 12, cap, length.out = 12)

  # Reach-rate (survival) and median best-so-far at each slice, paired.
  curve <- lapply(slices, function(t) {
    bOn  <- vapply(on,  FeBestAt, numeric(1), t = t)
    bOff <- vapply(off, FeBestAt, numeric(1), t = t)
    data.frame(t = t,
               reachOn  = mean(bOn  <= target + 1e-6),
               reachOff = mean(bOff <= target + 1e-6),
               medBestOn  = stats::median(bOn[is.finite(bOn)]),
               medBestOff = stats::median(bOff[is.finite(bOff)]))
  })
  curve <- do.call(rbind, curve)

  # Paired time-to-optimum: censor NA at cap for a conservative comparison.
  cOn  <- ifelse(is.na(tOn),  cap, tOn)
  cOff <- ifelse(is.na(tOff), cap, tOff)
  list(nShared = length(shared), target = target, cap = cap,
       medTimeOn = stats::median(cOn), medTimeOff = stats::median(cOff),
       reachOnFinal = mean(!is.na(tOn)), reachOffFinal = mean(!is.na(tOff)),
       fasterFrac = mean(cOn < cOff), slowerFrac = mean(cOn > cOff),
       medSpeedup = stats::median(cOff - cOn),  # >0 => fuse reaches sooner
       curve = curve)
}

if (sys.nframe() == 0L && !nzchar(Sys.getenv("FE_NO_MAIN"))) .feMain()
