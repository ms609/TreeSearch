#!/usr/bin/env Rscript
# A/B: adaptive-start bandit credit for user-warm-started replicates.
#
#   before = bandit credits WAGNER_RANDOM for every `tree = ` warm-started rep
#   after  = warm-started reps are excluded from the bandit (decay unchanged)
#
# Only the `adaptiveStart = TRUE` + `tree = <multiPhylo>` combination can
# differ: with no supplied tree every rep builds its own start, the guard is
# always true, and the two builds are identical.  The change bites only when
# spare COLD replicates remain after the pool is consumed, because those are
# the reps whose arm is drawn from the (previously polluted) posterior.
# Hence nRep >> nPool below; with nRep == nPool there is nothing to measure.
#
# Pre-registered decision rule (fixed before looking at any output):
#   ADOPT unless time-to-optimum regresses — i.e. unless `after` is worse on
#   wall-clock at equal score, or worse on score.  A tie adopts, because the
#   phantom credit is wrong on its own terms.
#
# Score/wall at pilot scale is underpowered (integer scores, frequent ties).
# The high-signal readout is the arm-selection distribution over the cold
# replicates, which directly shows whether the pollution is gone.
#
# Usage:
#   Rscript dev/benchmarks/ab_bandit_warmstart.R [nSeed] [outCsv]
#
# Requires two installed libraries (see the header of the report in
# dev/benchmarks/ for how they were built):
#   .ab-bandit-before  .ab-bandit-after

args    <- commandArgs(trailingOnly = TRUE)
nSeed   <- if (length(args) >= 1L) as.integer(args[[1]]) else 5L
outCsv  <- if (length(args) >= 2L) args[[2]] else
  "dev/benchmarks/ab_bandit_warmstart.csv"

beforeLib <- ".ab-bandit-before"
afterLib  <- ".ab-bandit-after"
datasets  <- c("Wortley2006", "Griswold1999", "Agnarsson2004")
nPool     <- 5L    # trees supplied via tree =
nRep      <- 20L   # total replicates; nRep - nPool are cold (the ones that matter)
poolDir   <- file.path(tempdir(), "abBanditPools")
dir.create(poolDir, showWarnings = FALSE, recursive = TRUE)

# Phase 1: build one starting pool per (dataset, seed), shared verbatim by both
# arms so no phase-1 noise leaks into the comparison.  Generated with the
# bandit off, so the arm being tested cannot influence its own input.
makePool <- function(dsName, seed, path) {
  if (file.exists(path)) return(invisible(path))
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', afterLib),
    'suppressMessages(library(TreeSearch))',
    sprintf('ds <- inapplicable.phyData[["%s"]]', dsName),
    sprintf('set.seed(%d)', seed),
    sprintf('p <- MaximizeParsimony(ds, maxReplicates = %dL, targetHits = 99L,',
            nPool),
    '  adaptiveStart = FALSE, verbosity = 0L, nThreads = 1L)',
    sprintf('saveRDS(p[seq_len(min(%dL, length(p)))], "%s")',
            nPool, gsub("\\\\", "/", path))
  ), tmp)
  system2("Rscript", c("--no-save", tmp), stdout = FALSE, stderr = FALSE)
  unlink(tmp)
  invisible(path)
}

# Phase 2: resume from that pool with the bandit on, under each build.
runArm <- function(lib, label, dsName, seed, poolPath) {
  tmp <- tempfile(fileext = ".R")
  writeLines(c(
    sprintf('.libPaths(c("%s", .libPaths()))', lib),
    'suppressMessages(library(TreeSearch))',
    sprintf('ds <- inapplicable.phyData[["%s"]]', dsName),
    sprintf('pool <- readRDS("%s")', gsub("\\\\", "/", poolPath)),
    sprintf('set.seed(%d)', seed),
    't0 <- proc.time()',
    sprintf('res <- MaximizeParsimony(ds, tree = pool, maxReplicates = %dL,',
            nRep),
    '  targetHits = 99L, adaptiveStart = TRUE, verbosity = 0L, nThreads = 1L)',
    'wall <- (proc.time() - t0)[[3]]',
    'att <- attr(res, "strategy_diagnostics")$attempts',
    sprintf(
      'cat(sprintf("%s|%s|%d|%%g|%%d|%%.3f|%%s\\n",', label, dsName, seed),
    '  attr(res, "score"), attr(res, "replicates"), wall,',
    '  paste(att, collapse = ",")))'
  ), tmp)
  out <- system2("Rscript", c("--no-save", tmp), stdout = TRUE, stderr = FALSE)
  unlink(tmp)
  trimws(tail(out, 1))
}

rows <- list()
for (dsName in datasets) {
  cat(sprintf("\n=== %s ===\n", dsName))
  for (seed in seq_len(nSeed)) {
    poolPath <- file.path(poolDir, sprintf("%s_%d.rds", dsName, seed))
    makePool(dsName, seed, poolPath)
    if (!file.exists(poolPath)) {
      cat(sprintf("  seed %d: pool generation FAILED, skipping\n", seed))
      next
    }
    for (arm in list(c("before", beforeLib), c("after", afterLib))) {
      line  <- runArm(arm[[2]], arm[[1]], dsName, seed, poolPath)
      parts <- strsplit(line, "\\|")[[1]]
      if (length(parts) != 7L) {
        cat(sprintf("  seed %d %s: UNPARSEABLE (%s)\n", seed, arm[[1]], line))
        next
      }
      cat(" ", line, "\n")
      rows[[length(rows) + 1L]] <- data.frame(
        arm = parts[1], dataset = parts[2], seed = as.integer(parts[3]),
        score = as.numeric(parts[4]), reps = as.integer(parts[5]),
        wall = as.numeric(parts[6]), attempts = parts[7],
        stringsAsFactors = FALSE
      )
    }
  }
}

results <- do.call(rbind, rows)
write.csv(results, outCsv, row.names = FALSE)
cat(sprintf("\nWrote %s (%d rows)\n", outCsv, nrow(results)))

# ---- Summary -----------------------------------------------------------
if (!is.null(results) && nrow(results) > 0L) {
  cat("\n== Score / wall by arm ==\n")
  print(aggregate(cbind(score, wall) ~ arm + dataset, results, mean))

  # Paired comparison: same (dataset, seed) under both builds.
  wide <- merge(
    results[results$arm == "before", c("dataset", "seed", "score", "wall")],
    results[results$arm == "after",  c("dataset", "seed", "score", "wall")],
    by = c("dataset", "seed"), suffixes = c("Before", "After")
  )
  if (nrow(wide) > 0L) {
    cat("\n== Paired (after - before) ==\n")
    cat(sprintf("  score: mean %+.3f | after better %d, worse %d, tied %d\n",
                mean(wide$scoreAfter - wide$scoreBefore),
                sum(wide$scoreAfter < wide$scoreBefore),
                sum(wide$scoreAfter > wide$scoreBefore),
                sum(wide$scoreAfter == wide$scoreBefore)))
    cat(sprintf("  wall : mean %+.3f s | after faster %d, slower %d\n",
                mean(wide$wallAfter - wide$wallBefore),
                sum(wide$wallAfter < wide$wallBefore),
                sum(wide$wallAfter > wide$wallBefore)))
  }

  # Mechanism: total attempts per arm. `before` should sum to reps (every
  # warm rep votes); `after` should sum to reps - nPool.
  cat("\n== Bandit attempts by strategy (summed over runs) ==\n")
  attMat <- function(sub) {
    m <- do.call(rbind, lapply(strsplit(sub$attempts, ","),
                               function(x) as.numeric(x)))
    colSums(m)
  }
  for (a in c("before", "after")) {
    sub <- results[results$arm == a, ]
    if (nrow(sub) == 0L) next
    cs <- attMat(sub)
    cat(sprintf("  %-6s total %5.0f  per-arm %s\n", a, sum(cs),
                paste(sprintf("%.0f", cs), collapse = " ")))
  }
  cat("\n  (arm order: Wagner(random) Wagner(Goloboff) Wagner(entropy)",
      "Random tree)\n")
}
