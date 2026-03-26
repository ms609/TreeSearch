#!/usr/bin/env Rscript
# Quick IW (XPIWE) benchmark — TreeSearch vs TNT 1.6
# Uses extended implied weights (k=10, r=0.5, max_f=5) on both engines.
# Skips the largest datasets; focuses on small-to-medium for quick results.
#
# Usage: Rscript bench_iw_quick.R
# Requires: TreeSearch.F built in .builds/TreeSearch-F/

cat("=== IW (XPIWE, k=10) benchmark ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

AGENT_LIB <- normalizePath("../../../.builds/TreeSearch-F", mustWork = TRUE)
.libPaths(c(AGENT_LIB, .libPaths()))

library(TreeSearch.F)
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))
library(TreeTools)

TNT_EXE <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
stopifnot("TNT not found" = file.exists(TNT_EXE))

STAGING_DIR <- normalizePath("../../.tnt-bench", mustWork = FALSE)
dir.create(STAGING_DIR, showWarnings = FALSE)

OUTFILE <- sprintf("iw_xpiwe_k10_%s.csv", format(Sys.time(), "%Y%m%d_%H%M"))
CONCAVITY <- 10

# --- Helpers ---
prepare_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(dataset),
    n_chars = sum(at$weight)
  )
}

obs_count <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  levels <- at$levels
  is_missing <- apply(contrast, 1, function(row) all(row == 1))
  inapp_col <- match("-", levels)
  if (!is.na(inapp_col)) {
    is_inapp <- apply(contrast, 1, function(row) {
      row[inapp_col] == 1 && sum(row) == 1
    })
    is_missing <- is_missing | is_inapp
  }
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  vapply(seq_len(ncol(tip_data)), function(p) {
    sum(!is_missing[tip_data[, p]])
  }, integer(1))
}

min_steps_for <- function(dataset) {
  MinimumLength(dataset, compress = TRUE)
}

run_tnt <- function(data_file, timeout_s, seed, hits, reps, concavity) {
  hh <- timeout_s %/% 3600
  mm <- (timeout_s %% 3600) %/% 60
  ss <- timeout_s %% 60
  # piwe must be set BEFORE proc so data is read with IW on.
  # xpiwe( activates per-character concavity (Goloboff 2014 XPIWE).
  # Do NOT use xpiwe(...) — the ) is a separate sub-option that reverts to
  # single concavity, silently undoing the correction.
  commands <- c(
    "mxram 1024;",
    sprintf("piwe=%d;", concavity),
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed),
    "xpiwe(;",
    sprintf("timeout %d:%02d:%02d;", hh, mm, ss),
    sprintf("xmult=hits %d replic %d;", hits, reps),
    "best;", "quit;"
  )
  script_path <- file.path(STAGING_DIR, "tntbench.run")
  writeLines(commands, script_path)

  old_wd <- setwd(STAGING_DIR)
  on.exit(setwd(old_wd), add = TRUE)

  t0 <- proc.time()
  output <- withCallingHandlers(
    system2(TNT_EXE, args = "tntbench.run;",
            stdout = TRUE, stderr = TRUE, timeout = timeout_s + 60),
    warning = function(w) invokeRestart("muffleWarning")
  )
  wall_s <- as.double((proc.time() - t0)[3])

  output <- iconv(output, from = "", to = "UTF-8", sub = "")
  out_text <- paste(output, collapse = "\n")

  score <- NA_real_
  # [0-9]+[.][0-9]+ avoids capturing the trailing period in "3.80000."
  m <- regmatches(out_text, regexpr("Best score:\\s+([0-9]+[.][0-9]+)", out_text))
  if (length(m) == 1) score <- as.numeric(sub("Best score:\\s+", "", m))

  n_trees <- NA_integer_
  m <- regmatches(out_text, regexpr("([0-9]+) trees? retained", out_text))
  if (length(m) == 1) n_trees <- as.integer(sub(" trees? retained", "", m))

  list(score = score, n_trees = n_trees, wall_s = wall_s)
}

run_treesearch <- function(ds, timeout_s, seed, hits, reps, concavity,
                           min_steps, obs_count) {
  set.seed(seed)
  t0 <- proc.time()
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = as.integer(reps),
    targetHits = as.integer(hits),
    maxSeconds = as.double(timeout_s),
    verbosity = 0L, nThreads = 1L,
    concavity = as.double(concavity),
    min_steps = min_steps,
    xpiwe = TRUE,
    xpiwe_r = 0.5,
    xpiwe_max_f = 5.0,
    obs_count = obs_count,
    # Tuned "default" strategy
    ratchetCycles = 12L,
    ratchetPerturbProb = 0.25,
    driftCycles = 2L,
    nniFirst = TRUE,
    outerCycles = 1L,
    maxOuterResets = 2L,
    adaptiveLevel = TRUE
  )
  wall_s <- as.double((proc.time() - t0)[3])

  list(score = result$best_score, n_trees = result$pool_size,
       wall_s = wall_s, replicates = result$replicates,
       hits = result$hits_to_best)
}

# --- Parameters ---
# Quick set: skip the very largest (Dikow 88t, Zhu 75t, Giles 78t, Zanol 74t)
BENCHMARK_NAMES <- c(
  "Longrich2010", "Vinther2008", "Sansom2010", "DeAssis2011",
  "Aria2015", "Wortley2006", "Griswold1999", "Schulze2007",
  "Eklund2004", "Agnarsson2004",
  # Medium-hard from round 2
  "Rougier2012", "Shultz2007", "Wilson2003", "Wetterer2000",
  "Capa2011", "Conrad2008", "OMeara2014", "Geisler2001",
  "Liljeblad2008", "Aguado2009"
)
TIMEOUT <- 30
SEEDS <- 1:3
HITS <- 5L
REPS <- 20L

# --- Prepare datasets ---
# IW uses original data (with inapplicable tokens), not cleaned.
# Both XPIWE implementations handle inapplicable/missing in their
# obs_count calculation.
cat("Preparing", length(BENCHMARK_NAMES), "datasets (original, not cleaned)...\n")
datasets <- list()
ds_min_steps <- list()
ds_obs_count <- list()
for (nm in BENCHMARK_NAMES) {
  raw <- inapplicable.phyData[[nm]]
  datasets[[nm]] <- prepare_ts_data(raw)
  ds_min_steps[[nm]] <- min_steps_for(raw)
  ds_obs_count[[nm]] <- obs_count(raw)
  WriteTntCharacters(raw,
                     filepath = file.path(STAGING_DIR, paste0(nm, ".tnt")))
}

# --- Run benchmark ---
rows <- list()
idx <- 0L
n_total <- length(BENCHMARK_NAMES) * length(SEEDS)

cat(sprintf("Running %d combinations (k=%d, %ds timeout)\n\n",
            n_total, CONCAVITY, TIMEOUT))

for (nm in BENCHMARK_NAMES) {
  ds <- datasets[[nm]]
  data_file <- paste0(nm, ".tnt")

  for (seed in SEEDS) {
    idx <- idx + 1L
    cat(sprintf("[%d/%d] %s (seed=%d, %dt) ... ",
                idx, n_total, nm, seed, ds$n_taxa))

    tnt <- run_tnt(data_file, timeout_s = TIMEOUT, seed = seed,
                   hits = HITS, reps = REPS, concavity = CONCAVITY)
    ts <- run_treesearch(ds, timeout_s = TIMEOUT, seed = seed,
                         hits = HITS, reps = REPS, concavity = CONCAVITY,
                         min_steps = ds_min_steps[[nm]],
                         obs_count = ds_obs_count[[nm]])

    rows[[idx]] <- data.frame(
      dataset = nm, n_taxa = ds$n_taxa, n_chars = ds$n_chars,
      seed = seed, timeout_s = TIMEOUT, concavity = CONCAVITY,
      tnt_score = tnt$score, tnt_trees = tnt$n_trees,
      tnt_wall_s = round(tnt$wall_s, 3),
      ts_score = ts$score, ts_trees = ts$n_trees,
      ts_wall_s = round(ts$wall_s, 3),
      ts_reps = ts$replicates, ts_hits = ts$hits,
      stringsAsFactors = FALSE
    )

    cat(sprintf("TNT=%.5g (%.1fs)  TS=%.5g (%.1fs)\n",
                tnt$score, tnt$wall_s, ts$score, ts$wall_s))

    results <- do.call(rbind, rows)
    write.csv(results, OUTFILE, row.names = FALSE)
  }
}

cat("\n=== IW benchmark complete ===\n")
cat("Results saved to:", OUTFILE, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

results <- do.call(rbind, rows)
cat("\nScore comparison (best across seeds):\n")
agg <- aggregate(cbind(tnt_score, ts_score) ~ dataset + n_taxa, data = results, FUN = min)
agg$delta <- agg$ts_score - agg$tnt_score
agg <- agg[order(agg$n_taxa), ]
print(agg, row.names = FALSE)
cat("\nTS wins:", sum(agg$delta < 0),
    " Ties:", sum(agg$delta == 0),
    " TNT wins:", sum(agg$delta > 0), "\n")
