#!/usr/bin/env Rscript
# T-290b: Phase profiling across scoring modes and weightings
#
# Runs 6 representative Brazeau-sample datasets under all 4 scoring
# conditions ({Brazeau, Fitch} x {EW, IW k=10}) with both strategy
# presets, capturing per-phase timing and per-replicate scores.
#
# Usage:
#   Rscript t290b_phase_profile.R <output_dir>

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "."

repo_root <- getwd()
if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error = function(e) getwd()
  )
  repo_root <- normalizePath(file.path(script_dir, "..", ".."),
                             mustWork = FALSE)
}
setwd(repo_root)

cat("=== T-290b: Phase Profiling ===\n")
cat("Repo root:", repo_root, "\n")
cat("Output dir:", outdir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")

# ---- Configuration ----
SUBSET <- c("project2346", "project2359", "project4146_(3)",
            "project2086", "project2771", "project804")

TIMING_NAMES <- c("wagner_ms", "tbr_ms", "xss_ms", "rss_ms", "css_ms",
                   "ratchet_ms", "drift_ms", "nni_perturb_ms", "anneal_ms",
                   "prune_reinsert_ms", "final_tbr_ms", "fuse_ms")

TIMEOUT <- 30
N_SEEDS <- 3L
STRATEGIES <- c("default", "thorough")
WEIGHTINGS <- list(EW = Inf, IW10 = 10)
SCORING <- c("brazeau", "fitch")

# ---- Load datasets ----
cat <- load_mbank_catalogue()

raw_datasets <- list()
for (k in SUBSET) {
  row <- cat[cat$key == k, ]
  nex_path <- file.path(NEOTRANS_MATRICES_DIR, row$filename)
  raw_datasets[[k]] <- suppressWarnings(TreeTools::ReadAsPhyDat(nex_path))
}
cat(sprintf("Loaded %d datasets\n\n", length(raw_datasets)))

# ---- Run profiling grid ----
conditions <- expand.grid(
  dataset = SUBSET,
  scoring = SCORING,
  weighting = names(WEIGHTINGS),
  strategy = STRATEGIES,
  seed = seq_len(N_SEEDS),
  stringsAsFactors = FALSE
)
n_total <- nrow(conditions)
cat(sprintf("Phase profiling grid: %d runs\n", n_total))
cat(sprintf("  6 datasets x 2 scoring x 2 weightings x 2 strategies x %d seeds\n\n",
            N_SEEDS))

rows <- vector("list", n_total)

for (i in seq_len(n_total)) {
  cond <- conditions[i, ]
  ds_raw <- raw_datasets[[cond$dataset]]
  ds_info <- cat[cat$key == cond$dataset, ]

  # Apply fitch_mode if needed
  ds <- if (cond$scoring == "fitch") fitch_mode(ds_raw) else ds_raw
  conc <- WEIGHTINGS[[cond$weighting]]

  cat(sprintf("[%d/%d] %s %s %s %s seed=%d ... ",
              i, n_total, cond$dataset, cond$scoring,
              cond$weighting, cond$strategy, cond$seed))

  set.seed(cond$seed)
  t0 <- proc.time()
  result <- tryCatch(
    MaximizeParsimony(
      ds,
      concavity = conc,
      maxReplicates = 50L,
      targetHits = 10L,
      maxSeconds = as.double(TIMEOUT),
      strategy = cond$strategy,
      verbosity = 0L,
      nThreads = 1L
    ),
    error = function(e) {
      warning("Search error: ", conditionMessage(e))
      NULL
    }
  )
  wall_s <- as.double((proc.time() - t0)[3])

  if (is.null(result)) {
    cat("ERROR\n")
    timings <- setNames(rep(NA_real_, length(TIMING_NAMES)), TIMING_NAMES)
    row <- data.frame(
      dataset = cond$dataset,
      n_taxa = ds_info$ntax,
      n_chars = ds_info$nchar,
      pct_inapp = ds_info$pct_inapp,
      scoring = cond$scoring,
      weighting = cond$weighting,
      concavity = conc,
      strategy = cond$strategy,
      seed = cond$seed,
      score = NA_real_,
      n_trees = NA_integer_,
      wall_s = wall_s,
      replicates = NA_integer_,
      hits = NA_integer_,
      stringsAsFactors = FALSE
    )
    for (tn in TIMING_NAMES) row[[tn]] <- NA_real_
    row$rep_scores <- NA_character_
    rows[[i]] <- row
    next
  }

  # Extract timing
  timings_vec <- attr(result, "timings")
  rep_scores <- attr(result, "replicate_scores")
  score <- attr(result, "score")
  reps <- attr(result, "replicates")
  hits <- attr(result, "hits_to_best")

  cat(sprintf("score=%.2f reps=%d time=%.1fs\n", score, reps, wall_s))

  row <- data.frame(
    dataset = cond$dataset,
    n_taxa = ds_info$ntax,
    n_chars = ds_info$nchar,
    pct_inapp = ds_info$pct_inapp,
    scoring = cond$scoring,
    weighting = cond$weighting,
    concavity = conc,
    strategy = cond$strategy,
    seed = cond$seed,
    score = score,
    n_trees = length(result),
    wall_s = round(wall_s, 3),
    replicates = reps,
    hits = hits,
    stringsAsFactors = FALSE
  )

  # Add per-phase timings
  for (tn in TIMING_NAMES) {
    row[[tn]] <- if (tn %in% names(timings_vec)) timings_vec[[tn]] else 0
  }

  # Serialize replicate scores as comma-separated string
  row$rep_scores <- paste(round(rep_scores, 4), collapse = ",")

  rows[[i]] <- row
}

results <- do.call(rbind, rows)

# ---- Save raw results ----
ts <- format(Sys.time(), "%Y%m%d_%H%M")
outfile <- file.path(outdir, sprintf("t290b_phase_%s.csv", ts))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nRaw results saved: %s (%d rows)\n\n", outfile, nrow(results)))

# ---- Phase summary ----
cat("=== Phase Time Distribution (% of wall time, by scoring x weighting) ===\n\n")

phase_cols <- TIMING_NAMES

for (sc in SCORING) {
  for (wt in names(WEIGHTINGS)) {
    sub <- results[results$scoring == sc & results$weighting == wt, ]
    if (nrow(sub) == 0) next

    total_ms <- rowSums(sub[, phase_cols], na.rm = TRUE)
    pct <- colMeans(
      sweep(sub[, phase_cols], 1, total_ms, "/") * 100,
      na.rm = TRUE
    )
    pct <- pct[pct > 0.5]  # skip negligible phases

    cat(sprintf("--- %s x %s (n=%d runs, mean wall=%.1fs) ---\n",
                sc, wt, nrow(sub), mean(sub$wall_s)))
    for (p in names(sort(pct, decreasing = TRUE))) {
      cat(sprintf("  %-20s %5.1f%%\n", p, pct[p]))
    }
    cat("\n")
  }
}

# ---- Brazeau vs Fitch cost ratio per phase ----
cat("=== Brazeau / Fitch Cost Ratio by Phase (EW, matched datasets) ===\n\n")

for (strat in STRATEGIES) {
  cat(sprintf("Strategy: %s\n", strat))
  ew_braz <- results[results$scoring == "brazeau" & results$weighting == "EW" &
                       results$strategy == strat, ]
  ew_fitch <- results[results$scoring == "fitch" & results$weighting == "EW" &
                        results$strategy == strat, ]

  # Aggregate by dataset (mean across seeds)
  braz_agg <- aggregate(. ~ dataset, data = ew_braz[, c("dataset", phase_cols)],
                        FUN = mean)
  fitch_agg <- aggregate(. ~ dataset, data = ew_fitch[, c("dataset", phase_cols)],
                          FUN = mean)

  merged <- merge(braz_agg, fitch_agg, by = "dataset", suffixes = c(".braz", ".fitch"))
  for (p in phase_cols) {
    braz_total <- sum(merged[[paste0(p, ".braz")]], na.rm = TRUE)
    fitch_total <- sum(merged[[paste0(p, ".fitch")]], na.rm = TRUE)
    if (fitch_total > 100) {
      cat(sprintf("  %-20s Braz: %7.0f ms  Fitch: %7.0f ms  ratio: %.2fx\n",
                  p, braz_total, fitch_total, braz_total / fitch_total))
    }
  }
  cat("\n")
}

cat(sprintf("\n=== Completed: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
