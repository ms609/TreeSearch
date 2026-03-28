#!/usr/bin/env Rscript
# T-290c: wagnerStarts = 1 vs 3 under Brazeau scoring
#
# Tests whether the thorough preset's wagnerStarts=3 (each start 3.6–5.2×
# more expensive than Fitch) is worthwhile under Brazeau scoring compared
# to wagnerStarts=1 with the same time budget.
#
# Holds all thorough preset params fixed except wagnerStarts.
# Runs on 6 representative Brazeau-sample datasets spanning 23–173 tips.
# EW + IW(k=10), 30s + 60s, 5 seeds each.
#
# Usage:
#   Rscript bench_wagner_starts_brazeau.R <output_dir>

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

cat("=== T-290c: wagnerStarts 1 vs 3 (Brazeau) ===\n")
cat("Repo root:", repo_root, "\n")
cat("Output dir:", outdir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")

# ---- Configuration ----
# Same 6 datasets as T-290b phase profiling
SUBSET <- c("project2346", "project2359", "project4146_(3)",
            "project2086", "project2771", "project804")

TIMEOUTS <- c(30, 60)
N_SEEDS  <- 5L
WEIGHTINGS <- list(EW = Inf, IW10 = 10)

# Thorough preset with wagnerStarts as the only variable
THOROUGH_BASE <- list(
  tbrMaxHits             = 3L,
  tabuSize               = 200L,
  ratchetCycles          = 20L,
  ratchetPerturbProb     = 0.25,
  ratchetPerturbMode     = 2L,
  ratchetPerturbMaxMoves = 5L,
  ratchetAdaptive        = FALSE,
  driftCycles            = 0L,
  xssRounds              = 3L,
  xssPartitions          = 4L,
  rssRounds              = 2L,
  cssRounds              = 1L,
  sectorMinSize          = 6L,
  sectorMaxSize          = 80L,
  fuseInterval           = 2L,
  fuseAcceptEqual        = TRUE,
  outerCycles            = 2L,
  nniFirst               = TRUE
)

CONDITIONS <- list(
  w1 = c(THOROUGH_BASE, list(wagnerStarts = 1L)),
  w3 = c(THOROUGH_BASE, list(wagnerStarts = 3L))
)

# ---- Load datasets ----
cat_df <- load_mbank_catalogue()
datasets <- list()
for (k in SUBSET) {
  row <- cat_df[cat_df$key == k, ]
  nex_path <- file.path(NEOTRANS_MATRICES_DIR, row$filename)
  datasets[[k]] <- suppressWarnings(TreeTools::ReadAsPhyDat(nex_path))
}
cat(sprintf("Loaded %d datasets\n\n", length(datasets)))

# ---- Run grid ----
grid <- expand.grid(
  dataset    = SUBSET,
  condition  = names(CONDITIONS),
  weighting  = names(WEIGHTINGS),
  timeout    = TIMEOUTS,
  seed       = seq_len(N_SEEDS),
  stringsAsFactors = FALSE
)
n_total <- nrow(grid)
cat(sprintf("Grid: %d runs (%d datasets × %d conditions × %d weights × %d timeouts × %d seeds)\n\n",
            n_total, length(SUBSET), length(CONDITIONS),
            length(WEIGHTINGS), length(TIMEOUTS), N_SEEDS))

rows <- vector("list", n_total)

for (i in seq_len(n_total)) {
  g     <- grid[i, ]
  ds    <- datasets[[g$dataset]]
  ctrl  <- CONDITIONS[[g$condition]]
  conc  <- WEIGHTINGS[[g$weighting]]
  ds_info <- cat_df[cat_df$key == g$dataset, ]

  cat(sprintf("[%d/%d] %s wStarts=%s %s %ds seed=%d ... ",
              i, n_total, g$dataset,
              ctrl$wagnerStarts, g$weighting, g$timeout, g$seed))

  set.seed(g$seed)
  t0 <- proc.time()
  result <- tryCatch(
    do.call(MaximizeParsimony, c(
      list(dataset     = ds,
           concavity   = conc,
           maxReplicates = 200L,
           targetHits  = 10L,
           maxSeconds  = as.double(g$timeout),
           verbosity   = 0L,
           nThreads    = 1L),
      ctrl
    )),
    error = function(e) {
      warning("Search error: ", conditionMessage(e))
      NULL
    }
  )
  wall_s <- as.double((proc.time() - t0)[3])

  if (is.null(result)) {
    cat("ERROR\n")
    rows[[i]] <- data.frame(
      dataset    = g$dataset,
      n_taxa     = ds_info$ntax,
      condition  = g$condition,
      weighting  = g$weighting,
      timeout    = g$timeout,
      seed       = g$seed,
      score      = NA_real_,
      replicates = NA_integer_,
      hits       = NA_integer_,
      wall_s     = wall_s,
      stringsAsFactors = FALSE
    )
    next
  }

  score <- attr(result, "score")
  reps  <- attr(result, "replicates")
  hits  <- attr(result, "hits_to_best")
  cat(sprintf("score=%.2f reps=%d time=%.1fs\n", score, reps, wall_s))

  rows[[i]] <- data.frame(
    dataset    = g$dataset,
    n_taxa     = ds_info$ntax,
    condition  = g$condition,
    weighting  = g$weighting,
    timeout    = g$timeout,
    seed       = g$seed,
    score      = score,
    replicates = reps,
    hits       = hits,
    wall_s     = round(wall_s, 3),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)

# ---- Save ----
ts <- format(Sys.time(), "%Y%m%d_%H%M")
outfile <- file.path(outdir, sprintf("t290c_wagner_starts_%s.csv", ts))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nSaved: %s (%d rows)\n\n", outfile, nrow(results)))

# ---- Summary ----
cat("=== Score Comparison: wagnerStarts = 1 vs 3 (Brazeau) ===\n\n")

library(dplyr)

for (wt in names(WEIGHTINGS)) {
  for (to in TIMEOUTS) {
    sub <- results |>
      filter(weighting == wt, timeout == to, !is.na(score))

    summ <- sub |>
      group_by(dataset, n_taxa, condition) |>
      summarise(
        best   = min(score),
        median = median(score),
        reps   = median(replicates),
        .groups = "drop"
      ) |>
      tidyr::pivot_wider(names_from = condition,
                         values_from = c(best, median, reps)) |>
      mutate(
        best_delta   = best_w1 - best_w3,    # positive = w3 better
        median_delta = median_w1 - median_w3
      ) |>
      arrange(n_taxa)

    cat(sprintf("--- %s, %ds ---\n", wt, to))
    print(summ, n = 20)
    cat("\n")
    cat(sprintf("  w3 better (best): %d/%d  mean delta: %.3f\n",
                sum(summ$best_delta > 0, na.rm = TRUE), nrow(summ),
                mean(summ$best_delta, na.rm = TRUE)))
    cat(sprintf("  Replicate ratio (w1/w3): %.2f\n\n",
                mean(summ$reps_w1 / summ$reps_w3, na.rm = TRUE)))
  }
}

cat(sprintf("\n=== Completed: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
