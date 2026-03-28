#!/usr/bin/env Rscript
# T-290c (quick): wagnerStarts = 1 vs 3 under Brazeau — two-regime validation
#
# Tests the analytical prediction from T-290 that:
#  - Regime A (0 reps in 30s, e.g. project2084_(1) 86t/3660c):
#    wagnerStarts=3 provides better starting topology → better score
#  - Regime B (multiple reps, e.g. project2086 91t/453c):
#    wagnerStarts=1 and 3 produce equivalent scores
#
# EW only, 30s and 60s, 5 seeds.
#
# Usage:
#   Rscript bench_t290c_quick.R <output_dir>

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args) >= 1) args[1] else "dev/benchmarks/t290_results"

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
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("=== T-290c (quick): wagnerStarts validation ===\n")
cat("Repo root:", repo_root, "\n")
cat("Output dir:", outdir, "\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

library(TreeSearch)
library(TreeTools)
source("dev/benchmarks/bench_datasets.R")

# ---- Two representative datasets ----
# Regime A: won't complete a rep in 30s → starting quality dominates
# Regime B: completes ~7 reps in 30s → ratchet dominates
DATASETS <- c("project2084_(1)",  # 86t, 3660c, Regime A
              "project2086")       # 91t,  453c, Regime B

TIMEOUTS <- c(30, 60)
N_SEEDS  <- 5L

# Thorough preset (current), varying only wagnerStarts
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
for (k in DATASETS) {
  row <- cat_df[cat_df$key == k, , drop = FALSE]
  nex_path <- file.path(NEOTRANS_MATRICES_DIR, row$filename[1])
  datasets[[k]] <- suppressWarnings(TreeTools::ReadAsPhyDat(nex_path))
  cat(sprintf("Loaded %s: %dt %dc %.0f%%inapp\n",
              k, row$ntax, row$nchar, row$pct_inapp))
}
cat("\n")

# ---- Run ----
grid <- expand.grid(
  dataset   = DATASETS,
  condition = names(CONDITIONS),
  timeout   = TIMEOUTS,
  seed      = seq_len(N_SEEDS),
  stringsAsFactors = FALSE
)
n_total <- nrow(grid)
cat(sprintf("Grid: %d runs (2 datasets × 2 conditions × 2 timeouts × %d seeds)\n\n",
            n_total, N_SEEDS))

rows <- vector("list", n_total)

for (i in seq_len(n_total)) {
  g    <- grid[i, ]
  ds   <- datasets[[g$dataset]]
  ctrl <- CONDITIONS[[g$condition]]
  row_info <- cat_df[cat_df$key == g$dataset, , drop = FALSE]

  cat(sprintf("[%d/%d] %s wStarts=%d %ds seed=%d ... ",
              i, n_total, g$dataset, ctrl$wagnerStarts, g$timeout, g$seed))

  set.seed(g$seed)
  t0 <- proc.time()
  result <- tryCatch(
    do.call(MaximizeParsimony, c(
      list(dataset       = ds,
           concavity     = Inf,   # EW
           maxReplicates = 200L,
           targetHits    = 10L,
           maxSeconds    = as.double(g$timeout),
           verbosity     = 0L,
           nThreads      = 1L),
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
      dataset = g$dataset, n_taxa = row_info$ntax,
      n_chars = row_info$nchar,
      condition = g$condition, timeout = g$timeout, seed = g$seed,
      score = NA_real_, replicates = NA_integer_, hits = NA_integer_,
      wall_s = wall_s, stringsAsFactors = FALSE
    )
    next
  }

  score <- attr(result, "score")
  reps  <- attr(result, "replicates")
  hits  <- attr(result, "hits_to_best")
  cat(sprintf("score=%.0f reps=%d time=%.1fs\n", score, reps, wall_s))

  rows[[i]] <- data.frame(
    dataset = g$dataset, n_taxa = row_info$ntax,
    n_chars = row_info$nchar,
    condition = g$condition, timeout = g$timeout, seed = g$seed,
    score = score, replicates = reps, hits = hits,
    wall_s = round(wall_s, 3), stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)

# ---- Save ----
ts <- format(Sys.time(), "%Y%m%d_%H%M")
outfile <- file.path(outdir, sprintf("t290c_quick_%s.csv", ts))
write.csv(results, outfile, row.names = FALSE)
cat(sprintf("\nSaved: %s\n\n", outfile))

# ---- Summary ----
cat("=== wagnerStarts 1 vs 3: Two-regime comparison (EW) ===\n\n")

for (to in TIMEOUTS) {
  cat(sprintf("--- %ds timeout ---\n", to))
  sub <- results[results$timeout == to & !is.na(results$score), ]

  by_cond <- aggregate(
    cbind(score, replicates) ~ dataset + n_taxa + n_chars + condition,
    data = sub, FUN = function(x) c(best = min(x), med = median(x))
  )

  # Print per-dataset comparison
  for (ds in DATASETS) {
    d <- sub[sub$dataset == ds, ]
    w1 <- d[d$condition == "w1", ]
    w3 <- d[d$condition == "w3", ]
    if (nrow(w1) == 0 || nrow(w3) == 0) next
    cat(sprintf("  %s (%dt %dc):\n", ds,
                row_info$ntax[1],
                results[results$dataset == ds, "n_chars"][1]))
    cat(sprintf("    wagnerStarts=1: best=%.0f  median=%.0f  reps=%.1f\n",
                min(w1$score, na.rm = TRUE),
                median(w1$score, na.rm = TRUE),
                median(w1$replicates, na.rm = TRUE)))
    cat(sprintf("    wagnerStarts=3: best=%.0f  median=%.0f  reps=%.1f\n",
                min(w3$score, na.rm = TRUE),
                median(w3$score, na.rm = TRUE),
                median(w3$replicates, na.rm = TRUE)))
    cat(sprintf("    Delta best (w1-w3): %.0f  [positive = w3 better]\n",
                min(w1$score, na.rm = TRUE) - min(w3$score, na.rm = TRUE)))
  }
  cat("\n")
}

cat(sprintf("=== Completed: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
