# bench_sector_scaling_ab.R
#
# Matched-wall A/B for the TNT-style sector-size SCALING knobs added to
# rss_search() (branch cpp-search; env-gated, default-OFF):
#   (a) nTip-scaled cap  -- TS_SECT_MAXSIZE / TS_SECT_MAXFRAC / TS_SECT_THRESHOLD
#   (b) adaptive growth  -- TS_SECT_GROW(+_INC/_SELFACT/_MOVEON/_START)
# See dev/plans/2026-07-08-sector-size-scaling.md for the design + env reference.
#
# METHOD (per the mission's evidence-first rule):
#   * TRAINING split only.  The validation / project175 split is SEQUESTERED --
#     never load split = "validation" here, and never tune knob values on it.
#   * Matched wall.  Every arm gets the same maxSeconds budget; the reported
#     best_score at that budget IS the anytime metric.  Sweep several budgets
#     to trace the anytime curve.
#   * rasStarts >= 2, so the RAS-retention sector channel is not inert
#     (thorough/large default rasStarts = 1 makes it inert).
#   * Large matrices (nTip >= 120: tiers large + xlarge), the `large` preset's
#     home and where feature (a) is designed to bind.
#
# Run one cell (Hamilton array element) or the whole grid from package root:
#   Rscript dev/benchmarks/bench_sector_scaling_ab.R [budget_s] [n_matrices] [lib]
# Env for a clean measurement: TERM=dumb OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
#
# Results appended to dev/benchmarks/results_sector_scaling_ab.csv

args     <- commandArgs(trailingOnly = TRUE)
BUDGET_S <- if (length(args) >= 1L) as.double(args[[1]]) else 120
N_MATRIX <- if (length(args) >= 2L) as.integer(args[[2]]) else 6L
LIB      <- if (length(args) >= 3L) args[[3]] else ".agent-X"
SEEDS    <- c(1031L, 2847L, 7193L, 4561L, 8822L)

.libPaths(c(LIB, .libPaths()))
library(TreeSearch)
library(TreeTools)

SRC <- getwd()
source(file.path(SRC, "dev/benchmarks/bench_datasets.R"))
# Pull presets from source (pure-R; no rebuild needed to adjust the recipe).
source(file.path(SRC, "R/SearchControl.R"))
source(file.path(SRC, "R/MaximizeParsimony.R"))

OUT_FILE <- file.path(SRC, "dev/benchmarks/results_sector_scaling_ab.csv")

# ---- Env-arm matrix.  These knob VALUES are the tuning surface; edit them and
# re-run on the TRAINING split.  An arm is a named character vector of env vars
# (empty = baseline / feature off).  min() means the preset cap still bounds the
# fraction, so `capraise` must accompany `frac` for it to bind on large trees.
ARMS <- list(
  baseline = character(0),
  # (a) raised fixed ceiling alone -- isolates "bigger sectors" from the scaling
  capraise = c(TS_SECT_MAXSIZE = "200"),
  # (a) nTip-scaled cap: ~nTip/2 once nTip >= 88, bounded by the raised ceiling
  frac     = c(TS_SECT_MAXSIZE = "200", TS_SECT_MAXFRAC = "0.5",
               TS_SECT_THRESHOLD = "88"),
  # (b) adaptive growth from small sectors up to the (preset) cap
  grow     = c(TS_SECT_GROW = "1", TS_SECT_GROW_INC = "50",
               TS_SECT_GROW_SELFACT = "40", TS_SECT_GROW_START = "6",
               TS_SECT_GROW_MOVEON = "0"),
  # (a)+(b) combined: grow up to the nTip-scaled ceiling
  both     = c(TS_SECT_MAXSIZE = "200", TS_SECT_MAXFRAC = "0.5",
               TS_SECT_THRESHOLD = "88", TS_SECT_GROW = "1",
               TS_SECT_GROW_INC = "50", TS_SECT_GROW_SELFACT = "40",
               TS_SECT_GROW_START = "6", TS_SECT_GROW_MOVEON = "0")
)

set_env <- function(vars) {
  clear_scaling_env()
  if (length(vars)) do.call(Sys.setenv, as.list(vars))
}
clear_scaling_env <- function() {
  Sys.unsetenv(c("TS_SECT_MAXSIZE", "TS_SECT_MAXFRAC", "TS_SECT_THRESHOLD",
                 "TS_SECT_GROW", "TS_SECT_GROW_INC", "TS_SECT_GROW_SELFACT",
                 "TS_SECT_GROW_MOVEON", "TS_SECT_GROW_START"))
}
on.exit(clear_scaling_env(), add = TRUE)

# ---- Datasets: TRAINING split, large + xlarge tiers (nTip >= 120). ----
catalogue <- load_mbank_catalogue()
pool <- catalogue[catalogue$split == "training" &
                    catalogue$tier %in% c("large", "xlarge"), ]
if (nrow(pool) == 0) stop("No large/xlarge training matrices in catalogue")
pool <- pool[order(-pool$ntax), ]
keys <- head(pool$key, N_MATRIX)
datasets <- load_mbank_datasets(catalogue, keys, verbose = TRUE)
cat(sprintf("\nTreeSearch %s | budget %ds | %d matrices x %d arms x %d seeds\n\n",
            as.character(packageVersion("TreeSearch")), BUDGET_S,
            length(datasets), length(ARMS), length(SEEDS)))

# `large` preset with rasStarts forced >= 2 (see header).
strat <- unclass(.StrategyPresets()[["large"]])
attr(strat, "class") <- NULL
strat$rasStarts <- 2L

rows <- list(); idx <- 0L
for (dkey in names(datasets)) {
  ds <- datasets[[dkey]]
  for (arm in names(ARMS)) {
    for (seed in SEEDS) {
      idx <- idx + 1L
      set_env(ARMS[[arm]])
      t0 <- proc.time(); set.seed(seed)
      res <- tryCatch(
        do.call(TreeSearch:::ts_driven_search,
                c(list(contrast = ds$contrast, tip_data = ds$tip_data,
                       weight = ds$weight, levels = ds$levels,
                       maxReplicates = 1000L,
                       targetHits = max(10L, ds$n_taxa %/% 5L),
                       maxSeconds = as.double(BUDGET_S), verbosity = 0L),
                  strat)),
        error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
      wall <- as.double((proc.time() - t0)[3])
      clear_scaling_env()
      if (is.null(res)) next
      cat(sprintf("[%d] %-22s %-8s seed %d: score=%.0f reps=%d wall=%.1fs\n",
                  idx, dkey, arm, seed, res$best_score, res$replicates, wall))
      rows[[idx]] <- data.frame(
        dataset = dkey, n_taxa = ds$n_taxa, arm = arm, seed = seed,
        budget_s = BUDGET_S, best_score = res$best_score,
        replicates = res$replicates, hits_to_best = res$hits_to_best,
        wall_s = wall, stringsAsFactors = FALSE)
    }
  }
}

df <- do.call(rbind, rows)
write.table(df, OUT_FILE, sep = ",", row.names = FALSE,
            col.names = !file.exists(OUT_FILE), append = file.exists(OUT_FILE))
cat("\nAppended", nrow(df), "rows to", OUT_FILE, "\n")

# ---- Summary: median best_score per arm (lower = better), delta vs baseline.
cat(sprintf("\n===== median best_score by arm (budget %ds) =====\n", BUDGET_S))
base_med <- median(df$best_score[df$arm == "baseline"], na.rm = TRUE)
for (a in names(ARMS)) {
  r <- df$best_score[df$arm == a]
  if (!length(r)) next
  cat(sprintf("  %-8s  median=%8.1f  delta_vs_base=%+7.1f  (n=%d)\n",
              a, median(r, na.rm = TRUE), median(r, na.rm = TRUE) - base_med,
              sum(!is.na(r))))
}
cat("\nNOTE: adopt an arm only on a matched-wall WIN with no regression, and\n",
    "only after confirming these knob values were tuned on the TRAINING split.\n")
