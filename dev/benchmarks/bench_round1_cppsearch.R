#!/usr/bin/env Rscript
# Round 1: 14-dataset benchmark — TreeSearch (cpp-search) vs TNT 1.6
# Uses tuned "default" strategy parameters from MaximizeParsimony.
# 10s timeout, 3 seeds, hits=5, reps=20, single-threaded.
#
# Usage: Rscript bench_round1_cppsearch.R
# Requires: TreeSearch.F built in .builds/TreeSearch-F/

cat("=== Round 1: 14-dataset benchmark (cpp-search) ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Setup ---
AGENT_LIB <- normalizePath("../../../.builds/TreeSearch-F", mustWork = TRUE)
.libPaths(c(AGENT_LIB, .libPaths()))

library(TreeSearch.F)
# Namespace shim so TreeSearch::: calls resolve
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))
library(TreeTools)

TNT_EXE <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
stopifnot("TNT not found" = file.exists(TNT_EXE))

STAGING_DIR <- normalizePath("../../.tnt-bench", mustWork = FALSE)
dir.create(STAGING_DIR, showWarnings = FALSE)

OUTFILE <- sprintf("round1_cppsearch_%s.csv", format(Sys.time(), "%Y%m%d_%H%M"))

# --- Helpers ---
clean_inapplicable <- function(phyDat_obj) {
  mat <- PhyDatToMatrix(phyDat_obj)
  mat[mat == "-"] <- "?"
  MatrixToPhyDat(mat)
}

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

run_tnt <- function(data_file, timeout_s, seed, hits, reps) {
  commands <- c(
    "mxram 1024;",
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed),
    sprintf("timeout %d:%02d:%02d;",
            timeout_s %/% 3600, (timeout_s %% 3600) %/% 60, timeout_s %% 60),
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
  m <- regmatches(out_text, regexpr("Best score:\\s+([0-9.]+)", out_text))
  if (length(m) == 1) score <- as.numeric(sub("Best score:\\s+", "", m))

  n_trees <- NA_integer_
  m <- regmatches(out_text, regexpr("([0-9]+) trees? retained", out_text))
  if (length(m) == 1) n_trees <- as.integer(sub(" trees? retained", "", m))

  list(score = score, n_trees = n_trees, wall_s = wall_s)
}

run_treesearch <- function(ds, timeout_s, seed, hits, reps) {
  set.seed(seed)
  t0 <- proc.time()
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = as.integer(reps),
    targetHits = as.integer(hits),
    maxSeconds = as.double(timeout_s),
    verbosity = 0L, nThreads = 1L,
    # Tuned "default" strategy from MaximizeParsimony
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
BENCHMARK_NAMES <- c(
  "Longrich2010", "Vinther2008", "Sansom2010", "DeAssis2011",
  "Aria2015", "Wortley2006", "Griswold1999", "Schulze2007",
  "Eklund2004", "Agnarsson2004", "Zanol2014", "Zhu2013",
  "Giles2015", "Dikow2009"
)
TIMEOUT <- 10
SEEDS <- 1:3
HITS <- 5L
REPS <- 20L

# --- Prepare datasets ---
cat("Preparing datasets...\n")
datasets <- list()
for (nm in BENCHMARK_NAMES) {
  raw <- inapplicable.phyData[[nm]]
  clean <- clean_inapplicable(raw)
  datasets[[nm]] <- prepare_ts_data(clean)
  WriteTntCharacters(clean,
                     filepath = file.path(STAGING_DIR, paste0(nm, ".tnt")))
}

# --- Run benchmark ---
rows <- list()
idx <- 0L
n_total <- length(BENCHMARK_NAMES) * length(SEEDS)

for (nm in BENCHMARK_NAMES) {
  ds <- datasets[[nm]]
  data_file <- paste0(nm, ".tnt")

  for (seed in SEEDS) {
    idx <- idx + 1L
    cat(sprintf("[%d/%d] %s seed=%d ... ", idx, n_total, nm, seed))

    tnt <- run_tnt(data_file, timeout_s = TIMEOUT, seed = seed,
                   hits = HITS, reps = REPS)
    ts <- run_treesearch(ds, timeout_s = TIMEOUT, seed = seed,
                         hits = HITS, reps = REPS)

    rows[[idx]] <- data.frame(
      dataset = nm, n_taxa = ds$n_taxa, n_chars = ds$n_chars,
      seed = seed, timeout_s = TIMEOUT,
      tnt_score = tnt$score, tnt_trees = tnt$n_trees,
      tnt_wall_s = round(tnt$wall_s, 3),
      ts_score = ts$score, ts_trees = ts$n_trees,
      ts_wall_s = round(ts$wall_s, 3),
      ts_reps = ts$replicates, ts_hits = ts$hits,
      stringsAsFactors = FALSE
    )

    cat(sprintf("TNT=%.0f (%.1fs)  TS=%.0f (%.1fs)\n",
                tnt$score, tnt$wall_s, ts$score, ts$wall_s))

    # Incremental save
    results <- do.call(rbind, rows)
    write.csv(results, OUTFILE, row.names = FALSE)
  }
}

cat("\n=== Round 1 complete ===\n")
cat("Results saved to:", OUTFILE, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- Quick summary ---
results <- do.call(rbind, rows)
cat("\nScore comparison (best across seeds):\n")
agg <- aggregate(cbind(tnt_score, ts_score) ~ dataset + n_taxa, data = results, FUN = min)
agg$delta <- agg$ts_score - agg$tnt_score
agg <- agg[order(agg$n_taxa), ]
print(agg, row.names = FALSE)
cat("\nTS wins:", sum(agg$delta < 0),
    " Ties:", sum(agg$delta == 0),
    " TNT wins:", sum(agg$delta > 0), "\n")
