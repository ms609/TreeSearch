# T-249: Rerun TNT comparison on current cpp-search HEAD
#
# Uses the agent-F build (TreeSearch.F) to benchmark the current engine
# against TNT on the round 2 dataset set.
#
# Run from the TS-TNT-bench directory:
#   Rscript dev/benchmarks/run_t249_comparison.R

# Load the agent build
agent_lib <- normalizePath("../.builds/TreeSearch-F", winslash = "/")
.libPaths(c(agent_lib, .libPaths()))
library(TreeSearch.F)

# Namespace shim: make TreeSearch::: resolve to TreeSearch.F
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))

library(TreeTools)

# Source benchmark helpers, skipping library() and nested source() calls
source_filtered <- function(path) {
  lines <- readLines(path)
  for (i in seq_along(lines)) {
    if (grepl("library(TreeSearch)", lines[i], fixed = TRUE) ||
        grepl("library(TreeTools)", lines[i], fixed = TRUE) ||
        grepl('source("dev/benchmarks/bench_datasets.R")', lines[i], fixed = TRUE)) {
      lines[i] <- paste0("# SKIP: ", lines[i])
    }
  }
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp))
  writeLines(lines, tmp)
  source(tmp, local = FALSE)
}
source_filtered("dev/benchmarks/bench_datasets.R")
source_filtered("dev/benchmarks/bench_tnt_compare.R")

# --- Round 2 datasets (same as original round2_hard.csv) ---
ROUND2_NAMES <- c(
  "Conrad2008",     # 64 tips
  "Liljeblad2008",  # 68 tips
  "Aguado2009",     # 43 tips
  "Capa2011",       # 61 tips
  "Eklund2004",     # 54 tips
  "Geisler2001",    # 73 tips  (gap: +5)
  "Giles2015",      # 78 tips
  "OMeara2014",     # 49 tips
  "Rougier2012",    # 75 tips
  "Shultz2007",     # 52 tips
  "Wetterer2000",   # 42 tips
  "Wilson2003",     # 53 tips
  "Wortley2006",    # 37 tips  (gap: +3)
  "Zanol2014",      # 74 tips  (gap: +2)
  "Zhu2013",        # 75 tips  (gap: +3)
  "Dikow2009"       # 88 tips
)

cat("=== T-249: TNT vs TreeSearch comparison (current cpp-search HEAD) ===\n")
cat("Engine: TreeSearch.F from agent build\n")
cat("Datasets:", length(ROUND2_NAMES), "\n")
cat("Timeouts: 120s, Seeds: 1:3, EW only\n")
cat("Start:", format(Sys.time()), "\n\n")

results <- run_comparison(
  dataset_names = ROUND2_NAMES,
  weightings = "EW",
  timeout_s = 120,
  seeds = 1:3,
  hits = 5L,
  reps = 50L
)

cat("\n=== Done:", format(Sys.time()), "===\n\n")

# Save results
outfile <- paste0("dev/benchmarks/round3_t249_",
                  format(Sys.time(), "%Y%m%d_%H%M"), ".csv")
write.csv(results, outfile, row.names = FALSE)
cat("Results saved to:", outfile, "\n")

# Quick summary: score gap per dataset
cat("\n--- Score gap summary (TS - TNT, best of 3 seeds) ---\n")
agg <- aggregate(cbind(tnt_score, ts_score) ~ dataset, data = results, FUN = min)
agg$gap <- agg$ts_score - agg$tnt_score
agg$pct <- round(100 * agg$gap / agg$tnt_score, 2)
agg <- agg[order(agg$pct, decreasing = TRUE), ]
print(agg, row.names = FALSE)

cat("\nDatasets with gap > 0:\n")
print(agg[agg$gap > 0, ], row.names = FALSE)
