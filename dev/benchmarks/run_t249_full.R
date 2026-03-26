# T-249: Full TNT vs TreeSearch comparison on current cpp-search HEAD
# 16 datasets from round 2, EW, 120s timeout, 3 seeds

agent_lib <- normalizePath("../.builds/TreeSearch-F", winslash = "/")
.libPaths(c(agent_lib, .libPaths()))
library(TreeSearch.F)
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))
library(TreeTools)

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

# Full round 2 dataset list (same order as original)
ROUND2_NAMES <- c(
  "Conrad2008",     # 64 taxa
  "Liljeblad2008",  # 68 taxa
  "OMeara2014",     # 63 taxa
  "Rougier2012",    # 58 taxa
  "Geisler2001",    # 68 taxa
  "Shultz2007",     # 59 taxa
  "Aguado2009",     # 76 taxa
  "Capa2011",       # 67 taxa
  "Wilson2003",     # 61 taxa
  "Wetterer2000",   # 63 taxa
  "Wortley2006",    # 37 taxa
  "Eklund2004",     # 54 taxa
  "Zanol2014",      # 74 taxa
  "Zhu2013",        # 75 taxa
  "Giles2015",      # 78 taxa
  "Dikow2009"       # 88 taxa
)

cat("=== T-249: Full round 3 comparison (current cpp-search HEAD) ===\n")
cat("Engine: TreeSearch.F via MaximizeParsimony (auto strategy)\n")
cat("Datasets:", length(ROUND2_NAMES), "\n")
cat("Timeouts: 120s, Seeds: 1:3, EW only, Fitch mode (fair TNT comparison)\n")
cat("Start:", format(Sys.time()), "\n\n")

results <- run_comparison(
  dataset_names = ROUND2_NAMES,
  weightings = "EW",
  timeout_s = 120,
  seeds = 1:3,
  hits = 10L,
  reps = 50L,
  use_fitch = TRUE
)

cat("\n=== Done:", format(Sys.time()), "===\n\n")

outfile <- "dev/benchmarks/round3_t249.csv"
write.csv(results, outfile, row.names = FALSE)
cat("Results saved to:", outfile, "\n")

# --- Analysis ---
r2 <- read.csv("dev/benchmarks/round2_hard.csv")

# Best score per dataset (min across seeds)
agg_new <- aggregate(cbind(tnt_score, ts_score) ~ dataset, data = results, FUN = min)
agg_new$gap_r3 <- agg_new$ts_score - agg_new$tnt_score
agg_new$pct_r3 <- round(100 * agg_new$gap_r3 / agg_new$tnt_score, 2)

agg_old <- aggregate(cbind(tnt_score, ts_score) ~ dataset, data = r2, FUN = min)
agg_old$gap_r2 <- agg_old$ts_score - agg_old$tnt_score

merged <- merge(agg_new[, c("dataset", "gap_r3", "pct_r3", "ts_score")],
                agg_old[, c("dataset", "gap_r2")],
                by = "dataset", all.x = TRUE)
merged$change <- merged$gap_r3 - merged$gap_r2
merged <- merged[order(merged$pct_r3, decreasing = TRUE), ]

cat("\n--- Round 2 vs Round 3 comparison ---\n")
cat("gap = TS best score - TNT best score (lower is better)\n\n")
print(merged, row.names = FALSE)

cat("\n--- Summary ---\n")
cat("Datasets where gap improved:", sum(merged$change < 0, na.rm = TRUE), "\n")
cat("Datasets where gap unchanged:", sum(merged$change == 0, na.rm = TRUE), "\n")
cat("Datasets where gap worsened:", sum(merged$change > 0, na.rm = TRUE), "\n")
cat("Mean gap (R3):", round(mean(merged$gap_r3), 1), "steps\n")
cat("Mean gap (R2):", round(mean(merged$gap_r2, na.rm = TRUE), 1), "steps\n")
cat("Mean pct gap (R3):", round(mean(merged$pct_r3), 2), "%\n")

# Wall-clock comparison
time_new <- aggregate(cbind(tnt_wall_s, ts_wall_s) ~ dataset, data = results,
                      FUN = median)
cat("\n--- Median wall-clock times ---\n")
print(time_new, row.names = FALSE)
cat("\nMedian TNT time:", round(median(time_new$tnt_wall_s), 1), "s\n")
cat("Median TS time:", round(median(time_new$ts_wall_s), 1), "s\n")
cat("Speed ratio (TS/TNT):", round(median(time_new$ts_wall_s / time_new$tnt_wall_s), 1), "x\n")
