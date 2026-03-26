# T-249: Focus on datasets with known score gaps from round 2
# Geisler2001 (+5), Wortley2006 (+3), Zhu2013 (+3), Zanol2014 (+2)
# Plus Conrad2008 (+1-2) as a borderline case

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
cat("bench_datasets.R loaded\n")
source_filtered("dev/benchmarks/bench_tnt_compare.R")
cat("bench_tnt_compare.R loaded\n")

GAP_DATASETS <- c(
  "Geisler2001",   # 73 tips, gap +5
  "Wortley2006",   # 37 tips, gap +3
  "Zhu2013",       # 75 tips, gap +3
  "Zanol2014",     # 74 tips, gap +2
  "Conrad2008"     # 64 tips, gap +1-2
)

cat("=== T-249: Gap dataset comparison (current cpp-search HEAD) ===\n")
cat("Start:", format(Sys.time()), "\n\n")

# Round 2 used clean_inapplicable() (=Fitch mode): inapplicable treated as
# missing for fair comparison with TNT, which always uses Fitch scoring.
results <- run_comparison(
  dataset_names = GAP_DATASETS,
  weightings = "EW",
  timeout_s = 120,
  seeds = 1:3,
  hits = 5L,
  reps = 50L,
  use_fitch = TRUE
)

cat("\n=== Done:", format(Sys.time()), "===\n\n")

outfile <- paste0("dev/benchmarks/t249_gap_",
                  format(Sys.time(), "%Y%m%d_%H%M"), ".csv")
write.csv(results, outfile, row.names = FALSE)
cat("Results saved to:", outfile, "\n")

# Load round 2 data for comparison
r2 <- read.csv("dev/benchmarks/round2_hard.csv")

# Summary: compare best scores per dataset
agg_new <- aggregate(cbind(tnt_score, ts_score) ~ dataset, data = results, FUN = min)
agg_new$gap_new <- agg_new$ts_score - agg_new$tnt_score

agg_old <- aggregate(cbind(tnt_score, ts_score) ~ dataset, data = r2, FUN = min)
agg_old$gap_old <- agg_old$ts_score - agg_old$tnt_score
agg_old <- agg_old[agg_old$dataset %in% GAP_DATASETS, ]

merged <- merge(agg_new[, c("dataset", "gap_new", "ts_score")],
                agg_old[, c("dataset", "gap_old")],
                by = "dataset", all.x = TRUE)
merged$improvement <- merged$gap_old - merged$gap_new
cat("\n--- Gap comparison: round 2 vs round 3 ---\n")
print(merged, row.names = FALSE)
