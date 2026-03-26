# Quick smoke test for the agent-F build + TNT comparison pipeline
agent_lib <- normalizePath("../.builds/TreeSearch-F", winslash = "/")
.libPaths(c(agent_lib, .libPaths()))
library(TreeSearch.F)
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))
library(TreeTools)

# Source helper files, skipping library() and nested source() calls
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

cat("Sourcing bench_datasets.R...\n")
source_filtered("dev/benchmarks/bench_datasets.R")
cat("Sourcing bench_tnt_compare.R...\n")
source_filtered("dev/benchmarks/bench_tnt_compare.R")

cat("Smoke test: Vinther2008, EW, 5s, seed=1\n")
res <- run_comparison(
  dataset_names = "Vinther2008",
  weightings = "EW",
  timeout_s = 5,
  seeds = 1L,
  hits = 3L,
  reps = 10L,
  use_fitch = TRUE
)
print(res[, c("dataset", "tnt_score", "ts_score", "tnt_wall_s", "ts_wall_s")])
