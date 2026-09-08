#!/usr/bin/env Rscript
# Benchmark: perturbStopFactor effectiveness across dataset sizes
#
# Compares search convergence with different perturbStopFactor settings.
# For each dataset, runs MaximizeParsimony with:
#   - Baseline (perturbStopFactor = 0, i.e. disabled)
#   - perturbStopFactor = 2
#   - perturbStopFactor = 5
#
# Measures: elapsed time, best score, replicates completed.

.libPaths(c(
  "C:/Users/pjjg18/GitHub/.builds/TreeSearch-Z",
  .libPaths()
))
library(TreeSearch.Z)
if (is.null(.Internal(getRegisteredNamespace("TreeSearch"))))
  .Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.Z")))
library(TreeTools)

# Select datasets across the size spectrum.
# Use the inst/datasets (inapplicable.phyData) for small/medium,
# plus morphobank datasets from neotrans for large/XL.
neotrans_dir <- system.file("matrices", package = "neotrans")

load_dataset <- function(name, source = "inapplicable") {
  if (source == "inapplicable") {
    return(TreeSearch::inapplicable.phyData[[name]])
  } else {
    path <- file.path(neotrans_dir, paste0(name, ".nex"))
    return(suppressWarnings(TreeTools::ReadAsPhyDat(path)))
  }
}

# Dataset selection: cover small (20-40), medium (41-80), large (81-150),
# XL (150+). Focus on medium-to-XL where the feature is most relevant.
datasets_spec <- list(
  # Small — expect quick convergence, perturb-stop shouldn't matter
  list(name = "Vinther2008", source = "inapplicable", ntip = 23),
  list(name = "Aria2015", source = "inapplicable", ntip = 35),

  # Medium — starts to get interesting
  list(name = "Griswold1999", source = "inapplicable", ntip = 43),
  list(name = "Eklund2004", source = "inapplicable", ntip = 54),

  # Medium-large — key range
  list(name = "Agnarsson2004", source = "inapplicable", ntip = 62),
  list(name = "Zhu2013", source = "inapplicable", ntip = 75),
  list(name = "Dikow2009", source = "inapplicable", ntip = 88),

  # Large — from morphobank/neotrans
  list(name = "project2086", source = "neotrans", ntip = 91),
  list(name = "project2769", source = "neotrans", ntip = 102),
  list(name = "project1013", source = "neotrans", ntip = 112),
  list(name = "project2286", source = "neotrans", ntip = 134),

  # XL
  list(name = "project1024", source = "neotrans", ntip = 163),
  list(name = "project2477", source = "neotrans", ntip = 213)
)

# perturbStopFactor values to test (0 = disabled = baseline)
psf_values <- c(0L, 2L, 5L)

# Per-dataset time budget: scale with tip count
# Small: 15s, Medium: 30s, Large: 60s, XL: 90s
time_budget <- function(ntip) {
  if (ntip <= 40) 15
  else if (ntip <= 80) 30
  else if (ntip <= 150) 60
  else 90
}

# Number of reps per condition
n_reps <- 2L

set.seed(7418)

results <- data.frame(
  dataset = character(),
  ntip = integer(),
  nchar = integer(),
  psf = integer(),
  rep = integer(),
  elapsed_s = numeric(),
  best_score = numeric(),
  n_replicates = integer(),
  stringsAsFactors = FALSE
)

cat("=== Perturbation-Stop Benchmark ===\n")
cat(sprintf("Datasets: %d, PSF values: %s, Reps: %d\n",
            length(datasets_spec),
            paste(psf_values, collapse = "/"),
            n_reps))

for (ds_spec in datasets_spec) {
  cat(sprintf("\n--- %s (%d tips) ---\n", ds_spec$name, ds_spec$ntip))

  dataset <- tryCatch(
    load_dataset(ds_spec$name, ds_spec$source),
    error = function(e) {
      cat("  SKIP: ", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(dataset)) next

  actual_ntip <- length(dataset)
  actual_nchar <- sum(attr(dataset, "weight"))
  budget <- time_budget(actual_ntip)

  cat(sprintf("  Actual: %d tips, %d chars, budget: %ds\n",
              actual_ntip, actual_nchar, budget))

  for (psf in psf_values) {
    for (r in seq_len(n_reps)) {
      seed <- sample.int(10000, 1)
      set.seed(seed)

      ctrl <- SearchControl(perturbStopFactor = psf)

      t0 <- proc.time()["elapsed"]
      res <- tryCatch(
        MaximizeParsimony(
          dataset,
          control = ctrl,
          maxSeconds = budget,
          maxReplicates = 500L,
          targetHits = max(10L, as.integer(actual_ntip / 5)),
          verbosity = 0L,
          nThreads = 2L
        ),
        error = function(e) {
          cat(sprintf("  ERROR (psf=%d, rep=%d): %s\n",
                      psf, r, conditionMessage(e)))
          NULL
        }
      )
      elapsed <- proc.time()["elapsed"] - t0

      if (!is.null(res)) {
        best <- attr(res, "score")
        if (is.null(best)) best <- TreeLength(res[[1]], dataset)
        n_reps_done <- attr(res, "replicates")
        if (is.null(n_reps_done)) n_reps_done <- NA_integer_

        results <- rbind(results, data.frame(
          dataset = ds_spec$name,
          ntip = actual_ntip,
          nchar = actual_nchar,
          psf = psf,
          rep = r,
          elapsed_s = round(elapsed, 2),
          best_score = best,
          n_replicates = n_reps_done,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("  psf=%d rep=%d: %.1fs, score=%.1f, reps=%s\n",
                    psf, r, elapsed,
                    best,
                    ifelse(is.na(n_reps_done), "?", as.character(n_reps_done))))
      }
    }
  }
}

cat("\n\n=== Summary ===\n")

# Aggregate by dataset x psf
agg <- aggregate(
  cbind(elapsed_s, best_score) ~ dataset + ntip + nchar + psf,
  data = results,
  FUN = mean
)
agg <- agg[order(agg$ntip, agg$psf), ]

# Reshape for comparison
baseline <- agg[agg$psf == 0, c("dataset", "ntip", "nchar",
                                 "elapsed_s", "best_score")]
names(baseline)[4:5] <- c("time_base", "score_base")

for (p in psf_values[psf_values > 0]) {
  psf_rows <- agg[agg$psf == p, c("dataset", "elapsed_s", "best_score")]
  names(psf_rows)[2:3] <- paste0(c("time_psf", "score_psf"), p)
  baseline <- merge(baseline, psf_rows, by = "dataset", all.x = TRUE)
}

baseline <- baseline[order(baseline$ntip), ]
cat("\n")
print(baseline, row.names = FALSE)

# Save
out_path <- "dev/benchmarks/results_perturb_stop.csv"
write.csv(results, out_path, row.names = FALSE)
cat(sprintf("\nRaw results saved to: %s\n", out_path))
