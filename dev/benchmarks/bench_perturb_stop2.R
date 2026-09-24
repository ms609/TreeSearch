#!/usr/bin/env Rscript
# Benchmark v2: perturb-stop with generous time, replicate-limited
#
# Goal: see if perturbStopFactor can terminate searches early
# (before maxReplicates) and whether the scores it finds are equivalent.
#
# Key change from v1: use maxSeconds = 600 (generous) so
# replicate-based criteria can fire. Cap maxReplicates = 200.

.libPaths(c(
  "C:/Users/pjjg18/GitHub/.builds/TreeSearch-Z",
  .libPaths()
))
library(TreeSearch.Z)
if (is.null(.Internal(getRegisteredNamespace("TreeSearch"))))
  .Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.Z")))
library(TreeTools)

neotrans_dir <- system.file("matrices", package = "neotrans")

load_dataset <- function(name, source = "inapplicable") {
  if (source == "inapplicable") {
    return(TreeSearch::inapplicable.phyData[[name]])
  } else {
    path <- file.path(neotrans_dir, paste0(name, ".nex"))
    return(suppressWarnings(TreeTools::ReadAsPhyDat(path)))
  }
}

# Focus on medium-to-large datasets where convergence behavior matters.
# Include small ones as controls.
datasets_spec <- list(
  # Small — should converge very quickly regardless
  list(name = "Vinther2008", source = "inapplicable", ntip = 23),
  list(name = "Aria2015", source = "inapplicable", ntip = 35),

  # Medium — may or may not converge
  list(name = "Griswold1999", source = "inapplicable", ntip = 43),
  list(name = "Eklund2004", source = "inapplicable", ntip = 54),
  list(name = "Agnarsson2004", source = "inapplicable", ntip = 62),
  list(name = "Zhu2013", source = "inapplicable", ntip = 75),
  list(name = "Dikow2009", source = "inapplicable", ntip = 88),

  # Large — from morphobank/neotrans
  list(name = "project2086", source = "neotrans", ntip = 91),
  list(name = "project2769", source = "neotrans", ntip = 102),
  list(name = "project1013", source = "neotrans", ntip = 112)
)

psf_values <- c(0L, 2L, 5L)

# Per-dataset: generous time, moderate maxReplicates to let
# convergence criteria fire. The question is whether PSF terminates
# before targetHits, and at what score.
max_reps_by_size <- function(ntip) {
  # Enough replicates that targetHits should fire for easy datasets,
  # but hard datasets won't hit targetHits within budget.
  if (ntip <= 40) 100L
  else if (ntip <= 80) 150L
  else 200L
}

max_seconds_by_size <- function(ntip) {
  # Generous: 5x what the first benchmark showed was needed
  if (ntip <= 40) 30
  else if (ntip <= 80) 120
  else 300
}

n_reps <- 2L
set.seed(4193)

results <- data.frame(
  dataset = character(), ntip = integer(), nchar = integer(),
  psf = integer(), rep = integer(), elapsed_s = numeric(),
  best_score = numeric(), n_replicates = integer(),
  stringsAsFactors = FALSE
)

cat("=== Perturbation-Stop Benchmark v2 ===\n")
cat("Focus: do stopping criteria fire before time limit?\n\n")

for (ds_spec in datasets_spec) {
  cat(sprintf("\n--- %s (%d tips) ---\n", ds_spec$name, ds_spec$ntip))
  dataset <- tryCatch(load_dataset(ds_spec$name, ds_spec$source),
                      error = function(e) { cat("SKIP:", conditionMessage(e), "\n"); NULL })
  if (is.null(dataset)) next

  actual_ntip <- length(dataset)
  actual_nchar <- sum(attr(dataset, "weight"))
  max_reps <- max_reps_by_size(actual_ntip)
  max_secs <- max_seconds_by_size(actual_ntip)
  target_hits <- max(10L, as.integer(actual_ntip / 5))

  cat(sprintf("  %d tips, %d chars | maxReps=%d, maxSec=%d, targetHits=%d\n",
              actual_ntip, actual_nchar, max_reps, max_secs, target_hits))

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
          maxSeconds = max_secs,
          maxReplicates = max_reps,
          targetHits = target_hits,
          verbosity = 0L,
          nThreads = 2L
        ),
        error = function(e) {
          cat(sprintf("  ERROR (psf=%d, rep=%d): %s\n", psf, r, conditionMessage(e)))
          NULL
        }
      )
      elapsed <- proc.time()["elapsed"] - t0

      if (!is.null(res)) {
        best <- attr(res, "score")
        if (is.null(best)) best <- TreeLength(res[[1]], dataset)
        n_reps_done <- attr(res, "replicates")
        if (is.null(n_reps_done)) n_reps_done <- NA_integer_

        # Determine which criterion likely fired
        stop_reason <- "?"
        if (!is.na(n_reps_done)) {
          if (n_reps_done >= max_reps) stop_reason <- "maxReps"
          else if (elapsed >= max_secs * 0.95) stop_reason <- "time"
          else stop_reason <- "converged"
        }

        results <- rbind(results, data.frame(
          dataset = ds_spec$name, ntip = actual_ntip, nchar = actual_nchar,
          psf = psf, rep = r, elapsed_s = round(elapsed, 2),
          best_score = best, n_replicates = n_reps_done,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("  psf=%d rep=%d: %.1fs, score=%.0f, reps=%s [%s]\n",
                    psf, r, elapsed, best,
                    ifelse(is.na(n_reps_done), "?", as.character(n_reps_done)),
                    stop_reason))
      }
    }
  }
}

cat("\n\n=== Summary Table ===\n")
agg <- aggregate(
  cbind(elapsed_s, best_score, n_replicates) ~ dataset + ntip + nchar + psf,
  data = results, FUN = mean
)
agg <- agg[order(agg$ntip, agg$psf), ]

# Print nicely
for (ds in unique(agg$dataset)) {
  rows <- agg[agg$dataset == ds, ]
  cat(sprintf("\n%s (%d tips, %d chars):\n", ds, rows$ntip[1], rows$nchar[1]))
  for (i in seq_len(nrow(rows))) {
    cat(sprintf("  psf=%d: %.1fs, score=%.1f, reps=%.0f\n",
                rows$psf[i], rows$elapsed_s[i],
                rows$best_score[i], rows$n_replicates[i]))
  }
}

write.csv(results, "dev/benchmarks/results_perturb_stop_v2.csv", row.names = FALSE)
cat("\nSaved to dev/benchmarks/results_perturb_stop_v2.csv\n")
