#!/usr/bin/env Rscript
# Benchmark v3: isolate PSF by disabling targetHits
#
# With targetHits effectively disabled (set to 999), the only
# stopping criteria are: maxReplicates, maxSeconds, or PSF.
# This shows whether PSF would ever fire as a pure convergence signal.

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
  if (source == "inapplicable") TreeSearch::inapplicable.phyData[[name]]
  else suppressWarnings(TreeTools::ReadAsPhyDat(
    file.path(neotrans_dir, paste0(name, ".nex"))))
}

# Focus on medium datasets where per-rep cost is low enough to get many reps
datasets_spec <- list(
  list(name = "Griswold1999", source = "inapplicable", ntip = 43),
  list(name = "Eklund2004", source = "inapplicable", ntip = 54),
  list(name = "Agnarsson2004", source = "inapplicable", ntip = 62),
  list(name = "Zhu2013", source = "inapplicable", ntip = 75),
  list(name = "Dikow2009", source = "inapplicable", ntip = 88)
)

psf_values <- c(0L, 2L, 5L)
n_reps <- 3L
set.seed(8321)

cat("=== PSF Isolation Test ===\n")
cat("targetHits=999 (disabled), maxReplicates=500, maxSeconds=300\n\n")

results <- list()

for (ds_spec in datasets_spec) {
  dataset <- load_dataset(ds_spec$name, ds_spec$source)
  if (is.null(dataset)) next
  ntip <- length(dataset)
  nchar <- sum(attr(dataset, "weight"))

  cat(sprintf("\n--- %s (%d tips, %d chars) ---\n", ds_spec$name, ntip, nchar))
  cat(sprintf("  PSF limits: psf=2 → %d reps, psf=5 → %d reps\n",
              ntip * 2L, ntip * 5L))

  for (psf in psf_values) {
    for (r in seq_len(n_reps)) {
      set.seed(sample.int(10000, 1))
      ctrl <- SearchControl(perturbStopFactor = psf)

      t0 <- proc.time()["elapsed"]
      res <- tryCatch(
        MaximizeParsimony(
          dataset, control = ctrl,
          maxSeconds = 300, maxReplicates = 500L,
          targetHits = 999L,
          verbosity = 0L, nThreads = 2L
        ),
        error = function(e) { cat("  ERR:", conditionMessage(e), "\n"); NULL }
      )
      elapsed <- proc.time()["elapsed"] - t0

      if (!is.null(res)) {
        best <- attr(res, "score")
        if (is.null(best)) best <- TreeLength(res[[1]], dataset)
        nrep <- attr(res, "replicates")
        if (is.null(nrep)) nrep <- NA_integer_

        stop_reason <- if (!is.na(nrep) && nrep >= 500) "maxReps"
          else if (elapsed >= 285) "time"
          else "PSF/converged"

        results[[length(results) + 1]] <- data.frame(
          dataset = ds_spec$name, ntip = ntip, nchar = nchar,
          psf = psf, rep = r, elapsed_s = round(elapsed, 1),
          best_score = best, n_replicates = nrep, stop = stop_reason,
          stringsAsFactors = FALSE
        )

        cat(sprintf("  psf=%d rep=%d: %.0fs, score=%.0f, reps=%s [%s]\n",
                    psf, r, elapsed, best,
                    ifelse(is.na(nrep), "?", as.character(nrep)),
                    stop_reason))
      }
    }
  }
}

results_df <- do.call(rbind, results)
cat("\n\n=== Did PSF ever fire? ===\n")
psf_fired <- results_df[results_df$psf > 0 & results_df$stop == "PSF/converged", ]
if (nrow(psf_fired) > 0) {
  cat("YES — PSF fired in these cases:\n")
  print(psf_fired, row.names = FALSE)
} else {
  cat("NO — PSF never fired. All runs ended by maxReplicates or time.\n")
}

cat("\n=== Full results ===\n")
print(results_df[order(results_df$ntip, results_df$psf, results_df$rep), ],
      row.names = FALSE)

write.csv(results_df, "dev/benchmarks/results_perturb_stop_v3.csv", row.names = FALSE)
