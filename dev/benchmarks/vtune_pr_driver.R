#!/usr/bin/env Rscript
# VTune driver: prune-reinsert hotspot profiling
#
# Exercises prune_reinsert_search heavily on Zhu2013 (75t) and Dikow2009 (88t).
# Target: ~30-60s of CPU time in the PR hot path.
#
# Usage:
#   Rscript dev/benchmarks/vtune_pr_driver.R
#   vtune -collect hotspots -result-dir vtune-pr-out -- Rscript dev/benchmarks/vtune_pr_driver.R

.libPaths(c(".vtune-lib", .libPaths()))
library(TreeSearch)
library(TreeTools)

cat("TreeSearch:", as.character(packageVersion("TreeSearch")), "\n")
cat("Dataset: Zhu2013 (75t) + Dikow2009 (88t)\n\n")

# Use both datasets for a more representative profile
datasets <- list(
  Zhu2013    = inapplicable.phyData[["Zhu2013"]],
  Dikow2009  = inapplicable.phyData[["Dikow2009"]]
)

t0 <- proc.time()

for (ds_name in names(datasets)) {
  ds <- datasets[[ds_name]]
  cat(sprintf("Running %s ...\n", ds_name))

  # Maximise PR time share: high cycle count, no ratchet/drift/NNI-perturb,
  # enough time for ~100+ replicates worth of PR work.
  set.seed(7531)
  MaximizeParsimony(
    ds,
    maxSeconds          = 40L,
    strategy            = "auto",
    pruneReinsertCycles = 5L,
    pruneReinsertDrop   = 0.10,
    driftCycles         = 0L,
    nniPerturbCycles    = 0L,
    verbosity           = 0L,
    nThreads            = 1L
  )

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("  done (%.1fs elapsed)\n", elapsed))
}

cat(sprintf("\nTotal: %.1fs\n", (proc.time() - t0)[3]))
