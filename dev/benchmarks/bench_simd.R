# Phase 3E SIMD benchmark: measure TBR search performance.
#
# This benchmark compares SIMD-enabled TBR performance across dataset sizes.
# Since SIMD is compiled in (no runtime toggle), we measure absolute timings
# and per-candidate costs to verify the Phase 3D profiling baseline is met
# or improved.
#
# Usage: Rscript dev/benchmarks/bench_simd.R

library(TreeSearch)
library(TreeTools)

cat("Phase 3E SIMD Benchmark\n")
cat("=======================\n\n")

# Helper: run TBR search and measure time
bench_tbr <- function(dataset, n_reps = 3, label = "") {
  ds <- list(
    contrast = attr(dataset, "contrast"),
    tip_data = t(vapply(dataset, I, dataset[[1]])),
    weight   = attr(dataset, "weight"),
    levels   = attr(dataset, "levels")
  )
  n_tip <- length(dataset)
  tree <- Preorder(PectinateTree(dataset))

  times <- vapply(seq_len(n_reps), function(i) {
    set.seed(4200 + i)
    t0 <- proc.time()
    TreeSearch:::ts_tbr_search(
      tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxHits = 5L
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]
    elapsed
  }, numeric(1))

  med <- median(times)
  cat(sprintf("  %-30s  tips=%3d  median=%.3fs  (%.3f, %.3f, %.3f)\n",
              label, n_tip, med, times[1], times[2], times[3]))
  data.frame(label = label, n_tip = n_tip, median_s = med,
             stringsAsFactors = FALSE)
}

# Helper: run driven search and measure time
bench_driven <- function(dataset, n_reps = 3, label = "") {
  n_tip <- length(dataset)
  ds <- list(
    contrast = attr(dataset, "contrast"),
    tip_data = t(vapply(dataset, I, dataset[[1]])),
    weight   = attr(dataset, "weight"),
    levels   = attr(dataset, "levels")
  )

  times <- vapply(seq_len(n_reps), function(i) {
    set.seed(4200 + i)
    t0 <- proc.time()
    TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = 2L, targetHits = 2L,
      ratchetCycles = 2L, driftCycles = 0L,
      xssPartitions = 2L, rssRounds = 0L,
      cssRounds = 0L, cssPartitions = 2L,
      fuseInterval = 0L, poolMaxSize = 2L,
      maxSeconds = 30, verbosity = 0L
    )
    elapsed <- (proc.time() - t0)[["elapsed"]]
    elapsed
  }, numeric(1))

  med <- median(times)
  cat(sprintf("  %-30s  tips=%3d  median=%.3fs  (%.3f, %.3f, %.3f)\n",
              label, n_tip, med, times[1], times[2], times[3]))
  data.frame(label = label, n_tip = n_tip, median_s = med,
             stringsAsFactors = FALSE)
}

# ---- TBR benchmarks ----
cat("TBR search (5 hits to best):\n")
results_tbr <- list()

for (ds_name in c("Vinther2008", "Agnarsson2004", "Wills2012",
                   "Aria2015", "Zhu2013")) {
  dataset <- inapplicable.phyData[[ds_name]]
  results_tbr[[ds_name]] <- bench_tbr(dataset, label = ds_name)
}

# DNA dataset
suppressWarnings(data("Laurasiatherian", package = "phangorn"))
results_tbr[["Laurasiatherian"]] <- bench_tbr(Laurasiatherian,
                                               label = "Laurasiatherian (DNA)")

cat("\nDriven search (2 replicates, 30s timeout):\n")
results_driven <- list()

for (ds_name in c("Vinther2008", "Agnarsson2004", "Zhu2013")) {
  dataset <- inapplicable.phyData[[ds_name]]
  results_driven[[ds_name]] <- bench_driven(dataset, label = ds_name)
}

# Phase benchmark diagnostic (if available)
cat("\nTBR phase timing (Phase 3D diagnostic):\n")
for (ds_name in c("Vinther2008", "Zhu2013")) {
  dataset <- inapplicable.phyData[[ds_name]]
  ds <- list(
    contrast = attr(dataset, "contrast"),
    tip_data = t(vapply(dataset, I, dataset[[1]])),
    weight   = attr(dataset, "weight"),
    levels   = attr(dataset, "levels")
  )
  tree <- Preorder(PectinateTree(dataset))
  set.seed(7777)
  ph <- tryCatch(
    TreeSearch:::ts_bench_tbr_phases(
      tree$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxHits = 3L
    ),
    error = function(e) NULL
  )
  if (!is.null(ph)) {
    cat(sprintf("  %s: clip=%.1fms indirect=%.1fms verify=%.1fms total=%.1fms\n",
                ds_name,
                ph$clip_us / 1000, ph$indirect_us / 1000,
                ph$verify_us / 1000, ph$total_us / 1000))
  }
}

cat("\nBenchmark complete.\n")
