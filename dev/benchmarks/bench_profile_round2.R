# Profiling round 2: Fresh baselines and detailed phase analysis
# Agent F, S-PROF, 2026-03-17
#
# Run via: Rscript -e "library(TreeSearch, lib.loc='.agent-f'); source('dev/benchmarks/bench_profile_round2.R')"

library(TreeSearch, lib.loc = ".agent-f")
library(TreeTools)

# Representative datasets spanning the size range
DATASETS <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")

prepare <- function(name) {
  ds <- TreeSearch::inapplicable.phyData[[name]]
  at <- attributes(ds)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(ds, use.names = FALSE),
                      nrow = length(ds), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(ds)
  )
}

# ---- Section 1: End-to-end with timings attribute ----

cat("=== Section 1: End-to-end with per-phase timings ===\n\n")

for (nm in DATASETS) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))

  # 3 runs, take medians
  timings_list <- list()
  wall_times <- numeric(3)

  for (run in 1:3) {
    set.seed(7300 + run)
    t0 <- proc.time()
    result <- TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = 5L,
      targetHits = 3L,
      ratchetCycles = 5L,
      driftCycles = 5L,
      xssRounds = 1L,
      rssRounds = 1L,
      cssRounds = 1L,
      cssPartitions = 3L,
      xssPartitions = 3L,
      fuseInterval = 5L,
      maxSeconds = 60,
      verbosity = 0L,
      nThreads = 1L
    )
    elapsed <- (proc.time() - t0)[3]
    wall_times[run] <- elapsed
    timings_list[[run]] <- result$timings
    cat(sprintf("  Run %d: %.3fs wall, score=%.0f, reps=%d\n",
                run, elapsed, result$best_score, result$replicates))
  }

  # Median wall time
  med_wall <- median(wall_times)
  # Median per-phase (element-wise)
  med_timings <- sapply(names(timings_list[[1]]), function(ph) {
    median(sapply(timings_list, function(t) t[[ph]]))
  })
  cpp_total <- sum(med_timings)
  r_overhead <- med_wall * 1000 - cpp_total

  cat(sprintf("\n  Median wall: %.3fs\n", med_wall))
  cat("  Per-phase (median ms):\n")
  for (ph in names(med_timings)) {
    pct <- if (cpp_total > 0) 100 * med_timings[[ph]] / cpp_total else 0
    cat(sprintf("    %-12s %8.1f ms  (%4.1f%%)\n", ph, med_timings[[ph]], pct))
  }
  cat(sprintf("    %-12s %8.1f ms  (C++ total)\n", "TOTAL", cpp_total))
  cat(sprintf("    %-12s %8.1f ms  (R overhead: %.1f%% of wall)\n\n",
              "R overhead", r_overhead, 100 * r_overhead / (med_wall * 1000)))
}

# ---- Section 2: IW comparison ----

cat("=== Section 2: IW mode comparison ===\n\n")

for (nm in c("Vinther2008", "Zhu2013")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips, IW k=10) ---\n", nm, ds$n_taxa))

  wall_times <- numeric(3)
  timings_list <- list()

  for (run in 1:3) {
    set.seed(7300 + run)
    t0 <- proc.time()
    result <- TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      concavity = 10.0,
      maxReplicates = 5L,
      targetHits = 3L,
      ratchetCycles = 5L,
      driftCycles = 5L,
      xssRounds = 1L,
      rssRounds = 1L,
      cssRounds = 1L,
      cssPartitions = 3L,
      xssPartitions = 3L,
      fuseInterval = 5L,
      maxSeconds = 60,
      verbosity = 0L,
      nThreads = 1L
    )
    elapsed <- (proc.time() - t0)[3]
    wall_times[run] <- elapsed
    timings_list[[run]] <- result$timings
    cat(sprintf("  Run %d: %.3fs wall, score=%.2f, reps=%d\n",
                run, elapsed, result$best_score, result$replicates))
  }

  med_wall <- median(wall_times)
  med_timings <- sapply(names(timings_list[[1]]), function(ph) {
    median(sapply(timings_list, function(t) t[[ph]]))
  })
  cpp_total <- sum(med_timings)

  cat(sprintf("\n  Median wall: %.3fs\n", med_wall))
  cat("  Per-phase (median ms):\n")
  for (ph in names(med_timings)) {
    pct <- if (cpp_total > 0) 100 * med_timings[[ph]] / cpp_total else 0
    cat(sprintf("    %-12s %8.1f ms  (%4.1f%%)\n", ph, med_timings[[ph]], pct))
  }
  cat(sprintf("    %-12s %8.1f ms  (C++ total)\n\n", "TOTAL", cpp_total))
}

# ---- Section 3: Scaling test ----

cat("=== Section 3: Scaling — single TBR pass timing ===\n\n")

for (nm in DATASETS) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))

  # Single replicate, no sectorial/ratchet/drift — just Wagner+TBR
  wall_times <- numeric(3)
  timings_list <- list()

  for (run in 1:3) {
    set.seed(7300 + run)
    t0 <- proc.time()
    result <- TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = 1L,
      targetHits = 1L,
      ratchetCycles = 0L,
      driftCycles = 0L,
      xssRounds = 0L,
      rssRounds = 0L,
      cssRounds = 0L,
      fuseInterval = 0L,
      maxSeconds = 60,
      verbosity = 0L,
      nThreads = 1L
    )
    elapsed <- (proc.time() - t0)[3]
    wall_times[run] <- elapsed
    timings_list[[run]] <- result$timings
  }

  med_wall <- median(wall_times)
  med_timings <- sapply(names(timings_list[[1]]), function(ph) {
    median(sapply(timings_list, function(t) t[[ph]]))
  })

  cat(sprintf("  Wagner:  %6.1f ms\n", med_timings[["wagner"]]))
  cat(sprintf("  TBR:     %6.1f ms\n", med_timings[["tbr"]]))
  cat(sprintf("  Wall:    %6.1f ms\n", med_wall * 1000))
  cat(sprintf("  R ovhd:  %6.1f ms\n\n", med_wall * 1000 - sum(med_timings)))
}

cat("=== Profiling complete ===\n")
