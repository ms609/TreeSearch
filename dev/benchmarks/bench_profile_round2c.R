# Profiling round 2c: Parallel scaling + quality impact of reduced cycles
# Agent F, S-PROF, 2026-03-17

library(TreeSearch, lib.loc = ".agent-f")
library(TreeTools)

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

# ---- Section 7: Quality impact with more statistical power ----

cat("=== Section 7: Drift/ratchet tuning — quality impact (10 seeds) ===\n\n")

run_config <- function(ds, drift, ratchet, seed) {
  set.seed(seed)
  t0 <- proc.time()
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = 5L,
    targetHits = 3L,
    ratchetCycles = ratchet,
    driftCycles = drift,
    xssRounds = 1L,
    rssRounds = 1L,
    cssRounds = 1L,
    cssPartitions = 3L,
    xssPartitions = 3L,
    fuseInterval = 5L,
    maxSeconds = 120,
    verbosity = 0L,
    nThreads = 1L
  )
  elapsed <- (proc.time() - t0)[3]
  c(score = unname(result$best_score), time = unname(elapsed), reps = unname(result$replicates))
}

configs <- list(
  "d5_r5"  = c(drift = 5, ratchet = 5),   # current default
  "d2_r2"  = c(drift = 2, ratchet = 2),   # reduced
  "d2_r5"  = c(drift = 2, ratchet = 5),   # drift only reduced
  "d5_r2"  = c(drift = 5, ratchet = 2),   # ratchet only reduced
  "d0_r5"  = c(drift = 0, ratchet = 5),   # no drift
  "d5_r0"  = c(drift = 5, ratchet = 0)    # no ratchet
)

seeds <- 7301:7310

for (nm in c("Zhu2013", "Dikow2009")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips, 10 seeds) ---\n", nm, ds$n_taxa))
  cat(sprintf("  %-8s %8s %8s %8s %8s %8s\n",
              "config", "med_scr", "mean_scr", "min_scr", "med_time", "mean_t"))

  for (cfg_name in names(configs)) {
    cfg <- configs[[cfg_name]]
    sc <- numeric(length(seeds))
    tm <- numeric(length(seeds))
    for (i in seq_along(seeds)) {
      r <- run_config(ds, cfg[["drift"]], cfg[["ratchet"]], seeds[i])
      sc[i] <- r[["score"]]
      tm[i] <- r[["time"]]
    }

    cat(sprintf("  %-8s %8.0f %8.1f %8.0f %8.1f %8.1f\n",
                cfg_name, median(sc), mean(sc), min(sc),
                median(tm), mean(tm)))
  }
  cat("\n")
}

# ---- Section 8: Parallel scaling ----

cat("=== Section 8: Parallel scaling ===\n\n")

for (nm in c("Zhu2013")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))
  cat(sprintf("  %-10s %8s %8s %8s\n", "nThreads", "time_ms", "score", "reps"))

  for (nt in c(1L, 2L)) {
    times <- numeric(3)
    scores <- numeric(3)
    reps <- numeric(3)

    for (run in 1:3) {
      set.seed(7300 + run)
      t0 <- proc.time()
      result <- TreeSearch:::ts_driven_search(
        ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxReplicates = 5L,
        targetHits = 5L,
        ratchetCycles = 5L,
        driftCycles = 5L,
        xssRounds = 1L,
        rssRounds = 1L,
        cssRounds = 1L,
        cssPartitions = 3L,
        xssPartitions = 3L,
        fuseInterval = 5L,
        maxSeconds = 120,
        verbosity = 0L,
        nThreads = nt
      )
      elapsed <- (proc.time() - t0)[3]
      times[run] <- elapsed
      scores[run] <- result$best_score
      reps[run] <- result$replicates
    }

    cat(sprintf("  %-10d %8.0f %8.0f %8.0f\n",
                nt, median(times) * 1000, median(scores), median(reps)))
  }
  cat("\n")
}

# ---- Section 9: Per-replicate cost breakdown ----

cat("=== Section 9: Per-replicate cost (ms/rep) ===\n\n")

DATASETS <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")
cat(sprintf("  %-15s %4s %8s %8s %8s %8s %8s %8s %8s\n",
            "dataset", "tips", "wagner", "tbr", "sect", "ratch", "drift", "fTBR", "TOTAL"))

for (nm in DATASETS) {
  ds <- prepare(nm)

  all_timings <- list()
  all_reps <- numeric(3)

  for (run in 1:3) {
    set.seed(7300 + run)
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
      maxSeconds = 120,
      verbosity = 0L,
      nThreads = 1L
    )
    all_timings[[run]] <- result$timings
    all_reps[run] <- result$replicates
  }

  med_reps <- median(all_reps)
  med_t <- sapply(names(all_timings[[1]]), function(ph) {
    median(sapply(all_timings, function(t) t[[ph]]))
  })

  # Per-rep average
  pr <- med_t / med_reps
  sect <- pr[["xss_ms"]] + pr[["rss_ms"]] + pr[["css_ms"]]
  total <- sum(pr)

  cat(sprintf("  %-15s %4d %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f\n",
              nm, ds$n_taxa,
              pr[["wagner_ms"]], pr[["tbr_ms"]], sect,
              pr[["ratchet_ms"]], pr[["drift_ms"]],
              pr[["final_tbr_ms"]], total))
}

cat("\n=== Profiling complete ===\n")
