# Profiling round 2b: Drift/ratchet deep dive + scaling
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

# ---- Section 3: Drift cycle count sensitivity ----

cat("=== Section 3: Drift cycle count sensitivity ===\n\n")

# How does drift time scale with cycle count?
# The question: are we doing too many drift cycles for the benefit?

for (nm in c("Zhu2013", "Dikow2009")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "dCycles", "drift_ms", "total_ms", "score", "reps"))

  for (dc in c(0L, 1L, 2L, 3L, 5L, 10L)) {
    scores <- numeric(3)
    drift_ms <- numeric(3)
    total_ms <- numeric(3)
    reps <- numeric(3)

    for (run in 1:3) {
      set.seed(7300 + run)
      t0 <- proc.time()
      result <- TreeSearch:::ts_driven_search(
        ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxReplicates = 3L,
        targetHits = 3L,
        ratchetCycles = 5L,
        driftCycles = dc,
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
      scores[run] <- result$best_score
      total_ms[run] <- elapsed * 1000
      drift_ms[run] <- result$timings[["drift_ms"]]
      reps[run] <- result$replicates
    }

    cat(sprintf("  %-8d %8.0f %8.0f %8.0f %8.0f\n",
                dc, median(drift_ms), median(total_ms),
                median(scores), median(reps)))
  }
  cat("\n")
}

# ---- Section 4: Ratchet cycle count sensitivity ----

cat("=== Section 4: Ratchet cycle count sensitivity ===\n\n")

for (nm in c("Zhu2013", "Dikow2009")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))
  cat(sprintf("  %-8s %8s %8s %8s %8s\n", "rCycles", "ratch_ms", "total_ms", "score", "reps"))

  for (rc in c(0L, 1L, 2L, 3L, 5L, 10L)) {
    scores <- numeric(3)
    ratch_ms <- numeric(3)
    total_ms <- numeric(3)
    reps <- numeric(3)

    for (run in 1:3) {
      set.seed(7300 + run)
      t0 <- proc.time()
      result <- TreeSearch:::ts_driven_search(
        ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxReplicates = 3L,
        targetHits = 3L,
        ratchetCycles = rc,
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
      elapsed <- (proc.time() - t0)[3]
      scores[run] <- result$best_score
      total_ms[run] <- elapsed * 1000
      ratch_ms[run] <- result$timings[["ratchet_ms"]]
      reps[run] <- result$replicates
    }

    cat(sprintf("  %-8d %8.0f %8.0f %8.0f %8.0f\n",
                rc, median(ratch_ms), median(total_ms),
                median(scores), median(reps)))
  }
  cat("\n")
}

# ---- Section 5: CSS effectiveness ----

cat("=== Section 5: CSS vs no CSS ===\n\n")

for (nm in c("Zhu2013", "Dikow2009")) {
  ds <- prepare(nm)
  cat(sprintf("--- %s (%d tips) ---\n", nm, ds$n_taxa))

  for (css in c(0L, 1L, 2L)) {
    scores <- numeric(3)
    css_ms <- numeric(3)
    total_ms <- numeric(3)

    for (run in 1:3) {
      set.seed(7300 + run)
      t0 <- proc.time()
      result <- TreeSearch:::ts_driven_search(
        ds$contrast, ds$tip_data, ds$weight, ds$levels,
        maxReplicates = 3L,
        targetHits = 3L,
        ratchetCycles = 5L,
        driftCycles = 5L,
        xssRounds = 1L,
        rssRounds = 1L,
        cssRounds = css,
        cssPartitions = 3L,
        xssPartitions = 3L,
        fuseInterval = 5L,
        maxSeconds = 120,
        verbosity = 0L,
        nThreads = 1L
      )
      elapsed <- (proc.time() - t0)[3]
      scores[run] <- result$best_score
      css_ms[run] <- result$timings[["css_ms"]]
      total_ms[run] <- elapsed * 1000
    }

    cat(sprintf("  cssRounds=%d: css_ms=%6.0f total_ms=%6.0f score=%6.0f\n",
                css, median(css_ms), median(total_ms), median(scores)))
  }
  cat("\n")
}

# ---- Section 6: Wagner + TBR-only (no perturbation) ----

cat("=== Section 6: Wagner + TBR only (scaling) ===\n\n")

DATASETS <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")

for (nm in DATASETS) {
  ds <- prepare(nm)
  wall_times <- numeric(5)
  tbr_ms <- numeric(5)
  wagner_ms <- numeric(5)

  for (run in 1:5) {
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
    tbr_ms[run] <- result$timings[["tbr_ms"]]
    wagner_ms[run] <- result$timings[["wagner_ms"]]
  }

  cat(sprintf("  %s (%2d tips): Wagner=%5.1f ms, TBR=%7.1f ms, Wall=%7.1f ms\n",
              nm, ds$n_taxa,
              median(wagner_ms), median(tbr_ms), median(wall_times) * 1000))
}

cat("\n=== Profiling complete ===\n")
