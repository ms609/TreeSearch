# T-069: Stress test at 150–225 taxa
# Agent F, 2026-03-18
#
# Three large neotrans matrices: project175 (165t), project3763 (205t), syab07204 (225t)
# Goals:
#   1. Per-phase timing and phase distribution at large size
#   2. TBR pass micro-benchmark via ts_bench_tbr_phases
#   3. Pool behaviour (pool_size, replicates, fuse events)
#   4. Scaling exponent for indirect scoring vs smaller datasets
#
# Run via:
#   Rscript --vanilla -e "library(TreeSearch, lib.loc='.agent-f'); source('inst/benchmarks/bench_stress_large.R')"

library(TreeSearch, lib.loc = ".agent-f")
library(TreeTools)

NEOTRANS_DIR <- "../neotrans/inst/matrices"
MATRICES <- c("project175.nex", "project3763.nex", "syab07204.nex")

# ---- Helpers ----------------------------------------------------------------

load_nex <- function(file) {
  path <- file.path(NEOTRANS_DIR, file)
  ReadAsPhyDat(path)
}

prep_ds <- function(phyDat) {
  at <- attributes(phyDat)
  contrast <- at$contrast
  storage.mode(contrast) <- "double"
  tip_data <- matrix(unlist(phyDat, use.names = FALSE),
                     nrow = length(phyDat), byrow = TRUE)
  storage.mode(tip_data) <- "integer"
  weight   <- at$weight
  levels   <- at$levels
  # min_steps: number of non-zero contrast entries minus 1, clamped to 0
  min_steps <- pmax(apply(contrast, 2, function(x) sum(x > 0)) - 1L, 0L)
  list(contrast = contrast, tip_data = tip_data, weight = weight,
       levels = levels, min_steps = min_steps,
       n_taxa = length(phyDat), n_chars = ncol(tip_data))
}

# ---- Section 1: Load matrices and summarise ---------------------------------

cat("=== T-069 Large-Matrix Stress Test ===\n\n")
cat("=== Section 1: Dataset summary ===\n\n")

datasets <- list()
for (f in MATRICES) {
  cat("  Loading", f, "...\n")
  pd <- load_nex(f)
  ds <- prep_ds(pd)
  datasets[[f]] <- ds
  inappl_pct <- if (!is.null(attributes(pd)$levels) &&
                     "-" %in% attributes(pd)$levels) {
    round(100 * mean(unlist(pd) == which(attributes(pd)$levels == "-")), 1)
  } else 0
  cat(sprintf("    %s: %d taxa, %d chars, inapplicable_pct=%.1f%%\n",
              f, ds$n_taxa, ds$n_chars, inappl_pct))
}

# ---- Section 2: TBR pass micro-benchmark ------------------------------------

cat("\n=== Section 2: TBR pass micro-benchmark (ts_bench_tbr_phases) ===\n\n")

tbr_results <- list()
for (f in MATRICES) {
  ds <- datasets[[f]]
  cat(sprintf("  %s (%d tips)...\n", f, ds$n_taxa))

  reps_raw <- vector("list", 3)
  for (i in 1:3) {
    set.seed(4100 + i)
    tree <- RandomTree(ds$n_taxa, root = TRUE)
    reps_raw[[i]] <- TreeSearch:::ts_bench_tbr_phases(
      tree$edge,
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      ds$min_steps
    )
  }

  avg <- function(field) mean(vapply(reps_raw, `[[`, numeric(1), field))
  row <- data.frame(
    file         = f,
    n_tips       = reps_raw[[1]]$n_tips,
    n_blocks     = reps_raw[[1]]$n_blocks,
    total_words  = reps_raw[[1]]$total_words,
    has_na       = reps_raw[[1]]$has_na,
    n_clips      = avg("n_clips"),
    n_candidates = avg("n_candidates"),
    full_rescore_us   = avg("time_full_rescore_us"),
    clip_incr_us      = avg("time_clip_incr_us"),
    indirect_us       = avg("time_indirect_us"),
    unclip_us         = avg("time_unclip_us"),
    snap_save_us      = avg("time_snapshot_save_us"),
    snap_restore_us   = avg("time_snapshot_restore_us"),
    snap_bytes        = avg("snapshot_bytes"),
    stringsAsFactors  = FALSE
  )
  tbr_results[[f]] <- row
  cat(sprintf("    clips=%.0f  cands=%.0f  indirect=%.0fms  snap=%.1fKB\n",
              row$n_clips, row$n_candidates,
              row$indirect_us / 1000, row$snap_bytes / 1024))
}

tbr_df <- do.call(rbind, tbr_results)
rownames(tbr_df) <- NULL

cat("\nTBR phase timing (μs, per pass):\n")
print(tbr_df[, c("file", "n_tips", "n_blocks", "full_rescore_us",
                  "clip_incr_us", "indirect_us", "unclip_us",
                  "snap_save_us", "snap_restore_us")], digits = 4)

cat("\nPer-candidate indirect timing (ns):\n")
ns_cand <- round(1000 * tbr_df$indirect_us / tbr_df$n_candidates, 1)
print(data.frame(file = tbr_df$file, n_tips = tbr_df$n_tips,
                 n_candidates = round(tbr_df$n_candidates),
                 indirect_total_ms = round(tbr_df$indirect_us / 1000, 1),
                 ns_per_candidate = ns_cand))

# ---- Section 3: Scaling vs smaller datasets --------------------------------
#
# Pull synthetic-series data from bench_memory.R baselines if available,
# otherwise run a quick synthetic series here.

cat("\n=== Section 3: Scaling analysis ===\n\n")

# Quick synthetic series: 20, 50, 100, 200, + new 225 point from tbr_df
make_synthetic <- function(n_tips, n_chars = 200, na_prob = 0.1) {
  tree <- RandomTree(n_tips, root = TRUE)
  mat  <- matrix(
    sample(c("0", "1", "-"), n_tips * n_chars, replace = TRUE,
           prob = c((1 - na_prob) / 2, (1 - na_prob) / 2, na_prob)),
    n_tips, n_chars,
    dimnames = list(tree$tip.label, NULL)
  )
  MatrixToPhyDat(mat)
}

bench_tbr_one <- function(n_tips, n_chars = 200, na_prob = 0.1, seed = 4200) {
  set.seed(seed)
  pd <- make_synthetic(n_tips, n_chars, na_prob)
  ds <- prep_ds(pd)
  set.seed(seed + 1)
  tree <- RandomTree(n_tips, root = TRUE)
  r <- TreeSearch:::ts_bench_tbr_phases(
    tree$edge,
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    ds$min_steps
  )
  data.frame(
    n_tips       = n_tips,
    n_candidates = r$n_candidates,
    indirect_us  = r$time_indirect_us,
    clip_incr_us = r$time_clip_incr_us
  )
}

synth_sizes <- c(20, 50, 100, 150, 200, 225)
cat("  Synthetic scaling series:", paste(synth_sizes, collapse = ", "), "tips...\n")
synth_rows <- lapply(synth_sizes, function(n) {
  cat("    n =", n, "\n")
  bench_tbr_one(n)
})
synth_df <- do.call(rbind, synth_rows)
print(synth_df)

# Fit scaling exponents
if (nrow(synth_df) >= 4) {
  fit_indirect   <- lm(log(indirect_us)  ~ log(n_tips), data = synth_df)
  fit_candidates <- lm(log(n_candidates) ~ log(n_tips), data = synth_df)
  fit_clip       <- lm(log(clip_incr_us) ~ log(n_tips), data = synth_df)
  cat(sprintf("\nScaling exponents (log-log fit):\n"))
  cat(sprintf("  indirect_us   ~ n^%.2f  (expected ~2.0)\n", coef(fit_indirect)[2]))
  cat(sprintf("  n_candidates  ~ n^%.2f  (expected ~2.0)\n", coef(fit_candidates)[2]))
  cat(sprintf("  clip_incr_us  ~ n^%.2f\n", coef(fit_clip)[2]))
}

# ---- Section 4: Full driven search (default params, 2 seeds) ---------------

cat("\n=== Section 4: Full driven search at default params ===\n\n")
cat("  (maxReplicates=2, nThreads=2, default strategy)\n\n")

driven_results <- list()
for (f in MATRICES) {
  ds <- datasets[[f]]
  cat(sprintf("--- %s (%d tips, %d chars) ---\n", f, ds$n_taxa, ds$n_chars))

  # Auto-select strategy: replicate what MaximizeParsimony() does
  # For large matrices, thorough if nChar < 100 AND nTip >= 65
  nTip  <- ds$n_taxa
  nChar <- ds$n_chars
  use_thorough <- (nTip >= 65) && (nChar < 100)
  if (use_thorough) {
    ratchet <- 20L; drift <- 12L; xss <- 1L; rss <- 1L; css <- 0L
    strat_name <- "thorough"
  } else {
    ratchet <- 5L; drift <- 2L; xss <- 1L; rss <- 1L; css <- 0L
    strat_name <- "default"
  }
  cat(sprintf("  Strategy: %s (ratchet=%d, drift=%d)\n", strat_name, ratchet, drift))

  run_list <- list()
  for (seed_i in 1:2) {
    set.seed(4300 + seed_i)
    t0 <- proc.time()
    result <- TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates  = 2L,
      targetHits     = 1L,
      ratchetCycles  = ratchet,
      driftCycles    = drift,
      xssRounds      = xss,
      rssRounds      = rss,
      cssRounds      = css,
      cssPartitions  = 3L,
      xssPartitions  = 3L,
      fuseInterval   = 5L,
      maxSeconds     = 300,
      verbosity      = 0L,
      nThreads       = 2L
    )
    elapsed <- (proc.time() - t0)[3]
    run_list[[seed_i]] <- list(result = result, elapsed = elapsed)
    cat(sprintf("  seed %d: %.2fs  score=%.1f  reps=%d  pool=%d\n",
                seed_i, elapsed, result$best_score,
                result$replicates, result$pool_size))
  }

  # Per-phase breakdown from first run
  r1 <- run_list[[1]]$result
  if (!is.null(r1$timings)) {
    timings <- r1$timings
    cpp_total <- sum(unlist(timings))
    cat(sprintf("  Per-phase breakdown (seed 1):\n"))
    for (ph in names(timings)) {
      pct <- if (cpp_total > 0) 100 * timings[[ph]] / cpp_total else 0
      cat(sprintf("    %-12s %7.0f ms  (%4.1f%%)\n", ph, timings[[ph]], pct))
    }
    cat(sprintf("    %-12s %7.0f ms  (C++ total)\n", "TOTAL", cpp_total))
  }

  driven_results[[f]] <- list(
    file     = f,
    n_tips   = ds$n_taxa,
    n_chars  = ds$n_chars,
    strategy = strat_name,
    score1   = run_list[[1]]$result$best_score,
    score2   = run_list[[2]]$result$best_score,
    time1    = run_list[[1]]$elapsed,
    time2    = run_list[[2]]$elapsed,
    pool1    = run_list[[1]]$result$pool_size,
    reps1    = run_list[[1]]$result$replicates
  )
  cat("\n")
}

cat("=== Summary table ===\n\n")
summary_df <- do.call(rbind, lapply(driven_results, as.data.frame))
rownames(summary_df) <- NULL
print(summary_df)

# Save results
out_path <- "inst/benchmarks/stress_large_results.csv"
write.csv(summary_df, out_path, row.names = FALSE)
cat("\nResults written to", out_path, "\n")

cat("\n=== T-069 complete ===\n")
