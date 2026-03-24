# Phase 3D: Memory layout profiling
#
# Measures TBR phase breakdown and scaling across tree sizes.
# Run with: source("inst/benchmarks/bench_memory.R")

library(TreeSearch)
library(TreeTools)

# --- Helper: prepare dataset args for Rcpp call ---
prep_ds <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  storage.mode(contrast) <- "double"
  # phyDat stores data as list of integer vectors (one per taxon)
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  storage.mode(tip_data) <- "integer"
  weight <- at$weight
  levels <- at$levels

  # min_steps from contrast matrix
  min_steps <- apply(contrast, 2, function(x) sum(x > 0)) - 1L
  min_steps <- pmax(min_steps, 0L)

  list(contrast = contrast, tip_data = tip_data, weight = weight,
       levels = levels, min_steps = min_steps)
}

# --- Helper: get random tree edge matrix for n tips ---
make_tree_edge <- function(dataset) {
  tree <- RandomTree(names(dataset), root = TRUE)
  tree$edge
}

# --- Helper: generate synthetic dataset ---
make_synthetic <- function(n_tips, n_chars = 200, na_prob = 0.1) {
  tree <- RandomTree(n_tips, root = TRUE)
  mat <- matrix(
    sample(c("0", "1", "-"), n_tips * n_chars, replace = TRUE,
           prob = c((1 - na_prob) / 2, (1 - na_prob) / 2, na_prob)),
    n_tips, n_chars,
    dimnames = list(tree$tip.label, NULL)
  )
  MatrixToPhyDat(mat)
}

# --- Benchmark one dataset ---
bench_one <- function(dataset, label, n_reps = 3) {
  ds_args <- prep_ds(dataset)
  edge <- make_tree_edge(dataset)

  results <- vector("list", n_reps)
  for (i in seq_len(n_reps)) {
    edge <- make_tree_edge(dataset) # different random tree each rep
    results[[i]] <- TreeSearch:::ts_bench_tbr_phases(
      edge, ds_args$contrast, ds_args$tip_data,
      ds_args$weight, ds_args$levels,
      ds_args$min_steps
    )
  }

  # Average across reps
  avg <- function(field) mean(vapply(results, `[[`, numeric(1), field))

  data.frame(
    label = label,
    n_tips = results[[1]]$n_tips,
    n_node = results[[1]]$n_node,
    n_blocks = results[[1]]$n_blocks,
    total_words = results[[1]]$total_words,
    total_chars = results[[1]]$total_chars,
    has_na = results[[1]]$has_na,
    score = avg("score"),
    n_clips = avg("n_clips"),
    n_candidates = avg("n_candidates"),
    # Timing (microseconds)
    full_rescore_us = avg("time_full_rescore_us"),
    clip_incr_us = avg("time_clip_incr_us"),
    indirect_us = avg("time_indirect_us"),
    unclip_us = avg("time_unclip_us"),
    snap_save_us = avg("time_snapshot_save_us"),
    snap_restore_us = avg("time_snapshot_restore_us"),
    snap_bytes = avg("snapshot_bytes"),
    stringsAsFactors = FALSE
  )
}

# --- Run benchmarks ---
cat("=== Phase 3D Memory Layout Profiling ===\n\n")

set.seed(7382)

# Empirical datasets
cat("Benchmarking empirical datasets...\n")
data("inapplicable.phyData", package = "TreeSearch")

empirical_results <- list()
for (name in c("Vinther2008", "Agnarsson2004")) {
  cat("  ", name, "...\n")
  empirical_results[[name]] <- bench_one(
    inapplicable.phyData[[name]], name, n_reps = 3
  )
}

# Synthetic datasets of increasing size
cat("Benchmarking synthetic datasets...\n")
sizes <- c(20, 50, 100, 200)
synthetic_results <- list()
for (n in sizes) {
  label <- paste0("synth_", n)
  cat("  ", label, "...\n")
  ds <- make_synthetic(n, n_chars = 200, na_prob = 0.1)
  synthetic_results[[label]] <- bench_one(ds, label, n_reps = 3)
}

# Combine results
all_results <- do.call(rbind, c(empirical_results, synthetic_results))

# --- Display ---
cat("\n=== Results ===\n\n")
print(all_results[, c("label", "n_tips", "n_blocks", "total_words",
                       "n_clips", "n_candidates")])

cat("\n=== Timing breakdown (microseconds, total across all clips) ===\n\n")
timing_cols <- c("label", "n_tips", "full_rescore_us", "clip_incr_us",
                 "indirect_us", "unclip_us", "snap_save_us", "snap_restore_us")
print(all_results[, timing_cols], digits = 3)

# Compute fractions
cat("\n=== Time fractions (clip+incr / indirect / unclip) ===\n\n")
total_pass <- all_results$clip_incr_us + all_results$indirect_us +
              all_results$unclip_us
fracs <- data.frame(
  label = all_results$label,
  n_tips = all_results$n_tips,
  pct_clip_incr = round(100 * all_results$clip_incr_us / total_pass, 1),
  pct_indirect = round(100 * all_results$indirect_us / total_pass, 1),
  pct_unclip = round(100 * all_results$unclip_us / total_pass, 1),
  snap_save_per_op_us = round(all_results$snap_save_us, 1),
  snap_restore_per_op_us = round(all_results$snap_restore_us, 1),
  snap_KB = round(all_results$snap_bytes / 1024, 1)
)
print(fracs)

# Per-candidate timing
cat("\n=== Per-candidate indirect timing ===\n\n")
per_cand <- data.frame(
  label = all_results$label,
  n_tips = all_results$n_tips,
  n_candidates = round(all_results$n_candidates),
  indirect_us_total = round(all_results$indirect_us),
  ns_per_candidate = round(1000 * all_results$indirect_us /
                           all_results$n_candidates, 1)
)
print(per_cand)

# Scaling analysis
cat("\n=== Scaling analysis (synthetic datasets) ===\n\n")
synth <- all_results[grepl("synth", all_results$label), ]
if (nrow(synth) >= 3) {
  fit <- lm(log(indirect_us) ~ log(n_tips), data = synth)
  cat("Indirect time scaling exponent:", round(coef(fit)[2], 2),
      "(expected ~2 for O(n^2))\n")
  fit2 <- lm(log(n_candidates) ~ log(n_tips), data = synth)
  cat("Candidate count scaling exponent:", round(coef(fit2)[2], 2), "\n")
}

cat("\nDone.\n")
