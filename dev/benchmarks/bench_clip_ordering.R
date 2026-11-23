# bench_clip_ordering.R
#
# Benchmark comparison of TBR clip ordering strategies.
#
# Compares six clip_order variants:
#   0 = RANDOM     (current default)
#   1 = INV_WEIGHT (w = 1/(1+s))
#   2 = TIPS_FIRST (tips first, then rest; shuffled within)
#   3 = BUCKET     (tips / small / large buckets; shuffled within)
#   4 = ANTI_TIP   (non-tips first, then tips)
#   5 = LARGE_FIRST (large > √n first, then small, then tips)
#
# Metric: time-adjusted expected best (TAEB) score — the expected minimum
# score from k independent replicates where k = floor(budget / rep_time).
# Estimated via bootstrap resampling of per-replicate scores.
#
# Usage: Rscript dev/benchmarks/bench_clip_ordering.R [lib_path] [n_seeds]
#   lib_path defaults to ".agent-wc"
#   n_seeds  defaults to 10

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else ".agent-wc"
n_seeds  <- if (length(args) >= 2) as.integer(args[2]) else 10L

library(TreeSearch, lib.loc = lib_path)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATASETS <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")
BUDGETS  <- c(30, 60)   # seconds
SEEDS    <- seq_len(n_seeds) * 1000L + 847L

ORDERS   <- c(
  RANDOM     = 0L,
  INV_WEIGHT = 1L,
  TIPS_FIRST = 2L,
  BUCKET     = 3L,
  ANTI_TIP   = 4L,
  LARGE_FIRST= 5L
)

# Expected-best bootstrap: given per-replicate scores and total wall time,
# estimate E[min score] at each time budget.
taeb <- function(scores, times_ms, budgets_s, n_boot = 2000L) {
  stopifnot(length(scores) == length(times_ms))
  n <- length(scores)
  if (n == 0) return(setNames(rep(NA_real_, length(budgets_s)), budgets_s))

  # Mean time per replicate (use median to be robust to outliers)
  med_time_s <- median(times_ms) / 1000
  if (med_time_s <= 0) med_time_s <- 1

  vapply(budgets_s, function(budget) {
    k <- max(1L, floor(budget / med_time_s))
    if (k >= n) {
      # Can't bootstrap more than n replicates; just return min
      return(min(scores))
    }
    boot_mins <- replicate(n_boot, min(sample(scores, k, replace = TRUE)))
    mean(boot_mins)
  }, numeric(1L))
}

# ---------------------------------------------------------------------------
# Prepare datasets
# ---------------------------------------------------------------------------

prepare <- function(name) {
  ds <- TreeSearch::inapplicable.phyData[[name]]
  at <- attributes(ds)
  list(
    name     = name,
    contrast = at$contrast,
    tip_data = matrix(unlist(ds, use.names = FALSE),
                      nrow = length(ds), byrow = TRUE),
    weight   = at$weight,
    levels   = at$levels,
    n_taxa   = length(ds)
  )
}

# ---------------------------------------------------------------------------
# Build a default SearchControl with preset based on n_tip
# ---------------------------------------------------------------------------

make_sc <- function(n_tip, clip_order_int = 0L) {
  # Mirror the "default" preset for datasets in 31-119 tip range,
  # with only the clip_order changed. This gives a realistic context
  # (same ratchet/XSS/RSS settings as normal use) for the comparison.
  #
  # NOTE: maxSeconds is set per run; runtimeConfig controls the budget.
  # Here we only set SearchControl parameters.
  SearchControl(
    tbrMaxHits      = 1L,
    clipOrder       = clip_order_int,
    nniFirst        = TRUE,
    sprFirst        = FALSE,
    wagnerStarts    = 1L,
    wagnerBias      = 0L,
    outerCycles     = 1L,
    maxOuterResets  = 0L,
    ratchetCycles   = 12L,
    ratchetPerturbProb   = 0.25,
    ratchetPerturbMaxMoves = 5L,
    ratchetAdaptive = FALSE,
    nniPerturbCycles = 0L,
    driftCycles     = 0L,
    xssRounds       = 3L, xssPartitions = 4L,
    rssRounds       = 1L,
    cssRounds       = 0L,
    fuseInterval    = 3L,
    adaptiveLevel   = TRUE,
    consensusStableReps = 0L
  )
}

make_runtime <- function(max_seconds) {
  list(maxReplicates = 9999L, targetHits = 9999L,
       maxSeconds = max_seconds, verbosity = 0L, nThreads = 1L)
}

# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

cat(sprintf("Benchmark: %d datasets × %d ordering variants × %d seeds\n",
    length(DATASETS), length(ORDERS), n_seeds))
cat(sprintf("Budgets: %s seconds\n\n", paste(BUDGETS, collapse = ", ")))

all_results <- list()

for (dname in DATASETS) {
  d <- prepare(dname)

  cat(sprintf("=== %s (n_tip=%d) ===\n", dname, d$n_taxa))

  ds_results <- list()

  for (oname in names(ORDERS)) {
    oint <- ORDERS[[oname]]
    sc   <- make_sc(d$n_taxa, oint)

    rep_scores <- numeric(n_seeds)
    rep_times  <- numeric(n_seeds)   # ms per replicate

    for (i in seq_along(SEEDS)) {
      set.seed(SEEDS[i])
      res <- TreeSearch:::ts_driven_search(
        d$contrast, d$tip_data, d$weight, d$levels,
        searchControl = sc,
        runtimeConfig = make_runtime(max(BUDGETS)),
        scoringConfig = list(min_steps = integer(), concavity = -1.0,
                             xpiwe = FALSE, xpiwe_r = 0.5, xpiwe_max_f = 5.0,
                             obs_count = integer(), infoAmounts = NULL)
      )
      rep_scores[i] <- res$best_score
      # Estimate per-replicate time from timings
      t_total_ms <- sum(unlist(res$timings))
      n_reps     <- max(1L, res$n_replicates)
      rep_times[i] <- t_total_ms / n_reps
    }

    taeb_vals <- taeb(rep_scores, rep_times, BUDGETS)
    cat(sprintf("  %-12s: scores %s  med_rep %.1fs  TAEB@%ds=%.1f  @%ds=%.1f\n",
        oname,
        paste(sprintf("%.0f", range(rep_scores)), collapse = "-"),
        median(rep_times)/1000,
        BUDGETS[1], taeb_vals[1],
        BUDGETS[2], taeb_vals[2]))

    ds_results[[oname]] <- list(
      order = oname, order_int = oint,
      dataset = dname, n_tip = d$n_taxa,
      scores = rep_scores, times_ms = rep_times,
      taeb = taeb_vals
    )
  }

  all_results[[dname]] <- ds_results
  cat("\n")
}

# ---------------------------------------------------------------------------
# Summary: Δ vs RANDOM baseline for each variant, averaged across datasets
# ---------------------------------------------------------------------------

cat("=== Summary: TAEB Δ vs RANDOM (lower = better) ===\n\n")

for (budget in BUDGETS) {
  cat(sprintf("--- Budget: %ds ---\n", budget))
  cat(sprintf("  %-15s", ""))
  for (dname in DATASETS) cat(sprintf(" %13s", dname))
  cat(sprintf(" %13s\n", "mean_delta"))

  for (oname in names(ORDERS)) {
    if (oname == "RANDOM") next
    cat(sprintf("  %-15s", oname))
    deltas <- numeric(length(DATASETS))
    for (j in seq_along(DATASETS)) {
      dname    <- DATASETS[j]
      ref      <- all_results[[dname]][["RANDOM"]]$taeb[[which(BUDGETS == budget)]]
      this_val <- all_results[[dname]][[oname]]$taeb[[which(BUDGETS == budget)]]
      delta    <- this_val - ref    # positive = worse (more steps)
      deltas[j] <- delta
      cat(sprintf(" %+13.2f", delta))
    }
    cat(sprintf(" %+13.2f\n", mean(deltas)))
  }
  cat("\n")
}

cat("Positive Δ = worse than RANDOM; negative Δ = better than RANDOM.\n")
cat("Done.\n")
