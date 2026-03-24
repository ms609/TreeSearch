# T-199: PT diagnostic profiling across size range
# Run ts_parallel_temper_diag() on 8 representative matrices with default
# temperature ladder {0, 3, 9, 27} to identify where swaps become useful.

.libPaths(c("../.builds/TreeSearch-C", .libPaths()))
library(TreeSearch.C)
library(TreeTools)

ts_parallel_temper_diag <- get("ts_parallel_temper_diag",
                               envir = asNamespace("TreeSearch.C"))

matrices <- c(
  small_1  = "project3419.nex",    # 40 tips, 368 chars
  small_2  = "project2604.nex",    # 43 tips, 307 chars
  medium_1 = "project4531.nex",    # 71 tips, 256 chars
  medium_2 = "project3741.nex",    # 86 tips, 110 chars
  large_1  = "project3253.nex",    # 125 tips, 394 chars
  large_2  = "project1221.nex",    # 150 tips, 252 chars
  xlarge_1 = "project3763.nex",    # 205 tips, 105 chars
  xlarge_2 = "project2722.nex"     # 385 tips, 520 chars
)

prepare_for_diag <- function(dataset) {
  at <- attributes(dataset)
  contrast <- at$contrast
  tip_data <- matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE)
  list(
    contrast  = contrast,
    tip_data  = tip_data,
    weight    = at$weight,
    levels    = at$levels,
    min_steps = integer(0)
  )
}

n_chains <- 4L
temps <- c(0, 3, 9, 27)
rounds <- 10L
moves_per_round <- 0L  # default = n_tip

results <- list()

for (nm in names(matrices)) {
  cat("\n=== ", nm, " (", matrices[nm], ") ===\n", sep = "")
  path <- file.path("../matrixes", matrices[nm])
  if (!file.exists(path)) {
    cat("  SKIP: file not found\n")
    next
  }

  res <- tryCatch({
  dataset <- ReadAsPhyDat(path)
  ntip <- length(dataset)
  nchar_total <- sum(attr(dataset, "weight"))
  cat("  Tips:", ntip, " Chars:", nchar_total, "\n")

  pd <- prepare_for_diag(dataset)

  set.seed(8127)
  t0 <- proc.time()
  diag <- ts_parallel_temper_diag(
    contrast   = pd$contrast,
    tip_data   = pd$tip_data,
    weight     = pd$weight,
    levels     = pd$levels,
    min_steps  = pd$min_steps,
    concavity  = -1,
    n_chains   = n_chains,
    temperatures = temps,
    rounds     = rounds,
    moves_per_round = moves_per_round
  )
  elapsed <- (proc.time() - t0)[["elapsed"]]

  cl <- diag$chain_log
  sl <- diag$swap_log

  cat("  Best score:", diag$best_score, "\n")
  cat("  Cold final:", diag$cold_final, "\n")
  cat("  Cold improvements from swaps:", diag$cold_improvements, "\n")
  cat("  Swap acceptance (overall):",
      diag$swaps_accepted, "/", diag$swaps_attempted, "\n")

  if (nrow(sl) > 0) {
    pair_rates <- tapply(sl$accepted, paste(sl$pair_lo, sl$pair_hi, sep = "-"),
                         mean)
    cat("  Pair acceptance rates:\n")
    for (p in names(pair_rates)) {
      cat("    Pair", p, ":", sprintf("%.1f%%", pair_rates[p] * 100), "\n")
    }
  }

  cat("  Timing: cold_tbr=", round(diag$cold_tbr_ms, 1), "ms",
      " hot=", round(diag$hot_stochastic_ms, 1), "ms",
      " total=", round(diag$total_pt_ms, 1), "ms",
      " (wall:", round(elapsed, 1), "s)\n")
  cat("  Hot overhead:", sprintf("%.1f%%",
      diag$hot_stochastic_ms / max(diag$total_pt_ms, 0.01) * 100), "\n")

  # Per-chain score trajectories
  cat("  Per-chain final scores:\n")
  for (ch in sort(unique(cl$chain))) {
    ch_rows <- cl[cl$chain == ch, ]
    cat("    Chain", ch, "(T=", temps[ch + 1], "):",
        "start=", ch_rows$score_before[1],
        "end=", tail(ch_rows$score_after, 1),
        "accept_rate=", sprintf("%.1f%%",
          sum(ch_rows$n_accepted) / max(sum(ch_rows$n_attempted), 1) * 100),
        "\n")
  }

  results[[nm]] <- list(
    diag = diag, ntip = ntip, nchar = nchar_total,
    chain_log = cl, swap_log = sl
  )
  NULL  # success
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
  })
}

# Summary table
cat("\n\n=== SUMMARY TABLE ===\n")
rows <- lapply(names(results), function(nm) {
  r <- results[[nm]]
  d <- r$diag
  sl <- r$swap_log

  # Pair-0-1 acceptance (cold-warm boundary)
  pair01 <- if (nrow(sl) > 0) {
    p01 <- sl[sl$pair_lo == 0 & sl$pair_hi == 1, ]
    if (nrow(p01) > 0) mean(p01$accepted) else NA
  } else NA

  data.frame(
    matrix = nm, ntip = r$ntip, nchar = r$nchar,
    best = d$best_score,
    swap_rate = if (d$swaps_attempted > 0) d$swaps_accepted / d$swaps_attempted else 0,
    pair01_rate = pair01,
    cold_improve = d$cold_improvements,
    cold_ms = d$cold_tbr_ms,
    hot_ms = d$hot_stochastic_ms,
    total_ms = d$total_pt_ms,
    stringsAsFactors = FALSE
  )
})
summary_df <- do.call(rbind, rows)
print(summary_df, digits = 3)
