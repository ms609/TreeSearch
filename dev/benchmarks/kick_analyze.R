#!/usr/bin/env Rscript
# Analyse the kick_anytime corpus A/B (fixed5 / auto / scaled).
#
# Reads the per-cell CSVs (improvement-event trace + summary fields) written by
# kick_anytime.R and answers the MISSION question: does ratchetPerturbMaxMoves
# 5 -> auto/scaled help or REGRESS wall-clock time-to-optimum, per size tier?
#
# Per (matrix, seed) cell the TARGET is the union-best final score across arms
# (self-contained -- no external best-known needed). For each arm:
#   * reached      = final_score <= target
#   * tt_hit       = first elapsed_s whose best_score <= target (NA if censored)
#   * wall_total_s = total wall to natural stop
# Anytime win  = smaller tt_hit at equal-or-better reach.
# Regression   = larger tt_hit / worse reach (watch the SMALL tier: auto=20 kick).
#
# Usage: Rscript kick_analyze.R <dir-of-cell-CSVs>

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else "dev/benchmarks/kick_anytime_out"
files <- list.files(dir, pattern = "^cell_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("no cell CSVs in ", dir)
D <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
cat(sprintf("Loaded %d rows from %d cells\n", nrow(D), length(files)))

cells <- unique(D[, c("dataset", "nTip", "nChar", "tier", "seed")])
arms <- c("fixed5", "auto", "scaled")

per_arm <- list()
for (i in seq_len(nrow(cells))) {
  c0 <- cells[i, ]
  sub <- D[D$dataset == c0$dataset & D$seed == c0$seed, ]
  finals <- sapply(arms, function(a) {
    s <- sub[sub$arm == a, ]
    if (!nrow(s)) NA_real_ else s$final_score[1]
  })
  target <- suppressWarnings(min(finals, na.rm = TRUE))
  for (a in arms) {
    s <- sub[sub$arm == a, ]
    if (!nrow(s)) next
    fin <- s$final_score[1]
    hit_events <- s[is.finite(s$best_score) & s$best_score <= target + 1e-9, ]
    tt <- if (nrow(hit_events)) min(hit_events$elapsed_s, na.rm = TRUE) else NA_real_
    per_arm[[length(per_arm) + 1L]] <- data.frame(
      dataset = c0$dataset, nTip = c0$nTip, tier = c0$tier, seed = c0$seed,
      arm = a, target = target, final = fin, reached = fin <= target + 1e-9,
      tt_hit = tt, wall_total = s$wall_total_s[1], reps = s$reps_done[1],
      stringsAsFactors = FALSE)
  }
}
PA <- do.call(rbind, per_arm)

tier_ord <- c("small", "medium", "large", "xlarge")
PA$tier <- factor(PA$tier, levels = tier_ord)

cat("\n=== Reach: cells where arm's final == union-best (per tier) ===\n")
reach <- tapply(PA$reached, list(PA$tier, PA$arm), mean)
print(round(reach, 3))

cat("\n=== Median time-to-first-hit of union-best (s), reached cells only ===\n")
r <- PA[PA$reached, ]
tt <- tapply(r$tt_hit, list(r$tier, r$arm), function(x) median(x, na.rm = TRUE))
print(round(tt, 3))

cat("\n=== Median total wall to natural stop (s) ===\n")
wt <- tapply(PA$wall_total, list(PA$tier, PA$arm), function(x) median(x, na.rm = TRUE))
print(round(wt, 2))

cat("\n=== Anytime ratio vs fixed5 (median tt_hit arm / tt_hit fixed5), paired per cell ===\n")
wide <- reshape(PA[, c("dataset", "seed", "tier", "arm", "tt_hit", "final", "reached")],
                idvar = c("dataset", "seed", "tier"), timevar = "arm", direction = "wide")
for (tt_col in c("tt_hit.auto", "tt_hit.scaled")) {
  ok <- is.finite(wide[[tt_col]]) & is.finite(wide$tt_hit.fixed5) & wide$tt_hit.fixed5 > 0
  ratio <- wide[[tt_col]][ok] / wide$tt_hit.fixed5[ok]
  cat(sprintf("  %-16s median ratio=%.3f  (n=%d; <1 => faster than fixed5)\n",
              sub("tt_hit\\.", "", tt_col), median(ratio, na.rm = TRUE), sum(ok)))
  by_t <- tapply(ratio, wide$tier[ok], function(x) round(median(x, na.rm = TRUE), 3))
  print(by_t)
}

cat("\n=== Reach wins/losses vs fixed5 (strictly better/worse FINAL score) ===\n")
for (a in c("auto", "scaled")) {
  fc <- paste0("final.", a)
  ok <- is.finite(wide[[fc]]) & is.finite(wide$final.fixed5)
  wins <- sum(wide[[fc]][ok] < wide$final.fixed5[ok] - 1e-9)
  loss <- sum(wide[[fc]][ok] > wide$final.fixed5[ok] + 1e-9)
  cat(sprintf("  %-7s: %d cells better, %d worse, %d tie (of %d)\n",
              a, wins, loss, sum(ok) - wins - loss, sum(ok)))
}

saveRDS(list(PA = PA, wide = wide), file.path(dir, "kick_analysis.rds"))
write.csv(PA, file.path(dir, "kick_per_arm.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s\n", file.path(dir, "kick_per_arm.csv")))
