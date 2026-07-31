#!/usr/bin/env Rscript
# Analyse the reach-escalation A/B. Usage: Rscript reach_escalation_analyze.R <dir>
#
# Per (matrix, seed) cell the TARGET is the union-best final score across arms (the
# established mbank convention -- no canonical best-known table exists for mbank).
# Per arm:
#   reached  = final_score <= target
#   tt_hit   = first elapsed_s   whose best_score <= target   (NA if never reached)
#   rep2hit  = first replicate   whose best_score <= target   (NA if never reached)
# Anytime win = smaller tt_hit at equal-or-better reach; regression = larger tt_hit /
# worse reach. rep2hit is reported ALONGSIDE tt_hit because a wall-only read
# misattributes "costlier per replicate" as "slower to the optimum" (the 2026-07-14
# kick study: wall2hit -7..-12% but rep2hit == 1.000).
#
# PRE-REGISTERED RULE (fixed before results):
#   SHIP    : reach(deltas) >= reach(base) overall AND no tier reach-regression.
#   NO-SHIP : reach(deltas) <  reach(base) on the general pool.
# Wall cost is accepted by construction (the gate is opt-in; the user raised
# targetHits meaning "time is no issue") but is reported honestly.

args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1L) args[[1]] else "."
files <- list.files(dir, pattern = "^cell_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop("no cell_*.csv in ", dir)
D <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
cat(sprintf("Loaded %d rows from %d cell files\n", nrow(D), length(files)))

# The engine's real stopping deadline, not the nominal cap (see the at_deadline comment).
ENUM_TIME_FRACTION_DEFAULT <- 0.1
deadline_of <- function(sa) {
  etf <- if ("enum_time_fraction" %in% names(sa) && !is.na(sa$enum_time_fraction[1])) {
    sa$enum_time_fraction[1]
  } else ENUM_TIME_FRACTION_DEFAULT
  sa$cap_s[1] * (1 - etf)
}
if (!"enum_time_fraction" %in% names(D))
  cat(sprintf("NOTE: no enum_time_fraction column; assuming the %.2f default for deadlines\n",
              ENUM_TIME_FRACTION_DEFAULT))

cells <- unique(D[, c("dataset", "nTip", "tier", "seed")])
out <- list()
for (i in seq_len(nrow(cells))) {
  cl <- cells[i, ]
  sub <- D[D$dataset == cl$dataset & D$seed == cl$seed, ]
  target <- min(sub$final_score, na.rm = TRUE)          # union-best across arms
  for (a in unique(sub$arm)) {
    sa <- sub[sub$arm == a, ]
    sa <- sa[order(sa$elapsed_s), ]
    hit <- which(sa$best_score <= target + 1e-9)
    out[[length(out) + 1L]] <- data.frame(
      dataset = cl$dataset, nTip = cl$nTip, tier = cl$tier, seed = cl$seed, arm = a,
      target = target, final = sa$final_score[1],
      reached = as.integer(sa$final_score[1] <= target + 1e-9),
      tt_hit  = if (length(hit)) sa$elapsed_s[hit[1]] else NA_real_,
      rep2hit = if (length(hit)) sa$replicate[hit[1]] else NA_integer_,
      wall_total = sa$wall_total_s[1], reps = sa$reps_done[1],
      cap_s = sa$cap_s[1],
      # DEADLINE FLAG. The deltas arm costs ~2.3x the wall per replicate, so at a fixed
      # budget it completes far fewer reps. An arm that stopped at the budget did not
      # converge -- reading its "reach" as a property of the config is exactly the error
      # that made 5432 arm B look replicate-capped when it was time-truncated.
      #
      # The budget is NOT maxSeconds. The engine stops the main search at
      #   main_deadline = maxSeconds * (1 - enumTimeFraction)         [src/ts_driven.cpp]
      # with enumTimeFraction defaulting to 0.1, so a deadline-bound cell lands at ~0.90 *
      # cap_s and NOT at cap_s. This flag originally tested `>= 0.95 * cap_s` and therefore
      # scored 59 of 110 genuinely deadline-bound ab6 cells as "converged" -- it missed the
      # very trap documented in reach_escalation_FINDINGS.md. Read enumTimeFraction from the
      # data when the harness records it, else assume the 0.1 default.
      deadline_s = deadline_of(sa),
      at_deadline = as.integer(!is.na(sa$wall_total_s[1]) &&
                               sa$wall_total_s[1] >= 0.98 * deadline_of(sa)),
      stringsAsFactors = FALSE)
  }
}
P <- do.call(rbind, out)
write.csv(P, file.path(dir, "reach_ab_per_arm.csv"), row.names = FALSE)

fmt <- function(x) if (all(is.na(x))) "NA" else sprintf("%.3g", median(x, na.rm = TRUE))
cat("\n=== REACH by arm (fraction of cells attaining the cell's union-best) ===\n")
for (a in sort(unique(P$arm)))
  cat(sprintf("  %-6s reach = %.3f (%d/%d)\n", a,
              mean(P$arm == a & P$reached == 1L) / mean(P$arm == a),
              sum(P$arm == a & P$reached == 1L), sum(P$arm == a)))

cat("\n=== REACH by tier (the tier-regression check) ===\n")
for (tr in c("small", "medium", "large", "xlarge")) {
  s <- P[P$tier == tr, ]
  if (!nrow(s)) next
  cat(sprintf("  %-7s", tr))
  for (a in sort(unique(P$arm)))
    cat(sprintf("  %s=%.3f (%d/%d)", a, mean(s$reached[s$arm == a]),
                sum(s$reached[s$arm == a]), sum(s$arm == a)))
  cat("\n")
}

cat("\n=== MEDIAN tt_hit (wall s) / rep2hit (replicates) / total wall ===\n")
for (a in sort(unique(P$arm)))
  cat(sprintf("  %-6s tt_hit=%-8s rep2hit=%-6s wall_total=%-8s reps=%s\n", a,
              fmt(P$tt_hit[P$arm == a]), fmt(P$rep2hit[P$arm == a]),
              fmt(P$wall_total[P$arm == a]), fmt(P$reps[P$arm == a])))

# Paired per-cell comparison (repo convention: counts + median ratios, not p-values).
# Arm names are READ FROM THE DATA rather than hard-coded, so a variant harness (e.g. a
# `deltas6` arm confirming the shipped six-lever form) analyses without an edit -- with
# fixed names the paired table silently comes back empty instead of erroring.
armNames <- sort(unique(P$arm))
baseArm <- if ("base" %in% armNames) "base" else armNames[[1]]
testArm <- setdiff(armNames, baseArm)[[1]]
b <- P[P$arm == baseArm, ];   d <- P[P$arm == testArm, ]
k <- intersect(paste(b$dataset, b$seed), paste(d$dataset, d$seed))
b <- b[match(k, paste(b$dataset, b$seed)), ]; d <- d[match(k, paste(d$dataset, d$seed)), ]
cat(sprintf("\n=== PAIRED (n = %d cells; baseline '%s' vs test '%s') ===\n",
            length(k), baseArm, testArm))
cat(sprintf("  final score : %d better, %d worse, %d tie  (%s vs %s)\n",
            sum(d$final < b$final), sum(d$final > b$final), sum(d$final == b$final),
            testArm, baseArm))
cat(sprintf("  reach       : %s %d, %s %d\n", baseArm, sum(b$reached),
            testArm, sum(d$reached)))
ok <- !is.na(b$tt_hit) & !is.na(d$tt_hit)
if (any(ok)) {
  cat(sprintf("  tt_hit  ratio (deltas/base), median = %.3f  [%d better, %d worse]\n",
              median(d$tt_hit[ok] / b$tt_hit[ok]),
              sum(d$tt_hit[ok] < b$tt_hit[ok]), sum(d$tt_hit[ok] > b$tt_hit[ok])))
  okr <- ok & !is.na(b$rep2hit) & !is.na(d$rep2hit) & b$rep2hit > 0
  if (any(okr))
    cat(sprintf("  rep2hit ratio (deltas/base), median = %.3f  <- if ~1.0, any wall gap is COST-PER-REP, not slower reach\n",
                median(d$rep2hit[okr] / b$rep2hit[okr])))
}
cat(sprintf("  wall_total ratio (deltas/base), median = %.3f\n",
            median(d$wall_total / b$wall_total, na.rm = TRUE)))

worse <- d$final > b$final
if (any(worse)) {
  cat("\n!! cells where deltas found a WORSE tree (the NO-SHIP signal):\n")
  print(data.frame(dataset = d$dataset[worse], tier = d$tier[worse], seed = d$seed[worse],
                   base = b$final[worse], deltas = d$final[worse]), row.names = FALSE)
} else cat("\nNo cell where deltas found a worse tree.\n")

# BUDGET REGIME. Being deadline-bound is not automatically a spoiled cell: when BOTH arms
# stop at the same deadline the cell is a valid EQUAL-WALL comparison, which is the stronger
# test (the cheaper-per-rep arm gets more replicates and still has to win). What invalidates
# a cell is ASYMMETRY -- one arm converged and the other was cut off -- because then the
# score gap may be purely budget. So classify cells three ways instead of dropping them.
cat(sprintf("\n=== BUDGET REGIME (deadline = cap_s x (1 - enumTimeFraction)) ===\n"))
cat(sprintf("  base   %d/%d cells stopped at the deadline\n", sum(b$at_deadline), nrow(b)))
cat(sprintf("  %-6s %d/%d cells stopped at the deadline\n", testArm,
            sum(d$at_deadline), nrow(d)))
bothDL <- b$at_deadline == 1L & d$at_deadline == 1L
neither <- b$at_deadline == 0L & d$at_deadline == 0L
asym <- !bothDL & !neither
cat(sprintf("  both at deadline (valid, EQUAL-WALL) : %d\n", sum(bothDL)))
cat(sprintf("  neither (both converged, valid)      : %d\n", sum(neither)))
cat(sprintf("  exactly one (ASYMMETRIC, suspect)    : %d\n", sum(asym)))
if (any(asym))
  print(data.frame(dataset = d$dataset[asym], seed = d$seed[asym],
                   baseAtDL = b$at_deadline[asym], testAtDL = d$at_deadline[asym],
                   base = b$final[asym], test = d$final[asym]), row.names = FALSE)
if (any(bothDL)) {
  rr <- b$reps[bothDL] / d$reps[bothDL]
  rr <- rr[is.finite(rr)]
  if (length(rr))
    cat(sprintf("  work per replicate on equal-wall cells: %s does %.2fx the replicates\n",
                baseArm, median(rr)))
}

# The verdict is computed on the VALID cells (both-at-deadline plus both-converged) and is
# driven by PAIRED SCORE COUNTS, not by reach. Reach here is measured against the union-best
# across arms, which is self-referential -- if one arm alone attains a score the other
# "misses" by construction -- so a reach fraction restates the paired counts with the losses
# inflated. Both are printed; the counts are the statistic.
valid <- bothDL | neither
nBetter <- sum(d$final[valid] < b$final[valid])
nWorse  <- sum(d$final[valid] > b$final[valid])
reachB <- mean(b$reached[valid]); reachD <- mean(d$reached[valid])
tierBad <- character(0)
for (tr in unique(P$tier)) {
  iv <- valid & b$tier == tr
  if (sum(iv) && sum(d$final[iv] > b$final[iv]) > sum(d$final[iv] < b$final[iv]))
    tierBad <- c(tierBad, tr)
}
cat(sprintf("\n=== PRE-REGISTERED VERDICT (valid cells, n = %d of %d) ===\n",
            sum(valid), length(valid)))
cat(sprintf("  paired score: %d better, %d worse; tier regressions: %s\n",
            nBetter, nWorse, if (length(tierBad)) paste(tierBad, collapse = ",") else "none"))
cat(sprintf("  reach (union-best, self-referential): base=%.3f %s=%.3f\n",
            reachB, testArm, reachD))
cat(sprintf("  --> %s\n",
            if (!sum(valid)) "INCONCLUSIVE -- no valid cells"
            else if (nBetter >= nWorse && !length(tierBad)) "SHIP"
            else "NO-SHIP (restrict or drop)"))

# How concentrated is the effect? A tier-level win can be a single matrix repeated across
# seeds -- that happened here (project4284) and was briefly written up as a tier property.
chg <- valid & d$final != b$final
if (any(chg)) {
  cat("\n=== WHERE THE EFFECT LIVES (per matrix; a 1-matrix effect is NOT a tier property) ===\n")
  for (ds in unique(P$dataset[P$dataset %in% b$dataset[chg]])) {
    i <- valid & b$dataset == ds
    cat(sprintf("  %-18s %5dt  win %d  loss %d  tie %d\n", ds, b$nTip[i][1],
                sum(d$final[i] < b$final[i]), sum(d$final[i] > b$final[i]),
                sum(d$final[i] == b$final[i])))
  }
  cat(sprintf("  matrices with any change: %d of %d in the battery\n",
              length(unique(b$dataset[chg])), length(unique(b$dataset))))
}
