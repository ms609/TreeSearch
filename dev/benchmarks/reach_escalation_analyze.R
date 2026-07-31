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
      # TRUNCATION FLAG. The deltas arm costs ~15x the candidates per replicate, so at a
      # fixed wall it completes far fewer reps. If an arm stopped AT the cap it did not
      # converge -- its "reach" is a budget artefact, not a property of the config. This
      # is exactly the error that made 5432 arm B look replicate-capped when it was
      # time-truncated. A reach comparison is only honest on non-truncated cells.
      truncated = as.integer(!is.na(sa$wall_total_s[1]) &&
                             sa$wall_total_s[1] >= 0.95 * sa$cap_s[1]),
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

cat(sprintf("\n=== TRUNCATION (stopped at the wall cap => did NOT converge) ===\n"))
cat(sprintf("  base   %d/%d cells truncated\n  deltas %d/%d cells truncated\n",
            sum(b$truncated), nrow(b), sum(d$truncated), nrow(d)))
clean <- b$truncated == 0L & d$truncated == 0L
cat(sprintf("  cells where NEITHER arm truncated (the honest reach comparison): %d/%d\n",
            sum(clean), length(clean)))
if (sum(d$truncated) > sum(b$truncated))
  cat("  NOTE: deltas truncated more often than base -- on those cells a reach gap is a\n",
      "        BUDGET artefact (deltas cost ~15x candidates/rep), not a config failure.\n")

# The verdict is computed on NON-TRUNCATED cells only: a truncated arm never converged,
# so scoring its reach would repeat the arm-B error of reading a budget cut as a result.
reachB <- mean(b$reached[clean]); reachD <- mean(d$reached[clean])
tierBad <- character(0)
for (tr in unique(P$tier)) {
  ib <- clean & b$tier == tr; id <- clean & d$tier == tr
  if (sum(ib) && mean(d$reached[id]) < mean(b$reached[ib])) tierBad <- c(tierBad, tr)
}
cat(sprintf("\n=== PRE-REGISTERED VERDICT (non-truncated cells, n = %d) ===\n", sum(clean)))
cat(sprintf("  reach base=%.3f deltas=%.3f; tier regressions: %s\n  --> %s\n",
            reachB, reachD, if (length(tierBad)) paste(tierBad, collapse = ",") else "none",
            if (sum(clean) < 0.5 * length(clean))
              "INCONCLUSIVE -- too few non-truncated cells; re-run with larger caps"
            else if (reachD >= reachB && !length(tierBad)) "SHIP v1"
            else "NO-SHIP (restrict or drop)"))
cat("\n(All-cells reach, for reference only -- confounded by truncation: ")
cat(sprintf("base=%.3f deltas=%.3f)\n", mean(b$reached), mean(d$reached)))
