#!/usr/bin/env Rscript
# T-364/T-370 battery, step 4: ANALYSIS.
#
# Reporting rules this script exists to enforce (they are house rules, learned the
# hard way, not stylistic preferences):
#  * MEDIAN + SIGN COUNT + count of cells >10% slower.  Never the mean: a -4.9s
#    MEAN once concealed a +2.1% MEDIAN regression with 24 of 44 matrices slower
#    (memory wall-median-not-mean).
#  * The pairing unit is the MATRIX, not the (matrix, seed) cell.  Signed-rank over
#    cells is pseudo-replicated -- one matrix once produced p = 0.00067 where the
#    matrix-level answer was p = 0.98 (memory pair-on-matrices-not-seeds).  So
#    seeds are collapsed to a per-matrix median FIRST, and the sign test runs over
#    matrices.
#  * Censored cells (maxSeconds reached) are excluded from ratio statistics and
#    reported separately: their wall is a setting, not a measurement.
#
# Env: BATT_DIR (meta_*.csv + events_*.csv), OUT (report path).

batt <- Sys.getenv("BATT_DIR", "batt")
out <- Sys.getenv("OUT", "t364_report.md")

meta_files <- list.files(batt, pattern = "^meta_.*\\.csv$", full.names = TRUE)
if (!length(meta_files)) stop("no meta_*.csv in ", batt)
M <- do.call(rbind, lapply(meta_files, read.csv, stringsAsFactors = FALSE))
M$arm <- as.character(M$arm)

cat(sprintf("Loaded %d arm-runs over %d matrices, %d shapes, %d seeds\n",
            nrow(M), length(unique(M$key)), length(unique(M$shape)),
            length(unique(M$seed))))

# ---- Completeness: every cell must have all three arms, or the pairing is a lie.
cellid <- paste(M$key, M$shape, M$seed, sep = "|")
tab <- table(cellid)
incomplete <- names(tab)[tab != 3L]
if (length(incomplete)) {
  cat(sprintf("WARNING: %d cells lack all 3 arms; dropped:\n  %s\n",
              length(incomplete), paste(head(incomplete, 10), collapse = "\n  ")))
  M <- M[!cellid %in% incomplete, ]
}

# ---- Library discrimination sanity: three identical libs would give a clean and
# ---- meaningless null.  Compliance must separate arm 1 from arms 2/3.
comp <- aggregate(cbind(n_compliant, n_mpt, n_canonical) ~ arm, M, sum)
comp$pct_compliant <- 100 * comp$n_compliant / comp$n_mpt
comp$pct_canonical <- 100 * comp$n_canonical / comp$n_mpt
cat("\n-- returned-tree compliance by arm (all trees, complement-aware) --\n")
print(comp[, c("arm", "n_mpt", "pct_compliant", "pct_canonical")], row.names = FALSE)

cens <- aggregate(censored ~ arm, M, sum)
cat("\n-- censored (maxSeconds reached) runs by arm --\n")
print(cens, row.names = FALSE)

# ---- Collapse seeds to a per-(key, shape, arm) median.
key_med <- function(d, val) {
  a <- aggregate(d[[val]], by = list(key = d$key, shape = d$shape, arm = d$arm),
                 FUN = stats::median)
  names(a)[4] <- val
  a
}
wide <- function(a, val) {
  w <- reshape(a, idvar = c("key", "shape"), timevar = "arm", direction = "wide")
  names(w) <- sub(paste0(val, "\\."), "arm", names(w))
  w
}

# Censored cells poison a wall ratio; drop them from wall stats only.
Mw <- M[!M$censored, ]
cellid_w <- paste(Mw$key, Mw$shape, Mw$seed, sep = "|")
tw <- table(cellid_w)
Mw <- Mw[cellid_w %in% names(tw)[tw == 3L], ]

W <- wide(key_med(Mw, "wall_s"), "wall_s")
S <- wide(key_med(M, "score"), "score")
R <- wide(key_med(M, "reps_done"), "reps_done")

# A run that returned a constraint-VIOLATING tree answered an easier question, so
# its score is not commensurable with a compliant arm's -- and it can be
# *spuriously better* (a violating tree may score below the constrained optimum,
# which is the tell T-390/T-391 were caught by).  So the score comparison is run
# twice: over everything, and restricted to cells where ALL THREE arms returned
# fully compliant trees.
M$fully_compliant <- M$n_compliant == M$n_mpt
cellid2 <- paste(M$key, M$shape, M$seed, sep = "|")
ok_cells <- names(which(tapply(M$fully_compliant, cellid2, all)))
Mc <- M[cellid2 %in% ok_cells, ]
cat(sprintf("\nFully-compliant-in-all-arms cells: %d of %d\n",
            length(ok_cells), length(unique(cellid2))))
bad <- aggregate(fully_compliant ~ arm, M, function(x) sum(!x))
names(bad)[2] <- "n_runs_with_a_violating_tree"
print(bad, row.names = FALSE)
Sc <- wide(key_med(Mc, "score"), "score")

pair_report <- function(W, num, den, label, strata) {
  lines <- character(0)
  for (st in strata) {
    sub <- if (identical(st, "ALL")) W else W[W$shape == st, ]
    a <- sub[[paste0("arm", num)]]
    b <- sub[[paste0("arm", den)]]
    ok <- is.finite(a) & is.finite(b) & b > 0
    a <- a[ok]; b <- b[ok]
    if (!length(a)) next
    ratio <- a / b
    n <- length(ratio)
    slower <- sum(ratio > 1)
    faster <- sum(ratio < 1)
    tie <- sum(ratio == 1)
    over10 <- sum(ratio > 1.10)
    p <- if (slower + faster > 0) {
      stats::binom.test(slower, slower + faster, 0.5)$p.value
    } else NA_real_
    lines <- c(lines, sprintf(
      "| %s | %s | %d | %.3f | %.3f | %d / %d / %d | %d | %s |",
      label, st, n, stats::median(ratio),
      stats::median(a) / stats::median(b),
      slower, faster, tie, over10,
      if (is.na(p)) "-" else formatC(p, format = "g", digits = 2)))
  }
  lines
}

strata <- c(sort(unique(W$shape)), "ALL")

# Timer-resolution honesty: a matrix whose whole search runs in a few clock ticks
# contributes a ratio built out of quantisation noise.  Report how many, so the
# sign counts can be read with that in mind rather than discovering it later.
RES_S <- 0.05
n_coarse <- sum(pmin(W$arm1, W$arm2, W$arm3, na.rm = TRUE) < RES_S, na.rm = TRUE)
cat(sprintf("\nMatrices whose fastest arm ran in < %.2fs (timer-resolution noise): %d of %d\n",
            RES_S, n_coarse, nrow(W)))

cat("\n-- WALL: arm ratios, per-matrix medians, matrix-level sign test --\n")
cat("| pair | stratum | nMatrix | medianRatio | medianOfMedians | slower/faster/tie | n>10% slower | signTest p |\n")
cat("|---|---|---|---|---|---|---|---|\n")
wall_lines <- c(
  pair_report(W, 2, 1, "arm2/arm1 (complement only vs pre-fix)", strata),
  pair_report(W, 3, 1, "arm3/arm1 (MERGED vs pre-fix)", strata),
  pair_report(W, 3, 2, "arm3/arm2 (reroot removes the cost)", strata))
cat(paste(wall_lines, collapse = "\n"), "\n")

score_report <- function(S, num, den, label, strata) {
  lines <- character(0)
  for (st in strata) {
    sub <- if (identical(st, "ALL")) S else S[S$shape == st, ]
    a <- sub[[paste0("arm", num)]]
    b <- sub[[paste0("arm", den)]]
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    if (!length(a)) next
    d <- a - b
    worse <- sum(d > 0)
    better <- sum(d < 0)
    tie <- sum(d == 0)
    p <- if (worse + better > 0) {
      stats::binom.test(worse, worse + better, 0.5)$p.value
    } else NA_real_
    lines <- c(lines, sprintf("| %s | %s | %d | %+.1f | %d / %d / %d | %s |",
      label, st, length(d), stats::median(d), worse, better, tie,
      if (is.na(p)) "-" else formatC(p, format = "g", digits = 2)))
  }
  lines
}
cat("\n-- SCORE (parsimony; lower is better): arm deltas --\n")
cat("| pair | stratum | nMatrix | medianDelta | worse/better/tie | signTest p |\n")
cat("|---|---|---|---|---|---|\n")
score_lines <- c(
  score_report(S, 2, 1, "arm2 - arm1", strata),
  score_report(S, 3, 1, "arm3 - arm1 (MERGED vs pre-fix)", strata),
  score_report(S, 3, 2, "arm3 - arm2", strata))
cat(paste(score_lines, collapse = "\n"), "\n")

cat("\n-- SCORE, restricted to cells fully compliant in ALL arms --\n")
cat("| pair | stratum | nMatrix | medianDelta | worse/better/tie | signTest p |\n")
cat("|---|---|---|---|---|---|\n")
cat(paste(c(
  score_report(Sc, 2, 1, "arm2 - arm1 [compliant only]", strata),
  score_report(Sc, 3, 1, "arm3 - arm1 [compliant only]", strata),
  score_report(Sc, 3, 2, "arm3 - arm2 [compliant only]", strata)),
  collapse = "\n"), "\n")

cat("\n-- REPLICATES RUN before the convergence rule stopped the search --\n")
rep_lines <- character(0)
for (st in strata) {
  sub <- if (identical(st, "ALL")) R else R[R$shape == st, ]
  rep_lines <- c(rep_lines, sprintf("| %s | %.1f | %.1f | %.1f |", st,
    stats::median(sub$arm1, na.rm = TRUE), stats::median(sub$arm2, na.rm = TRUE),
    stats::median(sub$arm3, na.rm = TRUE)))
}
cat("| stratum | arm1 | arm2 | arm3 |\n|---|---|---|---|\n")
cat(paste(rep_lines, collapse = "\n"), "\n")

# ---- Per-replicate wall: the CONDITIONAL claim.  If arm 2's cost is the
# ---- complement-rooted start, its per-replicate wall is a MIXTURE -- most
# ---- replicates normal, a minority much slower -- not a uniform inflation.
ev_files <- list.files(batt, pattern = "^events_.*\\.csv$", full.names = TRUE)
if (length(ev_files)) {
  E <- do.call(rbind, lapply(ev_files, read.csv, stringsAsFactors = FALSE))
  E <- E[E$phase == "replicate" & is.finite(E$wall_s), ]
  if (nrow(E)) {
    E <- E[order(E$task, E$arm, E$replicate), ]
    grp <- paste(E$task, E$arm)
    E$dt <- c(NA, diff(E$wall_s))
    E$dt[c(TRUE, grp[-1] != grp[-length(grp)])] <- NA
    E <- E[is.finite(E$dt) & E$dt >= 0, ]
    cat("\n-- per-replicate wall (s) by arm: median, p90, p99, max, ratio p90/median --\n")
    cat("| arm | nRep | median | p90 | p99 | max | p90/median |\n|---|---|---|---|---|---|---|\n")
    for (a in sort(unique(E$arm))) {
      x <- E$dt[E$arm == a]
      q <- stats::quantile(x, c(0.5, 0.9, 0.99), na.rm = TRUE)
      # A zero median means the clock cannot resolve one replicate on this
      # hardware (Windows proc.time() steps in ~15.6 ms), not that replicates are
      # free; print "-" rather than an Inf that reads like a finding.
      rat <- if (q[[1]] > 0) sprintf("%.2f", q[[2]] / q[[1]]) else "-"
      cat(sprintf("| %s | %d | %.4f | %.4f | %.4f | %.4f | %s |\n",
                  a, length(x), q[[1]], q[[2]], q[[3]], max(x), rat))
    }
    write.csv(E, file.path(batt, "per_replicate.csv"), row.names = FALSE)
  }
}

# ---- Per-PHASE timings: test the recorded T-384 signature on its own terms.
# The recorded claim is "every phase slower on a score-IDENTICAL trajectory =
# moves being rejected".  But cumulative phase totals are summed over however many
# replicates ran, so if a blocked arm runs MORE replicates its totals rise even
# when each replicate is cheaper.  Rejected moves should make a replicate do LESS
# work, not more.  So print both, and let the per-replicate column decide.
tm_files <- list.files(batt, pattern = "^timings_.*\\.csv$", full.names = TRUE)
if (length(tm_files)) {
  TM <- do.call(rbind, lapply(tm_files, read.csv, stringsAsFactors = FALSE))
  TM$arm <- as.character(TM$arm)
  TM <- TM[is.finite(TM$ms) & is.finite(TM$reps_done) & TM$reps_done > 0, ]
  TM$per_rep <- TM$ms / TM$reps_done
  live <- unique(TM$phase[TM$ms > 0])
  TM <- TM[TM$phase %in% live, ]
  phase_tab <- function(col, title) {
    a <- aggregate(TM[[col]], by = list(phase = TM$phase, arm = TM$arm),
                   FUN = stats::median)
    names(a)[3] <- "v"
    w <- reshape(a, idvar = "phase", timevar = "arm", direction = "wide")
    names(w) <- sub("v\\.", "arm", names(w))
    w <- w[order(-w$arm2), ]
    cat(sprintf("\n-- %s --\n", title))
    cat("| phase | arm1 | arm2 | arm3 | a2/a1 | a3/a1 |\n|---|---|---|---|---|---|\n")
    for (i in seq_len(nrow(w))) {
      r <- w[i, ]
      cat(sprintf("| %s | %.4g | %.4g | %.4g | %.3f | %.3f |\n", r$phase,
                  r$arm1, r$arm2, r$arm3, r$arm2 / r$arm1, r$arm3 / r$arm1))
    }
  }
  phase_tab("ms", "per-phase CUMULATIVE ms (median over cells)")
  phase_tab("per_rep", "per-phase PER-REPLICATE ms (cumulative / reps_done)")
  write.csv(TM, file.path(batt, "per_phase.csv"), row.names = FALSE)
}

# ---- Budget exhaustion: the budget-INDEPENDENT statement of arm 2's cost.
# A wall ratio partly reflects maxReplicates, which is a chosen setting; "did the
# search recognise convergence at all" does not.
ex <- M
ex$exhausted <- ex$reps_done >= ex$maxrep
cat("\n-- replicate-budget exhaustion (reps_done >= maxReplicates) --\n")
cat("| arm | stratum | nRuns | nExhausted | pct |\n|---|---|---|---|---|\n")
for (a in sort(unique(ex$arm))) {
  for (st in c(sort(unique(ex$shape)), "ALL")) {
    sub <- if (identical(st, "ALL")) ex[ex$arm == a, ] else ex[ex$arm == a & ex$shape == st, ]
    if (!nrow(sub)) next
    cat(sprintf("| %s | %s | %d | %d | %.1f%% |\n", a, st, nrow(sub),
                sum(sub$exhausted), 100 * mean(sub$exhausted)))
  }
}

# ---- Per-matrix gain table: never ship only the summary (pair-on-matrices).
per_key <- merge(W, S, by = c("key", "shape"), suffixes = c(".wall", ".score"))
per_key <- merge(per_key, R, by = c("key", "shape"))
nt <- unique(M[, c("key", "n_tip", "tier")])
per_key <- merge(per_key, nt, by = "key")
per_key <- per_key[order(per_key$n_tip, per_key$shape), ]
write.csv(per_key, file.path(batt, "per_matrix.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s and per_replicate.csv\n", file.path(batt, "per_matrix.csv")))

cat("\n-- per-matrix wall medians (s) and ratios --\n")
cat("| key | nTip | shape | arm1 | arm2 | arm3 | a2/a1 | a3/a1 | score a1/a2/a3 |\n")
cat("|---|---|---|---|---|---|---|---|---|\n")
for (i in seq_len(nrow(per_key))) {
  r <- per_key[i, ]
  cat(sprintf("| %s | %d | %s | %.3f | %.3f | %.3f | %.2f | %.2f | %.0f/%.0f/%.0f |\n",
              r$key, r$n_tip, r$shape, r$arm1.wall, r$arm2.wall, r$arm3.wall,
              r$arm2.wall / r$arm1.wall, r$arm3.wall / r$arm1.wall,
              r$arm1.score, r$arm2.score, r$arm3.score))
}
