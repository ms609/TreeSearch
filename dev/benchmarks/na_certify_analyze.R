# Analyse the NA-certification gate panel (bench_na_certify_cell.R partials).
#
# Rules this script exists to enforce, each of which has previously been got
# wrong in this project:
#
# * FLOOR ATTAINMENT is the quality metric -- the share of runs reaching the best
#   score any arm found for that matrix.  Never a rate whose numerator is
#   "improvements found": escapes/second once rated an arm that reached the best
#   tree 0% of the time a 1.57x win (reach-not-escape-rate).
# * ONE NUMBER PER MATRIX before any paired test.  (matrix x seed) pairing is
#   pseudo-replication and has manufactured p = 0.0007 from an effect that was
#   p = 0.98 at matrix level (pair-on-matrices-not-seeds).  So the sign test's n
#   is the number of MATRICES, and a per-matrix gain table is printed alongside.
# * WALL as median + sign count + how many matrices are >10% slower, never the
#   mean: a -4.9 s mean once hid a +2.1% median regression (wall-median-not-mean).
# * Firing evidence.  If nEvsSkipped == 0 the flag never reached a live call
#   site, which is a different finding from "it fired and bought nothing", and
#   the gate on the certifier (do_reroot) needs tabuSize == 0 -- which the
#   shipped presets do not set.  Reported per arm, per tabu.
#
# Usage: Rscript dev/benchmarks/na_certify_analyze.R [partialDir]

args <- commandArgs(trailingOnly = TRUE)
partdir <- if (length(args) >= 1L) args[[1]] else
  "dev/benchmarks/partials_na_certify"

files <- list.files(partdir, pattern = "^nacert_\\d+\\.csv$", full.names = TRUE)
if (!length(files)) stop("no partials in ", partdir)
dat <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
cat(sprintf("Cells: %d files, %d rows, %d matrices, %d seeds, tabu %s\n",
            length(files), nrow(dat), length(unique(dat$dataset)),
            length(unique(dat$seed)),
            paste(sort(unique(dat$tabu)), collapse = "/")))

baseArm <- "A_certify"
armOrder <- c("A_certify", "B_gate", "C_final", "B_gate_mw", "C_final_mw")
dat$arm <- factor(dat$arm, levels = intersect(armOrder, unique(dat$arm)))

signTest <- function(delta) {
  # Two-sided exact sign test on the nonzero differences.  Reported with the
  # counts, because at these panel sizes the counts are the finding and the
  # p-value is decoration.
  d <- delta[!is.na(delta) & delta != 0]
  if (!length(d)) return(list(nUp = 0L, nDown = 0L, p = NA_real_))
  nUp <- sum(d > 0); nDown <- sum(d < 0)
  list(nUp = nUp, nDown = nDown,
       p = stats::binom.test(nUp, length(d), 0.5)$p.value)
}

for (tb in sort(unique(dat$tabu))) {
  sub <- dat[dat$tabu == tb, ]
  cat(sprintf("\n\n================ tabuSize = %d ================\n", tb))
  cat(if (tb == 0L)
        "  do_reroot is LIVE at every whole-tree search (the `sprint` preset, and\n  the configuration the 97.7% profile was taken in).\n"
      else
        "  do_reroot requires tabu_size == 0, so only the sector sub-searches and\n  the fuse cleanup reach the certifier (the `default`/`thorough` presets).\n")

  # ---- one number per matrix ----
  best <- tapply(sub$score, sub$dataset, min)
  sub$isFloor <- sub$score <= best[sub$dataset] + 1e-8
  attain <- tapply(sub$isFloor, list(sub$dataset, sub$arm), mean)
  wallMed <- tapply(sub$wall, list(sub$dataset, sub$arm), stats::median)
  repsMed <- tapply(sub$reps_done, list(sub$dataset, sub$arm), stats::median)

  cat("\n-- Floor attainment (share of seeds reaching the matrix best) --\n")
  tab <- data.frame(matrix = rownames(attain), best = as.numeric(best[rownames(attain)]),
                    round(attain, 3), check.names = FALSE)
  print(tab, row.names = FALSE)

  cat("\n-- Wall, per-matrix median seconds --\n")
  print(data.frame(matrix = rownames(wallMed), round(wallMed, 2),
                   check.names = FALSE), row.names = FALSE)

  cat("\n-- Replicates completed, per-matrix median --\n")
  print(data.frame(matrix = rownames(repsMed), repsMed, check.names = FALSE),
        row.names = FALSE)

  cat("\n-- Certifier firing (summed over the whole arm) --\n")
  fire <- aggregate(cbind(nEvs, nEvsSkipped, nEvsImproved) ~ arm, sub, sum)
  print(fire, row.names = FALSE)
  cat("  nEvsSkipped == 0 for a gated arm => the flag never fired; a null\n",
      "  result from that arm says nothing about the lever.\n", sep = "")

  # ---- paired tests, matrix as the unit ----
  arms <- setdiff(colnames(attain), baseArm)
  cat("\n-- Paired vs ", baseArm, " (unit = MATRIX, n = ",
      nrow(attain), ") --\n", sep = "")
  res <- do.call(rbind, lapply(arms, function(a) {
    dAtt <- attain[, a] - attain[, baseArm]
    stAtt <- signTest(dAtt)
    ratio <- wallMed[, a] / wallMed[, baseArm]
    dWall <- wallMed[, a] - wallMed[, baseArm]
    stWall <- signTest(dWall)
    data.frame(
      arm = a,
      attainMeanDelta = round(mean(dAtt, na.rm = TRUE), 4),
      attainBetter = stAtt$nUp, attainWorse = stAtt$nDown,
      attainP = signif(stAtt$p, 3),
      wallMedRatio = round(stats::median(ratio, na.rm = TRUE), 3),
      wallSlower = stWall$nUp, wallFaster = stWall$nDown,
      wallP = signif(stWall$p, 3),
      nOver10pcSlower = sum(ratio > 1.1, na.rm = TRUE),
      stringsAsFactors = FALSE)
  }))
  print(res, row.names = FALSE)

  cat("\n-- Per-matrix floor-attainment gain table (matrices where arms differ) --\n")
  gain <- data.frame(matrix = rownames(attain),
                     round(attain[, arms, drop = FALSE] - attain[, baseArm], 3),
                     check.names = FALSE)
  changed <- rowSums(abs(as.matrix(gain[, -1, drop = FALSE])), na.rm = TRUE) > 0
  if (any(changed)) print(gain[changed, ], row.names = FALSE) else
    cat("  (none: every arm attained the same floor on every matrix)\n")
}

cat("\n\nVERDICT RULE: an arm ships only if floor attainment does not regress.\n")
cat("A wall win with an attainment loss is a quality-for-speed trade, and the\n")
cat("matched-wall arms (*_mw) are what decide whether the extra replicates buy\n")
cat("the lost reach back.\n")
