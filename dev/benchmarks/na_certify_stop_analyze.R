# Analyse the SHIPPED-STOPPING-RULES panel (bench_na_certify_stop_cell.R).
#
# Same rules as na_certify_analyze.R, and for the same reasons:
#   * FLOOR ATTAINMENT, never a rate whose numerator is improvements found.
#   * ONE NUMBER PER MATRIX before any paired test -- (matrix x seed) pairing is
#     pseudo-replication and has manufactured p = 0.0007 from an effect that was
#     p = 0.98 at matrix level.
#   * WALL as median + sign count + count >10% slower, never the mean.
#   * Firing evidence per arm, so a null cannot be read as "the lever is
#     worthless" when it means "the flag never fired".
#
# What is different here: the budget is NOT fixed.  Every arm runs under the
# shipped stopping rules, so `reps` and `hitsToBest` are OUTCOMES, and "where did
# it stop, and why" is half the finding.  An arm that scores the same while
# stopping at a quarter of the wall is a win that a score-only table hides.
#
# Usage: Rscript dev/benchmarks/na_certify_stop_analyze.R [partialDir]

args <- commandArgs(trailingOnly = TRUE)
partdir <- if (length(args) >= 1L) args[[1]] else
  "dev/benchmarks/na_certify_stop_partials"

files <- list.files(partdir, pattern = "^nastop_\\d+\\.csv$", full.names = TRUE)
if (!length(files)) stop("no partials in ", partdir)
dat <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))

baseArm <- "A_default"
armOrder <- c("A_default", "B_gate", "C_final", "B_gate_reps",
              "A_hits3", "B_gate_hits3")
dat$arm <- factor(dat$arm, levels = intersect(armOrder, unique(dat$arm)))

cat(sprintf("Cells: %d files, %d rows, %d matrices, %d seeds, %d arms\n",
            length(files), nrow(dat), length(unique(dat$dataset)),
            length(unique(dat$seed)), nlevels(dat$arm)))

signTest <- function(delta) {
  d <- delta[!is.na(delta) & delta != 0]
  if (!length(d)) return(list(nUp = 0L, nDown = 0L, p = NA_real_))
  nUp <- sum(d > 0); nDown <- sum(d < 0)
  list(nUp = nUp, nDown = nDown,
       p = stats::binom.test(nUp, length(d), 0.5)$p.value)
}

# ---- one number per matrix ----
best <- tapply(dat$score, dat$dataset, min)
dat$isFloor <- dat$score <= best[dat$dataset] + 1e-8
attain  <- tapply(dat$isFloor, list(dat$dataset, dat$arm), mean)
wallMed <- tapply(dat$wall,    list(dat$dataset, dat$arm), stats::median)
repsMed <- tapply(dat$reps,    list(dat$dataset, dat$arm), stats::median)
hitsMed <- tapply(dat$hitsToBest, list(dat$dataset, dat$arm), stats::median)
ntax    <- tapply(dat$ntax,    dat$dataset, max)

cat("\n-- Floor attainment (share of seeds reaching the matrix best) --\n")
print(data.frame(matrix = rownames(attain), ntax = as.integer(ntax[rownames(attain)]),
                 best = as.numeric(best[rownames(attain)]), round(attain, 2),
                 check.names = FALSE), row.names = FALSE)

cat("\n-- Wall, per-matrix median seconds --\n")
print(data.frame(matrix = rownames(wallMed), round(wallMed, 1),
                 check.names = FALSE), row.names = FALSE)

cat("\n-- Replicates completed, per-matrix median --\n")
print(data.frame(matrix = rownames(repsMed), repsMed, check.names = FALSE),
      row.names = FALSE)

# ---- which rule actually ended each run? ----
# The whole reason for this panel.  `targetHits` can only bind if a run
# accumulates that many hits before the replicate cap stops it; on hard matrices
# it never does, and raising it is then inert.
cat("\n-- What stopped the run? (share of cells) --\n")
dat$capBound <- dat$reps >= dat$maxReplicates
stopTab <- do.call(rbind, lapply(split(dat, dat$arm), function(x) data.frame(
  arm = x$arm[1],
  hitCapBound = round(mean(x$capBound), 3),
  medHitsToBest = stats::median(x$hitsToBest),
  medTargetHits = stats::median(x$targetHits),
  perturbStop = round(mean(x$perturbStop), 3))))
print(stopTab, row.names = FALSE)
cat("  hitCapBound = share of runs that hit `maxReplicates`; for those,\n",
    "  `targetHits` was never reached and raising it changes nothing.\n", sep = "")

cat("\n-- Certifier firing (summed over the whole arm) --\n")
print(aggregate(cbind(nEvs, nEvsSkipped, nEvsImproved) ~ arm, dat, sum),
      row.names = FALSE)

# ---- paired tests, matrix as the unit ----
arms <- setdiff(colnames(attain), baseArm)
cat("\n-- Paired vs ", baseArm, " (unit = MATRIX, n = ", nrow(attain), ") --\n",
    sep = "")
res <- do.call(rbind, lapply(arms, function(a) {
  dAtt <- attain[, a] - attain[, baseArm]
  stAtt <- signTest(dAtt)
  ratio <- wallMed[, a] / wallMed[, baseArm]
  stWall <- signTest(wallMed[, a] - wallMed[, baseArm])
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

cat("\n-- Per-matrix floor-attainment gain (matrices where arms differ) --\n")
gain <- data.frame(matrix = rownames(attain),
                   round(attain[, arms, drop = FALSE] - attain[, baseArm], 2),
                   check.names = FALSE)
changed <- rowSums(abs(as.matrix(gain[, -1, drop = FALSE])), na.rm = TRUE) > 0
if (any(changed)) print(gain[changed, ], row.names = FALSE) else
  cat("  (none)\n")

# ---- does raising targetHits do anything at all? ----
cat("\n-- targetHits x3: did it change the run? --\n")
for (pair in list(c("A_default", "A_hits3"), c("B_gate", "B_gate_hits3"))) {
  w <- merge(dat[dat$arm == pair[1], c("dataset", "seed", "score", "reps", "wall")],
             dat[dat$arm == pair[2], c("dataset", "seed", "score", "reps", "wall")],
             by = c("dataset", "seed"), suffixes = c(".base", ".x3"))
  same <- w$score.base == w$score.x3 & w$reps.base == w$reps.x3
  cat(sprintf("  %-13s -> %-13s : identical (score AND reps) in %d of %d cells; median wall ratio %.2f\n",
              pair[1], pair[2], sum(same), nrow(w),
              stats::median(w$wall.x3 / w$wall.base)))
}

cat("\nVERDICT RULE: an arm ships only if floor attainment does not regress.\n")
