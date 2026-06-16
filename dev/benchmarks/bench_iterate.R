# ITERATE tier — the pre-commit lever GATE, ~1-2 min, run POOL-DRAINED.
#
# Gap panel at a FIXED REPLICATE COUNT (NOT maxSeconds), nThreads=1, a few seeds.
# Fixed-replicate stopping is the only condition that makes candidates_evaluated
# machine-load-independent today (a true candidate-budget stop is the planned
# C++ refinement; see dev/plans). Reports per-dataset median candidates +
# median/best score, and a gap-to-TNT column from headtohead_phase0.csv targets.
#
# This is the signal a lever must move: a candidate-efficiency win shows as LOWER
# median candidates at equal-or-better score. ~0.7% seed spread on candidates
# (vs the +/-2-4 step score lottery), so 2-3 seeds resolve a real change.
#
# Env: TS_LIB, TS_DATASETS, TS_SEEDS (1 2 3), TS_REPS (20), OUT_CSV.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
dsN   <- strsplit(trimws(Sys.getenv("TS_DATASETS",
          "Wortley2006 Eklund2004 Zanol2014 Zhu2013 Giles2015 Dikow2009")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
reps  <- as.integer(Sys.getenv("TS_REPS", "20"))
out   <- Sys.getenv("OUT_CSV", "dev/benchmarks/iterate_latest.csv")
# TNT-best targets (apples-to-apples Fitch) from headtohead_phase0.csv.
tnt <- c(Wortley2006 = 479, Eklund2004 = 440, Zanol2014 = 1261,
         Zhu2013 = 624, Giles2015 = 670, Dikow2009 = 1606)

data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

t0 <- Sys.time()
rows <- list()
for (nm in dsN) {
  d <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    set.seed(sd)
    r <- suppressWarnings(MaximizeParsimony(d, maxReplicates = reps, targetHits = 999L,
                                            maxSeconds = 0, nThreads = 1L, verbosity = 0L))
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, seed = sd, score = attr(r, "score"),
      candidates = attr(r, "candidates_evaluated"), stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(res, res$dataset), function(d) {
  nm <- d$dataset[1]
  data.frame(dataset = nm, tips = length(inapplicable.phyData[[nm]]),
             score_best = min(d$score), score_med = median(d$score),
             gap = median(d$score) - (if (nm %in% names(tnt)) tnt[[nm]] else NA),
             cand_med = median(d$candidates),
             cand_spread_pct = round(100 * (max(d$candidates) - min(d$candidates)) /
                                       median(d$candidates), 2),
             stringsAsFactors = FALSE)
}))
agg <- agg[order(-agg$gap), ]
cat(sprintf("ITERATE | panel x %d seeds | %d reps | %.0fs\n", length(seeds), reps,
            as.double(difftime(Sys.time(), t0, units = "secs"))))
print(agg, row.names = FALSE)
write.csv(res, out, row.names = FALSE)
cat("rows ->", out, "\n")
