# Hamilton driver (task #29): full-search time-matched gate for rasStarts=3 in the
# AUTO-SELECTED `thorough` preset.  Runs the WHOLE thorough pipeline at matched
# wall-clock, varying ONLY rasStarts (explicit arg overrides the preset field).
# Authoritative wall-clock (vs the indicative local run diag_thorough_rasstarts_tm.R).
#
# Grid: datasets x rasStarts{1,3} x budgets{60,120}s x seeds{1..NSEED}.
# Decision: adopt rasStarts=3 in thorough iff it improves (or matches at lower
# variance) the median score at matched budget across datasets, without hurting
# replicate throughput enough to regress on any.
#
# Env: TS_LIB (installed pkg), OUTDIR, NSEED (default 10), BUDGETS, TS_DATASETS.
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("TS_LIB", .libPaths()[1]))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
target  <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
dsN     <- strsplit(trimws(Sys.getenv("TS_DATASETS",
             "Zanol2014 Zhu2013 Wortley2006 Giles2015")), "\\s+")[[1]]
budgets <- as.integer(strsplit(trimws(Sys.getenv("BUDGETS", "60 120")), "\\s+")[[1]])
nseed   <- as.integer(Sys.getenv("NSEED", "10"))
outdir  <- Sys.getenv("OUTDIR", "dev/benchmarks")
out_csv <- file.path(outdir, "thorough_rasstarts.csv")

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]]); tgt <- target[[nm]]; nt <- length(phy)
  for (secs in budgets) for (ras in c(1L, 3L)) for (s in seq_len(nseed)) {
    set.seed(s)
    t <- system.time(r <- suppressWarnings(MaximizeParsimony(phy,
           strategy = "thorough", rasStarts = ras, maxSeconds = secs,
           nThreads = 1L, verbosity = 0L)))
    sc <- min(as.double(attr(r, "score")))
    nrep <- length(as.double(attr(r, "score")))
    rows[[length(rows) + 1L]] <- data.frame(dataset = nm, nTip = nt,
      target = tgt, budget = secs, rasStarts = ras, seed = s,
      score = sc, over = sc - tgt, n_trees = nrep,
      elapsed = round(as.double(t["elapsed"]), 1))
    cat(sprintf("%-12s b=%3d ras=%d s=%2d -> %.0f (%+.0f)  [%.0fs]\n",
                nm, secs, ras, s, sc, sc - tgt, as.double(t["elapsed"])))
  }
}
df <- do.call(rbind, rows)
write.csv(df, out_csv, row.names = FALSE)

# Per (dataset,budget): median over vs target, by rasStarts.
cat("\n=== median (score - target) by dataset x budget x rasStarts ===\n")
agg <- aggregate(over ~ dataset + budget + rasStarts, df, median)
w <- reshape(agg, idvar = c("dataset", "budget"), timevar = "rasStarts",
             direction = "wide")
names(w) <- sub("over\\.", "ras", names(w))
w$delta <- w$ras3 - w$ras1   # negative = ras3 better
print(w, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
