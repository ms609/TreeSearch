#!/usr/bin/env Rscript
# Matched-wall A/B gate for the exact-directional EW drift scorer.
#
# Per (dataset, seed): build ONE local-optimum start (ts_tbr_search, scorer-
# independent), then run drift_search from that SAME start under the union
# (default) and exact (TS_DRIFT_EXACT) scorers with the SAME RNG stream
# (set.seed before each), for a sweep of nCycles (the replicate-like policy knob
# -- wall is the eval metric, per memory policy-in-replicates-not-seconds).
# Records wall_s + final score per run to a CSV; analysis compares score at
# matched wall between the two scorers.  Decision rule: flip default only if no
# time-to-optimum regression on ALL three datasets, else land opt-in.
#
# Env knobs (small local defaults; override for the Hamilton array):
#   GATE_SEEDS   (default 3)      number of seeds
#   GATE_NCYC    (default 8,24)   comma-separated nCycles sweep
#   GATE_DATA    (default all 3)  comma-separated dataset names
#   GATE_LIB     (default .agent-hj)  package library
#   GATE_OUT     (default dev/profiling/drift-exactness-gate-bench.csv)
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})

seeds  <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "3")))
ncyc   <- as.integer(strsplit(Sys.getenv("GATE_NCYC", "8,24"), ",")[[1]])
dnames <- strsplit(Sys.getenv("GATE_DATA", "Zhu2013,Zanol2014,Dikow2009"), ",")[[1]]
outfile <- Sys.getenv("GATE_OUT", "dev/profiling/drift-exactness-gate-bench.csv")

fitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m)
}
makeDs <- function(dataset) {
  at <- attributes(dataset)
  list(contrast = at$contrast,
       tip_data = matrix(unlist(dataset, use.names = FALSE),
                         nrow = length(dataset), byrow = TRUE),
       weight = at$weight, levels = at$levels)
}

runDrift <- function(startEdge, ds, seed, exact, nCycles) {
  if (exact) Sys.setenv(TS_DRIFT_EXACT = "1") else Sys.unsetenv("TS_DRIFT_EXACT")
  set.seed(seed)                                   # identical RNG stream both paths
  t0 <- Sys.time()
  r <- TreeSearch:::ts_drift_search(startEdge, ds$contrast, ds$tip_data,
                                    ds$weight, ds$levels, nCycles = nCycles,
                                    afdLimit = 3L, rfdLimit = 0.1, maxHits = 1L)
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  Sys.unsetenv("TS_DRIFT_EXACT")
  list(score = r$score, wall = wall)
}

rows <- list()
for (nm in dnames) {
  phy <- fitchPhy(ReadAsPhyDat(file.path("data-raw", paste0(nm, ".nex"))))
  ds <- makeDs(phy); nTip <- length(phy)
  for (s in seeds) {
    set.seed(s)
    start <- RandomTree(nTip, root = TRUE); start$tip.label <- names(phy)
    tbr <- TreeSearch:::ts_tbr_search(start$edge, ds$contrast, ds$tip_data,
                                      ds$weight, ds$levels, maxHits = 5L)
    startEdge <- tbr$edge
    for (nc in ncyc) {
      for (ex in c(FALSE, TRUE)) {
        res <- runDrift(startEdge, ds, s, ex, nc)
        rows[[length(rows) + 1L]] <- data.frame(
          dataset = nm, nTip = nTip, seed = s, nCycles = nc,
          scorer = if (ex) "exact" else "union",
          localopt = tbr$score, score = res$score, wall_s = round(res$wall, 3))
        cat(sprintf("%-10s seed=%d nc=%3d %-5s  score=%.0f  wall=%.2fs\n",
                    nm, s, nc, if (ex) "exact" else "union", res$score, res$wall))
      }
    }
  }
}
df <- do.call(rbind, rows)
write.csv(df, outfile, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(df), outfile))

# Directional summary: per (dataset, nCycles) mean score & wall by scorer.
cat("\n== mean score / wall by dataset x nCycles x scorer ==\n")
agg <- aggregate(cbind(score, wall_s) ~ dataset + nCycles + scorer, df, mean)
agg <- agg[order(agg$dataset, agg$nCycles, agg$scorer), ]
print(agg, row.names = FALSE)
