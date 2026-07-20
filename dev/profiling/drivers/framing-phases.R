# Phase-level attribution of the FULL EW MaximizeParsimony run (the mission
# workload), using the per-phase timings already accumulated in C++
# (result$timings, ms).  No rebuild: trace() captures the aggregated breakdown
# from the existing release DLL named by TS_LIB.  Drivers: Zhu/Zanol fitch.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Zhu2013 Zanol2014")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])
reps <- as.integer(Sys.getenv("TS_REPS", "3"))

# Capture the aggregated timings via the underlying Rcpp call's return value.
# returnValue() is the canonical trace-exit accessor; ts_driven_search's return
# list carries $timings before MaximizeParsimony post-processes it.
trace("ts_driven_search", where = asNamespace("TreeSearch"),
      exit = quote(assign("TS_TIMINGS", returnValue()$timings, envir = .GlobalEnv)),
      print = FALSE)

cat(sprintf("LIB=%s reps=%d\n", Sys.getenv("TS_LIB"), reps))
agg <- list()
for (nm in dsN) {
  phy <- fc(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    set.seed(sd); t0 <- Sys.time()
    r <- suppressWarnings(MaximizeParsimony(
      phy, maxReplicates = reps, nThreads = 1L, strategy = "auto", verbosity = 0L))
    w <- as.double(difftime(Sys.time(), t0, units = "secs"))
    tm <- unlist(get("TS_TIMINGS", envir = .GlobalEnv))
    cat(sprintf("\n=== %s seed%d  score=%s  wall=%.2fs  sum(timings)=%.0fms ===\n",
                nm, sd, attr(r, "score"), w, sum(tm, na.rm = TRUE)))
    o <- sort(tm[tm > 0], decreasing = TRUE)
    pct <- round(100 * o / sum(tm, na.rm = TRUE), 1)
    print(data.frame(phase = names(o), ms = round(as.numeric(o)), pct = as.numeric(pct)),
          row.names = FALSE)
    agg[[paste(nm, sd)]] <- tm
  }
}
cat("\nPHASES COMPLETE\n")
