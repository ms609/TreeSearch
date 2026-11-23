#!/usr/bin/env Rscript
# Sector-godrift isolation gate (route B). Full MaximizeParsimony "thorough"
# search with driftCycles=0 (main-tree drift OFF -- verified: gd=0 => 0 drift
# calls), so the ONLY drift is the in-sector godrift solve. Three arms:
#   none  : sectorGoDrift=0  (no sector drift at all -- does sector drift pay?)
#   union : sectorGoDrift=25, union-of-finals scorer (current default)
#   exact : sectorGoDrift=25, TS_DRIFT_EXACT exact directional scorer
# maxReplicates is the replicate-like policy knob; wall is the eval metric
# (policy-in-replicates-not-seconds). Same seed per arm => matched RNG. Writes
# (dataset, seed, reps, arm, wall_s, score) so the analysis reads score at
# matched wall.
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})
seeds   <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "15")))
reps    <- as.integer(strsplit(Sys.getenv("GATE_REPS", "1,2,4"), ",")[[1]])
dnames  <- strsplit(Sys.getenv("GATE_DATA", "Zhu2013,Zanol2014,Dikow2009"), ",")[[1]]
outfile <- Sys.getenv("GATE_OUT", "dev/profiling/drift-exactness-sector-bench.csv")
godrift <- as.integer(Sys.getenv("GATE_GODRIFT", "25"))

fitchPhy <- function(p){m<-PhyDatToMatrix(p,ambigNA=FALSE); m[m=="-"]<-"?"; MatrixToPhyDat(m)}
loadPhy <- function(nm){
  for (p in c(file.path("data-raw", paste0(nm, ".nex")),
              file.path("dev/benchmarks", paste0(nm, ".nex")),
              file.path("dev/benchmarks/data", paste0(nm, ".nex"))))
    if (file.exists(p)) return(fitchPhy(ReadAsPhyDat(p)))
  stop("dataset not found: ", nm)
}

rows <- list()
for (nm in dnames) {
  phy <- loadPhy(nm); nt <- length(phy)
  for (s in seeds) for (R in reps) {
    for (arm in c("none","union","exact")) {
      gd <- if (arm == "none") 0L else godrift
      if (arm == "exact") Sys.setenv(TS_DRIFT_EXACT="1") else Sys.unsetenv("TS_DRIFT_EXACT")
      set.seed(s)
      t0 <- Sys.time()
      res <- MaximizeParsimony(phy, strategy="thorough", driftCycles=0L,
                               sectorGoDrift=gd, maxReplicates=R,
                               nThreads=1L, verbosity=0L)
      wall <- as.double(difftime(Sys.time(), t0, units="secs"))
      Sys.unsetenv("TS_DRIFT_EXACT")
      sc <- as.double(attr(res, "score"))
      rows[[length(rows)+1L]] <- data.frame(dataset=nm, nTip=nt, seed=s, reps=R,
                                            arm=arm, wall_s=round(wall,3), score=sc)
      cat(sprintf("%-12s seed=%d reps=%d %-5s score=%.0f wall=%.2fs\n",
                  nm, s, R, arm, sc, wall))
    }
  }
}
df <- do.call(rbind, rows)
write.csv(df, outfile, row.names=FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(df), outfile))
