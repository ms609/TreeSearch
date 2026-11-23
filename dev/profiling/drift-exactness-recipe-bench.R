#!/usr/bin/env Rscript
# Recipe-level whole-search test (route B'): does the drift SCORER choice change
# a wall-limited thorough search on TRAINING-split MorphoBank data? One catalogue
# key per invocation (GATE_KEY), EW-recoded (-> Fitch). Three arms, sector drift
# OFF (sectorGoDrift=0) to isolate the MAIN-tree drift phase:
#   none  : driftCycles=0 (drift off; budget goes to replicates/other phases)
#   union : driftCycles=3, union-of-finals scorer (current default)
#   exact : driftCycles=3, TS_DRIFT_EXACT exact directional scorer
# maxReplicates is the wall knob; same seed per arm => matched RNG. TRAINING
# split only (validation-set-sequestered).
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})
key     <- Sys.getenv("GATE_KEY")
catpath <- Sys.getenv("GATE_CAT", "/nobackup/pjjg18/floor/mbank_catalogue.csv")
mdir    <- Sys.getenv("GATE_MATRICES", "")   # resolved below if empty
seeds   <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "6")))
reps    <- as.integer(strsplit(Sys.getenv("GATE_REPS", "2,4,8"), ",")[[1]])
outfile <- Sys.getenv("GATE_OUT", sprintf("dev/profiling/recipe_%s.csv", key))
driftN  <- as.integer(Sys.getenv("GATE_DRIFT", "3"))

cat <- read.csv(catpath, stringsAsFactors = FALSE)
row <- cat[cat$key == key, ]
if (nrow(row) != 1) stop("key not found / ambiguous: ", key)
if (row$split != "training") stop("REFUSING non-training split for tuning: ", key, " is ", row$split)
if (mdir == "") for (d in c("../neotrans/inst/matrices","../neotrans/matrices"))
  if (dir.exists(d)) { mdir <- normalizePath(d); break }
nex <- file.path(mdir, row$filename)
if (!file.exists(nex)) stop("matrix not found: ", nex)

fitchPhy <- function(p){m<-PhyDatToMatrix(p,ambigNA=FALSE); m[m=="-"]<-"?"; MatrixToPhyDat(m)}
phy <- fitchPhy(ReadAsPhyDat(nex)); nt <- length(phy)
cat(sprintf("KEY %s : %d tips (recoded EW), catalogue ntax=%d pct_inapp=%.1f\n",
            key, nt, row$ntax, row$pct_inapp))

rows <- list()
for (s in seeds) for (R in reps) for (arm in c("none","union","exact")) {
  dc <- if (arm == "none") 0L else driftN
  if (arm == "exact") Sys.setenv(TS_DRIFT_EXACT="1") else Sys.unsetenv("TS_DRIFT_EXACT")
  set.seed(s); t0 <- Sys.time()
  res <- MaximizeParsimony(phy, strategy="thorough", driftCycles=dc, sectorGoDrift=0L,
                           maxReplicates=R, nThreads=1L, verbosity=0L)
  wall <- as.double(difftime(Sys.time(), t0, units="secs"))
  Sys.unsetenv("TS_DRIFT_EXACT")
  sc <- as.double(attr(res, "score"))
  rows[[length(rows)+1L]] <- data.frame(key=key, nTip=nt, seed=s, reps=R, arm=arm,
                                        wall_s=round(wall,3), score=sc)
  cat(sprintf("%-16s seed=%d reps=%d %-5s score=%.0f wall=%.2fs\n", key, s, R, arm, sc, wall))
}
write.csv(do.call(rbind, rows), outfile, row.names=FALSE)
cat(sprintf("Wrote %s\n", outfile))
