#!/usr/bin/env Rscript
# Route-3 marginal-value experiment (advisor design). Question: in a wall-limited
# search, is a wall-second better spent on (cheap, exact) drift than on more of
# the incumbent diversifier (ratchet)? If yes, drift earns a place in fast/
# wall-limited recipes that currently skip it -> worth re-tuning.
#
# From a SHARED per-seed local optimum (ts_tbr_search to convergence; scorer-
# independent), three continuations, each swept over a cycles knob (wall is the
# eval metric, per policy-in-replicates-not-seconds):
#   ratchet     : ts_ratchet_search        (the incumbent marginal spend)
#   drift_exact : ts_drift_search + TS_DRIFT_EXACT
#   drift_union : ts_drift_search (union)  (control: shows the unlock is exact's
#                                           cheapness, not drift per se)
# Decision: if drift_exact's anytime frontier beats ratchet's over some wall
# range AND dominates drift_union -> install exact-drift in wall-limited recipes.
# If ratchet dominates everywhere -> route-3 no-go.
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})
seeds  <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "20")))
cyc    <- as.integer(strsplit(Sys.getenv("GATE_CYC", "1,2,4,8,16"), ",")[[1]])
dnames <- strsplit(Sys.getenv("GATE_DATA", "Zhu2013,Zanol2014,Dikow2009"), ",")[[1]]
outfile <- Sys.getenv("GATE_OUT", "dev/profiling/drift-exactness-route3-bench.csv")

fitchPhy <- function(p){m<-PhyDatToMatrix(p,ambigNA=FALSE); m[m=="-"]<-"?"; MatrixToPhyDat(m)}
loadPhy <- function(nm){
  for (p in c(file.path("data-raw", paste0(nm,".nex")),
              file.path("dev/benchmarks", paste0(nm,".nex"))))
    if (file.exists(p)) return(fitchPhy(ReadAsPhyDat(p)))
  stop("not found: ", nm)
}
timed <- function(expr){t0<-Sys.time(); v<-force(expr); list(v=v, w=as.double(difftime(Sys.time(),t0,units="secs")))}

rows <- list()
for (nm in dnames) {
  phy <- loadPhy(nm); nt <- length(phy); at <- attributes(phy)
  ds <- list(contrast=at$contrast, tip_data=matrix(unlist(phy,use.names=FALSE),nrow=nt,byrow=TRUE),
             weight=at$weight, levels=at$levels)
  for (s in seeds) {
    set.seed(s); st <- RandomTree(nt, root=TRUE); st$tip.label <- names(phy)
    L <- TreeSearch:::ts_tbr_search(st$edge, ds$contrast, ds$tip_data, ds$weight, ds$levels, maxHits=1L)$edge
    for (k in cyc) {
      arms <- list(
        ratchet     = function() TreeSearch:::ts_ratchet_search(L, ds$contrast, ds$tip_data, ds$weight, ds$levels, nCycles=k),
        drift_exact = function(){Sys.setenv(TS_DRIFT_EXACT="1"); on.exit(Sys.unsetenv("TS_DRIFT_EXACT"))
                                 TreeSearch:::ts_drift_search(L, ds$contrast, ds$tip_data, ds$weight, ds$levels, nCycles=k)},
        drift_union = function() TreeSearch:::ts_drift_search(L, ds$contrast, ds$tip_data, ds$weight, ds$levels, nCycles=k))
      for (arm in names(arms)) {
        set.seed(1000L + s)
        r <- timed(arms[[arm]]())
        rows[[length(rows)+1L]] <- data.frame(dataset=nm, nTip=nt, seed=s, cycles=k,
                                               arm=arm, wall_s=round(r$w,3), score=as.double(r$v$score))
        cat(sprintf("%-12s seed=%d k=%2d %-12s score=%.0f wall=%.2fs\n", nm, s, k, arm, r$v$score, r$w))
      }
    }
  }
}
df <- do.call(rbind, rows)
write.csv(df, outfile, row.names=FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(df), outfile))
