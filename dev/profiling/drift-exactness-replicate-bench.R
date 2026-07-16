#!/usr/bin/env Rscript
# Route-3 gap-closer: the +replicate arm. Marginal-value asks "is a wall-second
# better spent on exact-drift than on the incumbent lever?" -- and this project's
# real incumbent is MORE REPLICATES (policy-in-replicates-not-seconds:
# maxReplicates anytime-dominates), not ratchet. From the SAME per-seed local
# optimum L (identical construction to drift-exactness-route3-bench.R), spend the
# budget on k additional independent RAS+TBR restarts, tracking best-of(L, r1..rk)
# vs cumulative wall. Merge with route3_*.csv and compare frontiers: if
# drift_exact still beats `replicate` at matched wall -> drift earns recipe space.
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})
seeds  <- seq_len(as.integer(Sys.getenv("GATE_SEEDS", "20")))
kmax   <- as.integer(Sys.getenv("GATE_KMAX", "16"))   # max extra restarts
dnames <- strsplit(Sys.getenv("GATE_DATA", "Zhu2013,Zanol2014,Dikow2009"), ",")[[1]]
outfile <- Sys.getenv("GATE_OUT", "dev/profiling/drift-exactness-replicate-bench.csv")

fitchPhy <- function(p){m<-PhyDatToMatrix(p,ambigNA=FALSE); m[m=="-"]<-"?"; MatrixToPhyDat(m)}
loadPhy <- function(nm){
  for (p in c(file.path("data-raw",paste0(nm,".nex")), file.path("dev/benchmarks",paste0(nm,".nex"))))
    if (file.exists(p)) return(fitchPhy(ReadAsPhyDat(p)))
  stop("not found: ", nm)
}
tbrTo <- function(edge, ds) TreeSearch:::ts_tbr_search(edge, ds$contrast, ds$tip_data, ds$weight, ds$levels, maxHits=1L)

rows <- list()
for (nm in dnames) {
  phy <- loadPhy(nm); nt <- length(phy); at <- attributes(phy)
  ds <- list(contrast=at$contrast, tip_data=matrix(unlist(phy,use.names=FALSE),nrow=nt,byrow=TRUE),
             weight=at$weight, levels=at$levels)
  for (s in seeds) {
    # L: identical to the route-3 harness (same seed => same shared start).
    set.seed(s); st <- RandomTree(nt, root=TRUE); st$tip.label <- names(phy)
    t0 <- Sys.time(); Lr <- tbrTo(st$edge, ds); L_wall <- as.double(difftime(Sys.time(),t0,units="secs"))
    best <- Lr$score; cum <- 0.0
    # Continuation budget starts AFTER L (matching route-3, whose arms also
    # continue from L); cumulative wall is the continuation cost only.
    set.seed(1000L + s)
    kmark <- c(1,2,4,8,16); kmark <- kmark[kmark <= kmax]
    for (k in seq_len(kmax)) {
      r <- RandomTree(nt, root=TRUE); r$tip.label <- names(phy)
      tt <- Sys.time(); rr <- tbrTo(r$edge, ds); cum <- cum + as.double(difftime(Sys.time(),tt,units="secs"))
      if (rr$score < best) best <- rr$score
      if (k %in% kmark) {
        rows[[length(rows)+1L]] <- data.frame(dataset=nm, nTip=nt, seed=s, cycles=k,
                                              arm="replicate", wall_s=round(cum,3), score=as.double(best))
        cat(sprintf("%-12s seed=%d k=%2d replicate best=%.0f wall=%.2fs\n", nm, s, k, best, cum))
      }
    }
  }
}
df <- do.call(rbind, rows); write.csv(df, outfile, row.names=FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(df), outfile))
