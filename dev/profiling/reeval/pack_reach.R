# pack_reach.R -- BEHAVIORAL validation of TS_PACK_LOCAL (not byte-identity).
# The real bar (user 2026-07-15): fast to best score, good MPT spread, reach
# preserved. Per seed, OFF vs ON: best score, n MPT trees, wall. Aggregate:
# reach (min score over seeds), per-seed score deltas, MPT spread, total wall.
# Usage: TS_LIB=.agent-pack Rscript pack_reach.R <dataset> <mode ew|iw> [nSeed=8] [maxRep=]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB",".agent-pack"), winslash="/"))
  library(TreeTools)
})
a<-commandArgs(trailingOnly=TRUE); dset<-a[[1]]; mode<-if(length(a)>=2)a[[2]] else "ew"
nSeed<-if(length(a)>=3) as.integer(a[[3]]) else 8L
maxRep<-if(length(a)>=4) as.integer(a[[4]]) else NA_integer_
if (file.exists(dset)) { phy0<-ReadAsPhyDat(dset); nm<-basename(dset)
} else { data("inapplicable.phyData",package="TreeSearch"); phy0<-inapplicable.phyData[[dset]]; nm<-dset }
m<-PhyDatToMatrix(phy0,ambigNA=FALSE); m[m=="-"]<-"?"; phy<-MatrixToPhyDat(m)
conc<-if(mode=="iw") 10 else Inf
args<-list(dataset=phy, concavity=conc, nThreads=1L, verbosity=0L)
if(!is.na(maxRep)) args$maxReplicates<-maxRep
cat(sprintf("=== REACH %s mode=%s nTip=%d nChar=%d nSeed=%d ===\n", nm,mode,length(phy),sum(attr(phy,"weight")),nSeed))
run<-function(pk,s){ if(pk)Sys.setenv(TS_PACK_LOCAL="1") else Sys.setenv(TS_PACK_LOCAL="0"); set.seed(s)
  t0<-Sys.time(); r<-suppressWarnings(do.call(MaximizeParsimony,args)); w<-as.double(difftime(Sys.time(),t0,units="secs"))
  Sys.setenv(TS_PACK_LOCAL="0"); c(score=min(as.double(attr(r,"score"))), ntree=length(r), wall=w) }
O<-N<-matrix(NA_real_,nSeed,3,dimnames=list(NULL,c("score","ntree","wall")))
for(i in seq_len(nSeed)){ s<-i; O[i,]<-run(FALSE,s); N[i,]<-run(TRUE,s)
  cat(sprintf(" s%d OFF %.1f/%dt/%.1fs | ON %.1f/%dt/%.1fs\n",s,O[i,1],O[i,2],O[i,3],N[i,1],N[i,2],N[i,3])) }
cat(sprintf("\nREACH (min score): OFF %.1f  ON %.1f  %s\n", min(O[,1]),min(N[,1]),
  if(min(N[,1])<=min(O[,1])+1e-9)"ON reaches >= OFF" else "*** ON WORSE REACH ***"))
cat(sprintf("mean score: OFF %.2f ON %.2f (Δ %.2f)\n", mean(O[,1]),mean(N[,1]),mean(N[,1])-mean(O[,1])))
cat(sprintf("mean n MPT: OFF %.1f ON %.1f\n", mean(O[,2]),mean(N[,2])))
cat(sprintf("total wall: OFF %.1fs ON %.1fs  (ON/OFF %.2fx; <1 = ON faster)\n", sum(O[,3]),sum(N[,3]),sum(N[,3])/sum(O[,3])))
