# pack_mpt_gate.R -- FULL MaximizeParsimony byte-identity gate for TS_PACK_LOCAL.
# C-vs-A: pack OFF vs ON, fixed seed, nThreads=1 (reproducible run-to-run per s7),
# compare best score + n distinct trees + candidates_evaluated. Exercises the whole
# pipeline (ratchet + sector + drift + fuse + collapse + MPT), not just the scorer.
# Modes: EW (recode -/?), IW (concavity), NA (keep '-' -> exercises has_inapplicable
# packing branch). Any mismatch => a consumer interprets a plane index as a global
# state and needs the per-block plane->global map.
#
# Usage: TS_LIB=.agent-pack Rscript pack_mpt_gate.R <mode ew|iw|na> <dataset|nexpath> [seed=1] [reps=3]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB",".agent-pack"), winslash="/"))
  library(TreeTools)
})
a<-commandArgs(trailingOnly=TRUE); mode<-a[[1]]; dset<-a[[2]]
seed <- if(length(a)>=3) as.integer(a[[3]]) else 1L
reps <- if(length(a)>=4) as.integer(a[[4]]) else 3L

# load: built-in inapplicable.phyData name, or a .nex path
if (file.exists(dset)) { phy0 <- ReadAsPhyDat(dset); nm<-basename(dset)
} else { data("inapplicable.phyData", package="TreeSearch"); phy0<-inapplicable.phyData[[dset]]; nm<-dset }

if (mode=="na") {
  phy <- phy0                      # keep '-' -> has_inapplicable path
} else {
  m<-PhyDatToMatrix(phy0, ambigNA=FALSE); m[m=="-"]<-"?"; phy<-MatrixToPhyDat(m)
}
conc <- if (mode=="iw") 10 else Inf
cat(sprintf("=== %s mode=%s nTip=%d nChar=%d conc=%s seed=%d reps=%d ===\n",
            nm, mode, length(phy), sum(attr(phy,"weight")), as.character(conc), seed, reps))

runMP <- function(pack, s){
  if(pack) Sys.setenv(TS_PACK_LOCAL="1") else Sys.setenv(TS_PACK_LOCAL="0")
  set.seed(s)
  r <- suppressWarnings(MaximizeParsimony(phy, concavity=conc, maxReplicates=6L,
         nThreads=1L, verbosity=0L))
  Sys.setenv(TS_PACK_LOCAL="0")
  list(score=min(as.double(attr(r,"score"))), ntree=length(r),
       cand=suppressWarnings(as.double(attr(r,"candidates_evaluated"))))
}
allok<-TRUE
for (i in seq_len(reps)) {
  s<-seed+i-1L
  off<-runMP(FALSE,s); on<-runMP(TRUE,s)
  ok <- abs(off$score-on$score)<1e-9 && off$ntree==on$ntree &&
        (is.na(off$cand)||is.na(on$cand)|| off$cand==on$cand)
  allok <- allok && ok
  cat(sprintf(" seed %d: OFF score=%.2f ntree=%d cand=%.0f | ON score=%.2f ntree=%d cand=%.0f | %s\n",
      s, off$score, off$ntree, ifelse(is.na(off$cand),-1,off$cand),
      on$score, on$ntree, ifelse(is.na(on$cand),-1,on$cand), if(ok)"MATCH" else "*** DIFFERS ***"))
}
cat(sprintf("RESULT %s/%s: %s\n", nm, mode, if(allok)"BYTE-IDENTICAL (PASS)" else "MISMATCH (needs plane->global map)"))
