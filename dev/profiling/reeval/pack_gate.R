# pack_gate.R -- TS_PACK_LOCAL score-exactness gate + per-candidate ns A/B.
# (1) EXACTNESS: TreeLength of K fixed random trees must be IDENTICAL under
#     pack OFF vs ON (per-block-local alphabet is a bijective relabel + "?"->
#     char-alphabet, Fitch-exact). Fixed trees => no trajectory confound.
# (2) SPEED: per-candidate ns (TS_IW_TIMING) over a near-optimal round, pack
#     OFF vs ON, min-of-runs.
# Usage: TS_LIB=.agent-pack Rscript dev/profiling/reeval/pack_gate.R <matrix.nex> [K=40] [reps=12]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-pack"), winslash="/"))
  library(TreeTools)
})
a <- commandArgs(trailingOnly=TRUE); nex<-a[[1]]
K    <- if (length(a)>=2) as.integer(a[[2]]) else 40L
reps <- if (length(a)>=3) as.integer(a[[3]]) else 12L

phy <- ReadAsPhyDat(nex); m<-PhyDatToMatrix(phy, ambigNA=FALSE); m[m=="-"]<-"?"; phy<-MatrixToPhyDat(m)
labels<-names(phy); nTip<-length(phy); nChar<-sum(attr(phy,"weight"))
contrast<-attr(phy,"contrast"); levels<-attr(phy,"levels")
tip_data<-matrix(unlist(phy,use.names=FALSE),nrow=nTip,byrow=TRUE); weight<-attr(phy,"weight")
cat(sprintf("=== %s : nTip=%d nChar=%d nStates=%d ===\n", basename(nex),nTip,nChar,length(levels)))

# (1) exactness via TreeLength on fixed trees
Sys.setenv(TS_PACK_LOCAL="0")
set.seed(11); trees <- lapply(seq_len(K), function(i){ set.seed(100+i); RandomTree(phy, root=TRUE) })
lens_off <- vapply(trees, function(t) TreeLength(t, phy), double(1))
Sys.setenv(TS_PACK_LOCAL="1")
lens_on  <- vapply(trees, function(t) TreeLength(t, phy), double(1))
Sys.setenv(TS_PACK_LOCAL="0")
nmis <- sum(abs(lens_off - lens_on) > 1e-9)
cat(sprintf("EXACTNESS: %d/%d trees match (max abs diff %.3g)  [%s]\n",
            K-nmis, K, max(abs(lens_off-lens_on)), if(nmis==0) "PASS" else "FAIL"))
if (nmis>0) { bad<-which(abs(lens_off-lens_on)>1e-9)[1]; cat(sprintf("  first mismatch: off=%.1f on=%.1f\n", lens_off[bad], lens_on[bad])) }

# (2) per-candidate ns A/B (near-optimal round)
runDiag <- function(edge) TreeSearch:::ts_tbr_diagnostics(edge,contrast,tip_data,weight,levels,maxHits=1L,acceptEqual=FALSE,maxChanges=0L)
parseIWT<-function(msg){ line<-grep("^IWT ",msg,value=TRUE); if(!length(line))return(NULL); line<-line[[length(line)]]
  g<-function(p){mm<-regmatches(line,regexec(p,line))[[1]]; if(length(mm)>=2) as.numeric(mm[[2]]) else NA_real_}
  list(spr=g("SPR n=[0-9.]+ ([0-9.]+)ns"),rer=g("REROOT n=[0-9.]+ ([0-9.]+)ns"),nsp=g("SPR n=([0-9.]+)"),nre=g("REROOT n=([0-9.]+)")) }
timeit <- function(pack){
  if(pack) Sys.setenv(TS_PACK_LOCAL="1") else Sys.setenv(TS_PACK_LOCAL="0")
  set.seed(1); start<-RandomTree(phy,root=TRUE); e0<-Preorder(RenumberTips(start,labels))[["edge"]]
  Sys.unsetenv("TS_IW_TIMING"); conv<-runDiag(e0); ec<-conv$edge
  Sys.setenv(TS_IW_TIMING="1"); spr<-rer<-nsp<-nre<-rep(NA_real_,reps)
  for(r in seq_len(reps)){ msg<-capture.output(.t<-runDiag(ec),type="message"); p<-parseIWT(msg)
    if(!is.null(p)){spr[r]<-p$spr;rer[r]<-p$rer;nsp[r]<-p$nsp;nre[r]<-p$nre} }
  Sys.unsetenv("TS_IW_TIMING")
  sprM<-suppressWarnings(min(spr,na.rm=TRUE));rerM<-suppressWarnings(min(rer,na.rm=TRUE))
  nspM<-stats::median(nsp,na.rm=TRUE);nreM<-stats::median(nre,na.rm=TRUE)
  tot<-(sprM*nspM+rerM*nreM)/(nspM+nreM); c(spr=sprM,rer=rerM,tot=tot)
}
off<-timeit(FALSE); on<-timeit(TRUE); Sys.setenv(TS_PACK_LOCAL="0")
cat(sprintf("SPEED  pack OFF: SPR %.2f REROOT %.2f TOT %.2f ns/cand\n", off["spr"],off["rer"],off["tot"]))
cat(sprintf("SPEED  pack ON : SPR %.2f REROOT %.2f TOT %.2f ns/cand\n", on["spr"],on["rer"],on["tot"]))
cat(sprintf("SPEEDUP (OFF/ON): SPR %.2fx REROOT %.2fx TOT %.2fx\n",
            off["spr"]/on["spr"], off["rer"]/on["rer"], off["tot"]/on["tot"]))
