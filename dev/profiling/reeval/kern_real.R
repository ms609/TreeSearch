# kern_real.R -- per-candidate kernel timing on REAL large matrices, near-optimal
# regime (the ~47ns target). Load a .nex, fitch-convert (- -> ?), converge to a
# kernel local optimum, then time a near-optimal round of TBR via TS_IW_TIMING;
# min-of-runs. This is the A/B harness for kernel re-slicing (BASE vs a flagged
# variant): run twice with the toggle env set, compare SPR/REROOT/tot ns.
#
# Usage: TS_LIB=.agent-disc Rscript dev/profiling/reeval/kern_real.R <matrix.nex> [reps=15]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-disc"), winslash="/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly=TRUE)
nex  <- args[[1]]
reps <- if (length(args) >= 2) as.integer(args[[2]]) else 15L

phy <- ReadAsPhyDat(nex)
m <- PhyDatToMatrix(phy, ambigNA=FALSE); m[m=="-"] <- "?"; phy <- MatrixToPhyDat(m)
labels <- names(phy); nTip <- length(phy); nChar <- sum(attr(phy,"weight"))
contrast <- attr(phy,"contrast"); levels <- attr(phy,"levels")
tip_data <- matrix(unlist(phy, use.names=FALSE), nrow=nTip, byrow=TRUE)
weight <- attr(phy,"weight")
cat(sprintf("=== %s : nTip=%d nChar=%d nStates=%d nPat=%d ===\n",
            basename(nex), nTip, nChar, length(levels), length(weight)))

runDiag <- function(edge) TreeSearch:::ts_tbr_diagnostics(
  edge, contrast, tip_data, weight, levels, maxHits=1L, acceptEqual=FALSE, maxChanges=0L)

parseIWT <- function(msg){
  line <- grep("^IWT ", msg, value=TRUE); if(!length(line)) return(NULL)
  line <- line[[length(line)]]
  g <- function(p){ mm<-regmatches(line, regexec(p,line))[[1]]; if(length(mm)>=2) as.numeric(mm[[2]]) else NA_real_ }
  list(clips=g("clips=([0-9.]+)"), n_spr=g("SPR n=([0-9.]+)"), spr=g("SPR n=[0-9.]+ ([0-9.]+)ns"),
       n_rer=g("REROOT n=([0-9.]+)"), rer=g("REROOT n=[0-9.]+ ([0-9.]+)ns"))
}

set.seed(1); start <- RandomTree(phy, root=TRUE)
edge0 <- Preorder(RenumberTips(start, labels))[["edge"]]
Sys.unsetenv("TS_IW_TIMING")
conv <- runDiag(edge0); edgeC <- conv$edge

Sys.setenv(TS_IW_TIMING="1")
spr<-rer<-clp<-nsp<-nre<-rep(NA_real_, reps)
for (r in seq_len(reps)) {
  msg <- capture.output(.t <- runDiag(edgeC), type="message")
  p <- parseIWT(msg); if(!is.null(p)){spr[r]<-p$spr;rer[r]<-p$rer;clp[r]<-p$clips;nsp[r]<-p$n_spr;nre[r]<-p$n_rer}
}
Sys.unsetenv("TS_IW_TIMING")
sprM<-suppressWarnings(min(spr,na.rm=TRUE)); rerM<-suppressWarnings(min(rer,na.rm=TRUE))
nspM<-stats::median(nsp,na.rm=TRUE); nreM<-stats::median(nre,na.rm=TRUE)
totM<-(sprM*nspM+rerM*nreM)/(nspM+nreM)
cat(sprintf("clips=%.0f | SPR n=%.0f %.2fns | REROOT n=%.0f %.2fns | TOT/cand %.2fns\n",
            stats::median(clp,na.rm=TRUE), nspM, sprM, nreM, rerM, totM))
