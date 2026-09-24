# tbr_iw_multistart.R -- multi-start IW search quality on one dataset.  Runs the
# kernel IW TBR (unrooted) from many random starts and reports the BEST IW score
# reached, so single-start basin shifts (expected from the exact-scoring fix)
# wash out.  Use to check the fix did not degrade achievable IW quality.
#
# Usage: Rscript dev/benchmarks/tbr_iw_multistart.R [nTip] [idx] [concavity] [nStart]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
a <- commandArgs(trailingOnly = TRUE)
nTip <- if (length(a) >= 1) as.integer(a[[1]]) else 16L
idx  <- if (length(a) >= 2) as.integer(a[[2]]) else 19L
conc <- if (length(a) >= 3) as.numeric(a[[3]]) else 10
nStart <- if (length(a) >= 4) as.integer(a[[4]]) else 40L
nChar <- 60L; nState <- 3L

set.seed(1000L + idx)
tips <- paste0("t", seq_len(nTip))
m <- matrix(sample(0:(nState-1L), nTip*nChar, replace=TRUE), nrow=nTip, dimnames=list(tips,NULL))
phy <- phangorn::phyDat(m, type="USER", levels=as.character(0:(nState-1L)))
at <- attributes(phy)
d <- list(contrast=at$contrast,
          tip_data=matrix(unlist(phy,use.names=FALSE), nrow=length(phy), byrow=TRUE),
          weight=at$weight, levels=at$levels, labels=names(phy), phy=phy)

kbest <- function(seedBase) {
  best <- Inf
  for (s in seq_len(nStart)) {
    set.seed(seedBase + s); st <- RandomTree(d$phy, root=TRUE)
    edge <- Preorder(RenumberTips(st, d$labels))[["edge"]]
    set.seed(s)
    res <- TreeSearch:::ts_tbr_diagnostics(edge, d$contrast, d$tip_data, d$weight, d$levels,
             maxHits=1L, acceptEqual=FALSE, maxChanges=0L, concavity=conc, unrooted=TRUE)
    if (res$score < best) best <- res$score
  }
  best
}
b <- kbest(50000L)
cat(sprintf("nTip=%d #%d conc=%g  %d random starts  =>  BEST IW = %.5f\n",
            nTip, idx, conc, nStart, b))
