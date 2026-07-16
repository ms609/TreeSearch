#!/usr/bin/env Rscript
# Mission A angle #1 test (a), DEFINITIVE: does the concurrent session's EXACT verified
# clean-synapomorphy count DISCRIMINATE true deep splits from wrong deep splits?
# Their def (h2_and_synap.R): A = SMALLER side (the clade); a char is a clean synapomorphy
# if some single state appears in >=2 tips of A and 0 scored tips of B. "near" = <=1 in B.
# Applied identically to TRUE (1943-absent-from-TS) and WRONG (TS, absent-from-1943) splits.
# + a missing-data-controlled clean (require B adequately scored) per advisor.
# AUC=P(count(true)>count(wrong)); emphasise size>=60 (n=21).
suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(ape)}))
a <- commandArgs(trailingOnly = TRUE); nex <- a[[1]]; csv <- a[[2]]
d <- ReadAsPhyDat(nex); tips <- names(d); nt <- length(tips)
ch <- ReadCharacters(nex); rownames(ch) <- gsub(" ", "_", rownames(ch))
stopifnot(all(tips %in% rownames(ch))); ch <- ch[tips, , drop = FALSE]
miss <- c("?","-","N","n")
clean1 <- function(x){ x <- as.character(x); ifelse(!(x %in% miss) & nchar(x)==1L & grepl("^[0-9A-Za-z]$",x), x, NA_character_) }
M <- matrix(clean1(ch), nrow = nt)
code <- matrix(NA_integer_, nt, ncol(M)); K <- integer(ncol(M))
for (c in seq_len(ncol(M))) { f <- factor(M[,c]); code[,c] <- as.integer(f); K[c] <- nlevels(f) }
cat(sprintf("matrix %d tips x %d chars\n", nt, ncol(M)))

S <- read.csv(csv, stringsAsFactors = FALSE)
S <- S[S$class %in% c("true_missing","wrong") & S$size >= 30, ]
keyIdx <- strsplit(S$key, ",")

## their exact clean/near on the SMALLER side; + controlled clean
counts <- function(ix){
  ix <- as.integer(ix)
  A <- if (length(ix) <= nt/2) ix else setdiff(seq_len(nt), ix)   # smaller side = clade
  Bl <- setdiff(seq_len(nt), A); nB <- length(Bl)
  cl <- 0L; nr <- 0L; clc <- 0L
  thrB <- max(4, 0.5*nB)
  for (c in seq_len(ncol(M))) {
    k <- K[c]; if (k==0L) next
    ca <- code[A,c]; cb <- code[Bl,c]
    ca <- ca[!is.na(ca)]; cb <- cb[!is.na(cb)]; sb <- length(cb)
    if (length(ca) < 2L) next
    tA <- tabulate(ca, k); tB <- tabulate(cb, k)
    if (any(tA>=2L & tB==0L)) { cl <- cl+1L; if (sb >= thrB) clc <- clc+1L }
    if (any(tA>=2L & tB<=1L)) nr <- nr+1L
  }
  c(clean=cl, near=nr, clean_ctrl=clc)
}
CM <- t(vapply(keyIdx, counts, integer(3)))
S$clean <- CM[,"clean"]; S$near <- CM[,"near"]; S$clean_ctrl <- CM[,"clean_ctrl"]

auc <- function(pos,neg){ if(!length(pos)||!length(neg)) return(NA_real_)
  r<-rank(c(pos,neg)); (sum(r[seq_along(pos)])-length(pos)*(length(pos)+1)/2)/(length(pos)*length(neg)) }
band <- function(sig,lab){ cat(sprintf("\n=== %s  AUC=P(true>wrong) ===\n", lab))
  for (thr in c(30,60,120)){ tm<-S[S$class=="true_missing"&S$size>=thr,sig]; wr<-S[S$class=="wrong"&S$size>=thr,sig]
    cat(sprintf("  size>=%3d: true n=%2d med=%.1f | wrong n=%3d med=%.1f | AUC=%.3f\n",
      thr,length(tm),median(tm),length(wr),median(wr),auc(tm,wr))) } }
band("clean","CLEAN (their exact: >=2 in clade, 0 in scored B)")
band("near","NEAR (>=2 in clade, <=1 in B)")
band("clean_ctrl","CLEAN + missing-control (B scored >=50%)")

tmD<-S[S$class=="true_missing",]; dp<-tmD[which.max(tmD$size),]; comp<-S[S$size>=0.7*dp$size,]
cat(sprintf("\n238-split: size=%d clean=%d near=%d clean_ctrl=%d | clean pctile among %d deep = %.1f%%\n",
    dp$size,dp$clean,dp$near,dp$clean_ctrl,nrow(comp),100*mean(comp$clean<=dp$clean)))
cat(sprintf("sanity: true deep(min>=100) clean median=%.1f max=%d (concurrent node reported 4-10)\n",
    median(S$clean[S$class=="true_missing"&S$size>=100]), max(S$clean[S$class=="true_missing"&S$size>=100])))
write.csv(S, "/nobackup/pjjg18/reeval/detect5432_synap2.csv", row.names=FALSE)
cat("\nVERDICT: AUC>>0.5 at size>=60 -> discriminates -> synapomorphy-clique LIVE.",
    "\n         AUC~0.5 both -> per-split dead in proportion AND absolute -> angle #1 CLOSED.\n")
