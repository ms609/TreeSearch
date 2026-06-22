# Comprehensive byte-identity validation for the IW optimizations.
# Compares the pure-scalar baseline (both kill-switches on) against:
#   dirty-only (x4 off), x4-only (dirty off), both-on (default).
# Each must give identical score AND n_evaluated. SCANCHK on => per-accept
# predicted==actual (catches a wrong base_iw from the dirty-region).
# Pure-IW (recode -> ?). Same R seed before EACH call (make_rng <- unif_rand).
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  CONC <- if (length(args) >= 2) as.numeric(args[[2]]) else 10
  library(TreeSearch, lib.loc = LIB); library(TreeTools)
}))
data("inapplicable.phyData", package = "TreeSearch")
Sys.setenv(TS_IW_SCANCHK = "1")
recode <- function(phy){at<-attributes(phy);ct<-at$contrast;lv<-at$levels;al<-at$allLevels
 ic<-which(lv=="-");dr<-which(al=="-");if(length(dr))ct[dr,]<-1
 if(length(ic)){ct<-ct[,-ic,drop=FALSE];lv<-lv[-ic]};attr(phy,"contrast")<-ct;attr(phy,"levels")<-lv;phy}
prep <- function(phy){at<-attributes(phy)
 list(contrast=at$contrast,tip=matrix(unlist(phy,use.names=FALSE),nrow=length(phy),byrow=TRUE),
      weight=at$weight,levels=at$levels)}
setflags <- function(nox4, nodirty){
  if(nox4) Sys.setenv(TS_IW_NOX4="1") else Sys.unsetenv("TS_IW_NOX4")
  if(nodirty) Sys.setenv(TS_IW_NODIRTY="1") else Sys.unsetenv("TS_IW_NODIRTY")
}
run <- function(p, edge){ set.seed(424242L)
  TreeSearch:::ts_tbr_search(edge,p$contrast,p$tip,p$weight,p$levels,maxHits=1L,concavity=CONC)}

DATASETS <- c("Vinther2008","Wortley2006","Zanol2014","Giles2015","Dikow2009")
configs <- list(base=c(T,T), dirty=c(T,F), x4=c(F,T), both=c(F,F))
for (nm in DATASETS) {
  phy<-recode(inapplicable.phyData[[nm]]); p<-prep(phy)
  set.seed(7294); edge<-RandomTree(names(inapplicable.phyData[[nm]]),root=TRUE)$edge
  res<-list()
  for(cfg in names(configs)){ setflags(configs[[cfg]][1],configs[[cfg]][2]); res[[cfg]]<-run(p,edge) }
  b<-res$base
  ok<-TRUE; for(cfg in c("dirty","x4","both")){
    if(!isTRUE(all.equal(b$score,res[[cfg]]$score)) || b$n_evaluated!=res[[cfg]]$n_evaluated) ok<-FALSE }
  message(sprintf("%-12s base score=%.4f neval=%d | dirty %.4f/%d | x4 %.4f/%d | both %.4f/%d => %s",
    nm, b$score, b$n_evaluated,
    res$dirty$score,res$dirty$n_evaluated, res$x4$score,res$x4$n_evaluated,
    res$both$score,res$both$n_evaluated, if(ok)"IDENTICAL" else "*** DIFF ***"))
}
