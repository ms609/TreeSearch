# Clean multi-rep timing of IW reroot x4 vs scalar. Deterministic trajectory
# (fixed seed each call) => identical candidate sequence; only timing varies, so
# the median over reps is a clean per-candidate measurement. Emits TAG lines so
# the (dataset,config) of each IWT line is recoverable downstream.
# Usage: Rscript bench_iw_x4_timing.R [lib] [reps]
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  REPS <- if (length(args) >= 2) as.integer(args[[2]]) else 7L
  library(TreeSearch, lib.loc = LIB); library(TreeTools)
}))
data("inapplicable.phyData", package = "TreeSearch")
Sys.setenv(TS_IW_TIMING = "1")
recode <- function(phy){at<-attributes(phy);ct<-at$contrast;lv<-at$levels;al<-at$allLevels
 ic<-which(lv=="-");dr<-which(al=="-");if(length(dr))ct[dr,]<-1
 if(length(ic)){ct<-ct[,-ic,drop=FALSE];lv<-lv[-ic]};attr(phy,"contrast")<-ct;attr(phy,"levels")<-lv;phy}
DATASETS <- c("Wortley2006","Zanol2014","Giles2015","Dikow2009")
for (nm in DATASETS) {
  phy<-recode(inapplicable.phyData[[nm]]);at<-attributes(phy)
  tip<-matrix(unlist(phy,use.names=FALSE),nrow=length(phy),byrow=TRUE)
  set.seed(7294); edge<-RandomTree(names(inapplicable.phyData[[nm]]),root=TRUE)$edge
  for (cfg in c("scalar","x4")) {
    if (cfg=="scalar") Sys.setenv(TS_IW_NOX4="1") else Sys.unsetenv("TS_IW_NOX4")
    for (r in seq_len(REPS)) {
      message(sprintf("TAG %s %s", nm, cfg))
      set.seed(424242L)
      invisible(TreeSearch:::ts_tbr_search(edge, at$contrast, tip, at$weight,
                  at$levels, maxHits = 1L, concavity = 10))
    }
  }
}
