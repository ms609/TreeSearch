suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(ape)}))
task <- as.integer(commandArgs(trailingOnly=TRUE)[[1]])
NEX<-"/nobackup/pjjg18/lgsweep/matrices/project5432.nex"; OUT<-"/nobackup/pjjg18/reeval/clique_out"; dir.create(OUT,showWarnings=FALSE)
pd<-ReadAsPhyDat(NEX); tips<-names(pd)
TL<-function(t) TreeLength(Preorder(t),pd,concavity=Inf,inapplicable="missing")
D<-readRDS("/nobackup/pjjg18/reeval/clique_trees.rds"); trees<-D$trees; info<-D$info; NC<-length(trees); SEEDS<-5
if(task<=NC*SEEDS){ ci<-((task-1)%%NC)+1; seed<-((task-1)%/%NC)+1; ctree<-trees[[ci]]; ov<-info$overlap1943[ci]; lab<-paste0("cq",ci)
} else { r<-task-NC*SEEDS; ci<-0; seed<-r; ctree<-NULL; ov<-NA_real_; lab<-"RAS" }
set.seed(3000+task); t0<-Sys.time()
res<-tryCatch(MaximizeParsimony(pd,constraint=ctree,concavity=Inf,inapplicable="missing",strategy="thorough",
  maxReplicates=8,targetHits=NULL,maxSeconds=1800,nThreads=1L,verbosity=0L,collapse=FALSE),
  error=function(e){cat("ERR:",conditionMessage(e),"\n");NULL})
el<-as.numeric(difftime(Sys.time(),t0,units="secs")); if(is.null(res)){cat("NULL",task,"\n");quit(save="no")}
trs<-if(inherits(res,"phylo"))list(res) else res; sc<-vapply(trs,TL,numeric(1)); best<-min(sc); bt<-Preorder(trs[[which.min(sc)]])
held<-NA_real_; if(!is.null(ctree)){sf<-tryCatch(SplitFrequency(ctree,structure(list(bt),class="multiPhylo")),error=function(e)NULL); if(!is.null(sf)) held<-mean(!is.na(sf)&sf>=1)}
cat(sprintf("task=%d %s seed=%d best=%.0f hit1944=%d hit1943=%d ov=%.2f held=%.2f el=%.0fs\n",task,lab,seed,best,best<=1944,best<=1943,ov,held,el))
write.csv(data.frame(task=task,clique=lab,cidx=ci,seed=seed,best=best,hit1944=best<=1944,hit1943=best<=1943,overlap1943=ov,held=held,elapsed=el),
  file.path(OUT,sprintf("cell_%03d.csv",task)),row.names=FALSE)
