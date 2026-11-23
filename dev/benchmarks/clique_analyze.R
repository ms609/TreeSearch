d<-do.call(rbind,lapply(Sys.glob("/nobackup/pjjg18/reeval/clique_out/cell_*.csv"),read.csv,stringsAsFactors=FALSE))
cat(sprintf("cells=%d\n\n",nrow(d)))
agg<-aggregate(cbind(hit1944,hit1943,best,held)~clique+overlap1943,d,function(x)round(mean(x),3))
ras<-d[d$clique=="RAS",]
cat("RAS control: hit1944=",mean(ras$hit1944)," med best=",median(ras$best),"\n\n")
o<-d[d$clique!="RAS",]; ag<-do.call(rbind,lapply(split(o,o$clique),function(x)data.frame(clique=x$clique[1],ov=x$overlap1943[1],n=nrow(x),hit1944=mean(x$hit1944),hit1943=mean(x$hit1943),medBest=median(x$best),minBest=min(x$best),held=round(mean(x$held),2))))
ag<-ag[order(-ag$hit1944,-ag$ov),]; print(ag,row.names=FALSE)
cat(sprintf("\nANY clique-start hit<=1944: %d/%d cells | best overall=%d\n",sum(o$hit1944),nrow(o),min(o$best)))
cat(sprintf("corr(overlap1943, hit1944) = %.2f\n", suppressWarnings(cor(o$overlap1943,o$hit1944))))
