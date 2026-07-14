# Combined frontier: original nc{8..256} + extension nc{512,1024}, matched wall.
dir <- "dev/profiling"
read1 <- function(f) read.csv(file.path(dir, f))
scoreAtWall <- function(d, W) {
  sd <- unique(d$seed)
  v <- sapply(sd, function(s){x<-d[d$seed==s & d$wall_s<=W,]; if(nrow(x)==0) NA else min(x$score)})
  mean(v, na.rm=TRUE)
}
for (nm in c("Zhu2013","Zanol2014")) {
  df <- rbind(read1(sprintf("drift-exactness-gate-hamilton-%s.csv", nm)),
              read1(sprintf("drift-exactness-gate-hamilton-%s-ext.csv", nm)))
  cat(sprintf("\n== %s : mean score / mean wall by nCycles x scorer ==\n", nm))
  ag <- aggregate(cbind(score,wall_s)~nCycles+scorer, df, mean)
  ag <- ag[order(ag$scorer, ag$nCycles),]; ag$score<-round(ag$score,2); ag$wall_s<-round(ag$wall_s,2)
  print(ag, row.names=FALSE)
  cat(sprintf("\n-- %s matched-WALL frontier (mean over 30 seeds) --\n", nm))
  cat(sprintf("  %6s %9s %9s   %s\n","wall","union","exact","winner(margin)"))
  for (W in c(0.8,1.6,3.2,6.4,9,12,16,20)) {
    u<-scoreAtWall(df[df$scorer=="union",],W); e<-scoreAtWall(df[df$scorer=="exact",],W)
    win <- if(any(is.na(c(u,e)))) "-" else if(abs(u-e)<1e-9) "tie" else if(e<u) sprintf("EXACT (%.2f)",u-e) else sprintf("union (%.2f)",e-u)
    cat(sprintf("  %6.1f %9.2f %9.2f   %s\n",W,u,e,win))
  }
  du<-df[df$scorer=="union",]; de<-df[df$scorer=="exact",]
  cat(sprintf("  best-of-30 any budget: union=%.0f exact=%.0f ; mean@maxNc(%d): union=%.2f exact=%.2f\n",
      min(du$score),min(de$score),max(df$nCycles),
      mean(du$score[du$nCycles==max(du$nCycles)]),mean(de$score[de$nCycles==max(de$nCycles)])))
}
