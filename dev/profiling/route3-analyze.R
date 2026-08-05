dir <- "dev/profiling/route3-data"
sAtW <- function(d,W){s<-unique(d$seed); v<-sapply(s,function(x){r<-d[d$seed==x & d$wall_s<=W,]; if(nrow(r)==0) NA else min(r$score)}); mean(v,na.rm=TRUE)}

cat("################  RECIPE B (training MorphoBank, whole-search thorough)  ################\n")
cat("arms: none (driftCycles=0) / union / exact ; knob=maxReplicates{2,4,8}; matched wall\n")
rk <- Sys.glob(file.path(dir,"recipe_project*.csv"))
wintab <- c(none=0,union=0,exact=0,tie=0)
for(f in rk){ d<-read.csv(f); key<-sub("_.*","",sub("recipe_","",basename(f)))
  arms<-c("none","union","exact"); mw<-max(d$wall_s); Ws<-signif(mw*c(.15,.35,.6,1),2)
  # winner at the LARGEST matched wall (most budget) per key
  vE<-sAtW(d[d$arm=="exact",],mw); vU<-sAtW(d[d$arm=="union",],mw); vN<-sAtW(d[d$arm=="none",],mw)
  best<-which.min(c(none=vN,union=vU,exact=vE)); w<-names(best)
  if(max(abs(c(vN,vU,vE)-min(c(vN,vU,vE))))<1e-9) w<-"tie"
  wintab[w]<-wintab[w]+1
  cat(sprintf("\n%-14s nTip=%d  @maxwall %.1fs: none=%.1f union=%.1f exact=%.1f  -> %s\n",
              key, d$nTip[1], mw, vN,vU,vE, toupper(w)))
  # per-arm mean score at reps=8 and mean wall
  for(a in arms){da<-d[d$arm==a & d$reps==max(d$reps),]; cat(sprintf("    %-5s reps=8 mean score=%.1f wall=%.1fs\n",a,mean(da$score),mean(da$wall_s)))}
}
cat(sprintf("\n== RECIPE win-tally @max budget across %d keys: none=%d union=%d exact=%d tie=%d ==\n",
            length(rk), wintab["none"],wintab["union"],wintab["exact"],wintab["tie"]))

cat("\n\n################  ROUTE-3 + REPLICATE (marginal value from shared local-opt)  ################\n")
cat("arms: ratchet / drift_exact / drift_union / replicate ; knob=cycles or restarts; matched wall\n")
for(nm in c("Zhu2013","Zanol2014","Dikow2009")){
  r3<-file.path(dir,sprintf("route3_%s.csv",nm)); rp<-file.path(dir,sprintf("repl_%s.csv",nm))
  if(!file.exists(r3)||!file.exists(rp)) next
  d<-rbind(read.csv(r3), read.csv(rp)); arms<-c("ratchet","drift_exact","drift_union","replicate")
  mw<-min(sapply(arms,function(a) max(d$wall_s[d$arm==a]))); Ws<-signif(mw*c(.2,.5,1),2)
  cat(sprintf("\n== %s : matched-WALL frontier (mean score) ==\n",nm))
  cat(sprintf("  %6s %9s %11s %11s %10s   best\n","wall","ratchet","drift_exact","drift_union","replicate"))
  for(W in Ws){ v<-sapply(arms,function(a) sAtW(d[d$arm==a,],W))
    cat(sprintf("  %6.2f %9.2f %11.2f %11.2f %10.2f   %s\n",W,v["ratchet"],v["drift_exact"],v["drift_union"],v["replicate"],names(which.min(v)))) }
}
