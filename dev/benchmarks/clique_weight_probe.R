suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(igraph)}))
nex<-"/nobackup/pjjg18/lgsweep/matrices/project5432.nex"; FLOOR<-"/nobackup/pjjg18/floors/project5432_tnt_floor_one.tre"
pd<-ReadAsPhyDat(nex); tips<-names(pd); nt<-length(tips)
ch<-ReadCharacters(nex); rownames(ch)<-gsub(" ","_",rownames(ch)); ch<-ch[tips,,drop=FALSE]
miss<-c("?","-","N","n"); cl1<-function(x){x<-as.character(x); ifelse(!(x%in%miss)&nchar(x)==1L&grepl("^[0-9A-Za-z]$",x),x,NA_character_)}
M<-matrix(cl1(ch),nrow=nt)
codeL<-vector("list",ncol(M)); infm<-logical(ncol(M))
for(c in seq_len(ncol(M))){f<-factor(M[,c]); v<-as.integer(f); codeL[[c]]<-v; tb<-tabulate(v[!is.na(v)],nlevels(f)); infm[c]<-sum(tb>=2)>=2}
inf<-which(infm); ni<-length(inf)
compat<-function(a,b,mo=20){ok<-!is.na(a)&!is.na(b); if(sum(ok)<mo) return(FALSE); a<-a[ok];b<-b[ok]
  ed<-unique(cbind(paste0("A",a),paste0("B",b))); nodes<-unique(c(ed[,1],ed[,2])); idx<-setNames(seq_along(nodes),nodes); parent<-seq_along(nodes)
  fr<-function(i){while(parent[i]!=i)i<-parent[i];i}; cyc<-FALSE
  for(r in seq_len(nrow(ed))){u<-fr(idx[[ed[r,1]]]);v<-fr(idx[[ed[r,2]]]);if(u==v){cyc<-TRUE;break};parent[u]<-v}; !cyc}
deg<-integer(ni)
for(ii in seq_len(ni)){ai<-codeL[[inf[ii]]]; d<-0L; for(jj in seq_len(ni)) if(jj!=ii && compat(ai,codeL[[inf[jj]]])) d<-d+1L; deg[ii]<-d}
# effective sample size if sampling prob ~ degree  vs uniform (=ni)
w<-deg/sum(deg); ess<-1/sum(w^2)
cat(sprintf("informative chars=%d | degree: min=%d med=%.0f mean=%.1f max=%d CV=%.2f\n", ni, min(deg), median(deg), mean(deg), max(deg), sd(deg)/mean(deg)))
cat(sprintf("sampling prob ~ degree -> effective #chars = %.0f (uniform would be %d) -> ratio %.2f\n", ess, ni, ess/ni))
# 1943 deep clades (min-side>=100) + clean-synapomorphy support flag per char
tnt<-ReadTntTree(FLOOR,tipLabels=tips); if(inherits(tnt,"multiPhylo"))tnt<-tnt[[1]]; tnt<-Preorder(tnt)
sp<-as.Splits(tnt,tipLabels=tips); tis<-TipsInSplits(sp); ms<-pmin(tis,nt-tis); smat<-as.matrix(sp); cn<-colnames(smat); if(is.null(cn))cn<-tips
deepN<-as.integer(names(tis))[ms>=100]
clades<-lapply(deepN,function(nd){r<-as.logical(smat[as.character(nd),]); side<-cn[r]; if(length(side)>nt/2) setdiff(tips,side) else side})
supp<-logical(ni)
for(ii in seq_len(ni)){ v<-codeL[[inf[ii]]]; ok<-FALSE
  for(cl in clades){ A<-tips%in%cl; aT<-v[A]; bT<-v[!A]; aT<-aT[!is.na(aT)]; bT<-bT[!is.na(bT)]
    for(s in unique(aT)) if(sum(aT==s)>=2 && sum(bT==s)==0){ok<-TRUE;break}; if(ok)break }
  supp[ii]<-ok }
cat(sprintf("\n1943-deep-supporting chars: %d/%d\n", sum(supp), ni))
cat(sprintf("mean DEGREE: deep-supporting=%.1f  non-supporting=%.1f  (higher degree = more clique-central)\n", mean(deg[supp]), mean(deg[!supp])))
cat(sprintf("corr(degree, is-deep-supporting) = %.3f\n", suppressWarnings(cor(deg, as.numeric(supp)))))
# per-char steps on 1943 tree (homoplasy) if available
clen<-tryCatch(CharacterLength(tnt, pd, compress=FALSE), error=function(e) NULL)
if(!is.null(clen) && length(clen)==ncol(M)) cat(sprintf("corr(degree, steps@1943) = %.3f (n=%d)\n", suppressWarnings(cor(deg, clen[inf])), ni)) else cat("CharacterLength unavailable/misaligned\n")
