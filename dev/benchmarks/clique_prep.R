suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(igraph)}))
nex <- "/nobackup/pjjg18/lgsweep/matrices/project5432.nex"
FLOOR <- "/nobackup/pjjg18/floors/project5432_tnt_floor_one.tre"
pd <- ReadAsPhyDat(nex); tips <- names(pd); nt <- length(tips)
ch <- ReadCharacters(nex); rownames(ch) <- gsub(" ","_",rownames(ch)); ch <- ch[tips,,drop=FALSE]
miss <- c("?","-","N","n"); cl1<-function(x){x<-as.character(x); ifelse(!(x%in%miss)&nchar(x)==1L&grepl("^[0-9A-Za-z]$",x),x,NA_character_)}
M <- matrix(cl1(ch), nrow=nt); rownames(M) <- tips
codeL<-vector("list",ncol(M)); infm<-logical(ncol(M))
for(c in seq_len(ncol(M))){ f<-factor(M[,c]); v<-as.integer(f); codeL[[c]]<-v; tb<-tabulate(v[!is.na(v)],nlevels(f)); infm[c]<-sum(tb>=2)>=2 }
inf<-which(infm); ni<-length(inf)
compat<-function(a,b,minOv=20){ ok<-!is.na(a)&!is.na(b); if(sum(ok)<minOv) return(FALSE); a<-a[ok];b<-b[ok]
  ed<-unique(cbind(paste0("A",a),paste0("B",b))); nodes<-unique(c(ed[,1],ed[,2])); idx<-setNames(seq_along(nodes),nodes); parent<-seq_along(nodes)
  findR<-function(i){while(parent[i]!=i)i<-parent[i];i}; cyc<-FALSE
  for(r in seq_len(nrow(ed))){u<-findR(idx[[ed[r,1]]]);v<-findR(idx[[ed[r,2]]]);if(u==v){cyc<-TRUE;break};parent[u]<-v}; !cyc }
el<-list();k<-0L
for(ii in seq_len(ni-1)){ai<-codeL[[inf[ii]]]; for(jj in (ii+1):ni) if(compat(ai,codeL[[inf[jj]]])){k<-k+1L;el[[k]]<-c(ii,jj)}}
E<-do.call(rbind,el); G<-make_empty_graph(n=ni,directed=FALSE); G<-add_edges(G,as.vector(t(E)))
mc <- max_cliques(G, min=30); sz<-vapply(mc,length,integer(1)); ord<-order(sz,decreasing=TRUE)
picked<-integer(0); psets<-list()
for(i in ord){ s<-as.integer(mc[[i]]); ok<-TRUE
  for(p in psets){ if(length(intersect(s,p))/length(union(s,p))>0.6){ok<-FALSE;break} }
  if(ok){picked<-c(picked,i); psets<-c(psets,list(s))}; if(length(picked)>=12) break }
tnt<-ReadTntTree(FLOOR,tipLabels=tips); if(inherits(tnt,"multiPhylo"))tnt<-tnt[[1]]; tnt<-Preorder(tnt)
key <- function(spl){ m<-as.matrix(spl); apply(m,1,function(r){v<-as.logical(r); if(v[1])v<-!v; paste(which(v),collapse=",")}) }
k1943 <- unique(key(as.Splits(tnt,tipLabels=tips)))
buildCT<-function(cc){ sub<-M[,inf[cc],drop=FALSE]; sub[is.na(sub)]<-"?"; pdS<-MatrixToPhyDat(sub)
  r<-MaximizeParsimony(pdS,concavity=Inf,inapplicable="missing",strategy="thorough",maxReplicates=4,collapse=TRUE,verbosity=0L)
  if(inherits(r,"phylo")) r<-structure(list(r),class="multiPhylo")
  ct<-if(length(r)>1) Consensus(r,p=1) else r[[1]]; Preorder(ct) }
trees<-list(); info<-list()
for(pi in seq_along(picked)){ ct<-buildCT(as.integer(mc[[picked[pi]]])); trees[[pi]]<-ct
  kc<-unique(key(as.Splits(ct,tipLabels=tips))); ov<-mean(kc %in% k1943); nres<-length(kc)
  info[[pi]]<-data.frame(clique=pi, sz=sz[picked[pi]], nsplit=nres, overlap1943=round(ov,3))
  cat(sprintf("clique %2d: %d chars -> %d splits, 1943-overlap=%.3f\n", pi, sz[picked[pi]], nres, ov)) }
saveRDS(list(trees=trees, info=do.call(rbind,info)), "/nobackup/pjjg18/reeval/clique_trees.rds")
cat("SAVED", length(trees), "clique trees\n")
