suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(igraph)}))
nex <- "/nobackup/pjjg18/lgsweep/matrices/project5432.nex"
FLOOR <- "/nobackup/pjjg18/floors/project5432_tnt_floor_one.tre"
pd <- ReadAsPhyDat(nex); tips <- names(pd); nt <- length(tips)
ch <- ReadCharacters(nex); rownames(ch) <- gsub(" ","_",rownames(ch)); ch <- ch[tips,,drop=FALSE]
miss <- c("?","-","N","n")
cl1 <- function(x){x<-as.character(x); ifelse(!(x%in%miss)&nchar(x)==1L&grepl("^[0-9A-Za-z]$",x),x,NA_character_)}
M <- matrix(cl1(ch), nrow=nt)
codeL <- vector("list", ncol(M)); infm <- logical(ncol(M))
for (c in seq_len(ncol(M))){ f<-factor(M[,c]); v<-as.integer(f); codeL[[c]]<-v
  tb<-tabulate(v[!is.na(v)], nlevels(f)); infm[c] <- sum(tb>=2) >= 2 }
inf <- which(infm); ni <- length(inf)
compat <- function(a,b,minOv){ ok<-!is.na(a)&!is.na(b); if(sum(ok)<minOv) return(FALSE); a<-a[ok]; b<-b[ok]
  ed<-unique(cbind(paste0("A",a),paste0("B",b))); nodes<-unique(c(ed[,1],ed[,2]))
  idx<-setNames(seq_along(nodes),nodes); parent<-seq_along(nodes)
  findR<-function(i){ while(parent[i]!=i) i<-parent[i]; i }; cyc<-FALSE
  for(r in seq_len(nrow(ed))){ u<-findR(idx[[ed[r,1]]]); v<-findR(idx[[ed[r,2]]]); if(u==v){cyc<-TRUE;break}; parent[u]<-v }
  !cyc }
buildG <- function(minOv){ el<-list(); k<-0L
  for (ii in seq_len(ni-1)){ ai<-codeL[[inf[ii]]]
    for (jj in (ii+1):ni) if (compat(ai, codeL[[inf[jj]]], minOv)){ k<-k+1L; el[[k]]<-c(ii,jj) } }
  E <- if(k) do.call(rbind, el) else matrix(integer(0),ncol=2)
  G <- make_empty_graph(n=ni, directed=FALSE); if(nrow(E)) G<-add_edges(G, as.vector(t(E))); G }
# 1943 splits as membership sets (side not containing tip1)
tnt <- ReadTntTree(FLOOR, tipLabels=tips); if(inherits(tnt,"multiPhylo")) tnt<-tnt[[1]]; tnt<-Preorder(tnt)
sp <- as.Splits(tnt, tipLabels=tips); smat <- as.matrix(sp)
key1943 <- apply(smat, 1, function(r){ v<-as.logical(r); if(v[1]) v<-!v; paste(which(v),collapse=",") })
key1943 <- unique(key1943)
cliqueKeys <- function(cc){ ks<-character(0); for(c in inf[cc]){ v<-codeL[[c]]
   for(s in unique(v[!is.na(v)])){ side<-!is.na(v)&v==s; nsd<-sum(side); if(nsd>=2 && nsd<=nt-2){
     w<-side; if(w[1]) w<-!w; ks<-c(ks, paste(which(w),collapse=",")) } } }; unique(ks) }
for (mo in c(0, 20, 40)){
  G <- buildG(mo); ne <- ecount(G); dens <- ne/(ni*(ni-1)/2)
  cn <- clique_num(G); ncq <- tryCatch(count_max_cliques(G, min=2), error=function(e) NA)
  lc <- largest_cliques(G)
  ov <- vapply(lc, function(cl){ ks<-cliqueKeys(as.integer(cl)); if(!length(ks)) return(0); mean(ks %in% key1943) }, numeric(1))
  nsplit <- vapply(lc, function(cl) length(cliqueKeys(as.integer(cl))), integer(1))
  cat(sprintf("minOverlap=%2d: edges=%d dens=%.3f | largestClique=%d chars nMaxCliques=%s | maxCliques induce %d-%d splits, 1943-overlap max=%.2f mean=%.2f\n",
      mo, ne, dens, cn, as.character(ncq), min(nsplit), max(nsplit), max(ov), mean(ov)))
}
cat(sprintf("(1943 has %d nontrivial splits)\n", length(key1943)))
