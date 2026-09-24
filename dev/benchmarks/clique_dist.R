suppressWarnings(suppressMessages({library(TreeSearch); library(TreeTools); library(igraph)}))
nex <- commandArgs(trailingOnly=TRUE)[[1]]
pd <- ReadAsPhyDat(nex); tips <- names(pd); nt <- length(tips)
ch <- ReadCharacters(nex); rownames(ch) <- gsub(" ","_",rownames(ch)); ch <- ch[tips,,drop=FALSE]
miss <- c("?","-","N","n")
cl1 <- function(x){x<-as.character(x); ifelse(!(x%in%miss)&nchar(x)==1L&grepl("^[0-9A-Za-z]$",x),x,NA_character_)}
M <- matrix(cl1(ch), nrow=nt)
codeL <- vector("list", ncol(M)); infm <- logical(ncol(M))
for (c in seq_len(ncol(M))){ f<-factor(M[,c]); v<-as.integer(f); codeL[[c]]<-v
  tb<-tabulate(v[!is.na(v)], nlevels(f)); infm[c] <- sum(tb>=2) >= 2 }
inf <- which(infm); ni <- length(inf)
cat(sprintf("chars total=%d informative=%d\n", ncol(M), ni))
compat <- function(a,b){ ok<-!is.na(a)&!is.na(b); a<-a[ok]; b<-b[ok]; if(!length(a)) return(TRUE)
  ed<-unique(cbind(paste0("A",a),paste0("B",b))); nodes<-unique(c(ed[,1],ed[,2]))
  idx<-setNames(seq_along(nodes),nodes); parent<-seq_along(nodes)
  findR<-function(i){ while(parent[i]!=i) i<-parent[i]; i }; cyc<-FALSE
  for(r in seq_len(nrow(ed))){ u<-findR(idx[[ed[r,1]]]); v<-findR(idx[[ed[r,2]]]); if(u==v){cyc<-TRUE;break}; parent[u]<-v }
  !cyc }
el<-list(); k<-0L
for (ii in seq_len(ni-1)){ ai<-codeL[[inf[ii]]]
  for (jj in (ii+1):ni) if (compat(ai, codeL[[inf[jj]]])){ k<-k+1L; el[[k]]<-c(ii,jj) } }
E <- if(k) do.call(rbind, el) else matrix(integer(0), ncol=2)
G <- make_empty_graph(n=ni, directed=FALSE); if(nrow(E)) G <- add_edges(G, as.vector(t(E)))
cat(sprintf("compat graph: %d nodes %d edges density=%.3f\n", ni, nrow(E), if(ni>1) nrow(E)/(ni*(ni-1)/2) else 0))
mc <- max_cliques(G, min=2); sizes <- sort(vapply(mc, length, integer(1)), decreasing=TRUE)
cat("n maximal cliques(>=2):", length(mc), "\n")
cat("clique-size dist (top 30 chars/clique):", paste(head(sizes,30),collapse=" "), "\n")
cat(sprintf("largest=%d 2nd=%s 3rd=%s median=%.1f  n>=10=%d n>=20=%d n>=40=%d\n",
  sizes[1], ifelse(length(sizes)>1,sizes[2],NA), ifelse(length(sizes)>2,sizes[3],NA),
  median(sizes), sum(sizes>=10), sum(sizes>=20), sum(sizes>=40)))
big <- mc[[which.max(vapply(mc,length,integer(1)))]]; cc <- inf[as.integer(big)]; spl<-0L
for (c in cc){ v<-codeL[[c]]; for (s in unique(v[!is.na(v)])){ sd<-sum(!is.na(v)&v==s); if(sd>=2 && sd<=nt-2) spl<-spl+1L } }
cat(sprintf("largest clique = %d chars -> up to %d nontrivial splits (of %d needed for full resolution)\n", length(big), spl, nt-3))
