source("gfrank.R")   # reuse gf2rank (prints table; ignore)
# Represent the SAME decisive quartet set two ways:
#  - by the CHARACTER's pairing (4-cycle in K_{n,l}, I x O)      -> Delta_char
#  - by the SPLIT's pairing     (4-cycle in K_{m,t-m}, A x B)     -> Delta_split
# Concordant quartets pair the same either way; discordant differ.
delta_both <- function(p, q, r, s) {
  n <- p + q; l <- r + s; m <- p + r; b <- q + s   # split sizes
  # taxon coords: state (I=1..n grouped p then q), side; build 4 cells
  # cell membership: give each taxon (charState in {I,O}, side in {A,B})
  taxa <- rbind(
    cbind(I=1, A=1, id=seq_len(p)),                       # p: I&A
    cbind(I=1, A=0, id=seq_len(q)),                       # q: I&B
    cbind(I=0, A=1, id=seq_len(r)),                       # r: O&A
    cbind(I=0, A=0, id=seq_len(s)))                       # s: O&B
  N <- nrow(taxa)
  # index taxa within char-graph: I-vertices and O-vertices
  Iidx <- which(taxa[,"I"]==1); Oidx <- which(taxa[,"I"]==0)
  Aidx <- which(taxa[,"A"]==1); Bidx <- which(taxa[,"A"]==0)
  vc <- match(seq_len(N), Iidx); vo <- match(seq_len(N), Oidx)  # char graph coord
  va <- match(seq_len(N), Aidx); vb <- match(seq_len(N), Bidx)  # split graph coord
  edgeC <- function(iA, iB) (iA-1)*length(Oidx) + iB           # char edge id
  edgeS <- function(a, b) (a-1)*length(Bidx) + b               # split edge id
  rowsC <- list(); rowsS <- list()
  quads <- combn(N, 4)
  for (j in seq_len(ncol(quads))) {
    t4 <- quads[, j]
    isI <- taxa[t4,"I"]==1; isA <- taxa[t4,"A"]==1
    if (sum(isI)!=2 || sum(isA)!=2) next                  # decisive = both resolve 2-2
    Is <- t4[isI]; Os <- t4[!isI]; As <- t4[isA]; Bs <- t4[!isA]
    # char pairing: (Is)|(Os)
    ec <- c(edgeC(vc[Is[1]],vo[Os[1]]), edgeC(vc[Is[1]],vo[Os[2]]),
            edgeC(vc[Is[2]],vo[Os[1]]), edgeC(vc[Is[2]],vo[Os[2]]))
    es <- c(edgeS(va[As[1]],vb[Bs[1]]), edgeS(va[As[1]],vb[Bs[2]]),
            edgeS(va[As[2]],vb[Bs[1]]), edgeS(va[As[2]],vb[Bs[2]]))
    vC <- integer(length(Iidx)*length(Oidx)); vC[ec] <- 1; rowsC[[length(rowsC)+1]] <- vC
    vS <- integer(length(Aidx)*length(Bidx)); vS[es] <- 1; rowsS[[length(rowsS)+1]] <- vS
  }
  mk <- function(L, w) if (length(L)) do.call(rbind, L) else matrix(0,0,w)
  c(Dchar = gf2rank(mk(rowsC, length(Iidx)*length(Oidx))),
    Dsplit= gf2rank(mk(rowsS, length(Aidx)*length(Bidx))),
    char_wt = (n-1)*(l-1), split_wt = (m-1)*(b-1))
}
for (cfg in list(c(3,1,2,2), c(2,3,4,1), c(5,2,3,4), c(4,1,2,3))) {
  d <- delta_both(cfg[1],cfg[2],cfg[3],cfg[4])
  cat(sprintf("p,q,r,s=%-10s  Delta_char=%d (charWt=%d)   Delta_split=%d (splitWt=%d)\n",
      paste(cfg,collapse=","), d["Dchar"], d["char_wt"], d["Dsplit"], d["split_wt"]))
}
