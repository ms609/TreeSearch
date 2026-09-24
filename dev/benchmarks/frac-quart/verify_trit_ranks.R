# Independent-quartet counts as GF(2) ranks in the character's bipartite
# cycle space. A character quartet {i1,i2 | o1,o2} = a 4-cycle of K_{n,l}
# (n = I-taxa, l = O-taxa); entailment = cycle addition (Nelson-Ladiges).
gf2rank <- function(M) {                       # rank over GF(2)
  if (!nrow(M) || !ncol(M)) return(0L)
  M <- M %% 2; r <- 0L; ncol <- ncol(M); row <- 1L
  for (col in seq_len(ncol)) {
    piv <- which(M[row:nrow(M), col] == 1)
    if (!length(piv)) next
    piv <- piv[1] + row - 1L
    tmp <- M[row, ]; M[row, ] <- M[piv, ]; M[piv, ] <- tmp
    others <- which(M[, col] == 1); others <- others[others != row]
    for (o in others) M[o, ] <- (M[o, ] + M[row, ]) %% 2
    row <- row + 1L; r <- r + 1L
    if (row > nrow(M)) break
  }
  r
}

# taxa: I-vertices 1..n (in A iff <=p), O-vertices 1..l (in A iff <=r)
quartet_vec <- function(i1, i2, o1, o2, l) {
  v <- integer(0)  # edge id = (i-1)*l + o
  e <- c((i1-1)*l+o1, (i1-1)*l+o2, (i2-1)*l+o1, (i2-1)*l+o2)
  e
}
counts <- function(p, q, r, s) {
  n <- p + q; l <- r + s; nEdge <- n * l
  inA_I <- function(i) i <= p
  inA_O <- function(o) o <= r
  rows_conc <- list(); rows_disc <- list(); rows_dec <- list()
  for (i1 in 1:(n-1)) for (i2 in (i1+1):n) for (o1 in 1:(l-1)) for (o2 in (o1+1):l) {
    nA <- inA_I(i1) + inA_I(i2) + inA_O(o1) + inA_O(o2)
    if (nA != 2) next                              # split doesn't resolve -> not decisive
    e <- c((i1-1)*l+o1, (i1-1)*l+o2, (i2-1)*l+o1, (i2-1)*l+o2)
    v <- integer(nEdge); v[e] <- 1
    # concordant iff the two I-taxa are on the same split side
    conc <- (inA_I(i1) == inA_I(i2))
    rows_dec[[length(rows_dec)+1]] <- v
    if (conc) rows_conc[[length(rows_conc)+1]] <- v
    else       rows_disc[[length(rows_disc)+1]] <- v
  }
  mk <- function(L) if (length(L)) do.call(rbind, L) else matrix(0, 0, nEdge)
  list(conc = gf2rank(mk(rows_conc)),
       disc = gf2rank(mk(rows_disc)),
       dec  = gf2rank(mk(rows_dec)),
       n_conc = length(rows_conc), n_disc = length(rows_disc))
}
pf <- function(x) pmax(x - 1, 0)
cat(sprintf("%-14s %5s %5s %5s | %-22s %-10s\n",
            "p,q,r,s", "conc", "disc", "dec", "A=(p1)(s1)+(q1)(r1)", "n-1 l-1"))
for (cfg in list(c(2,1,1,2), c(3,2,2,3), c(3,1,2,2), c(2,2,2,2),
                 c(4,1,1,4), c(3,3,1,1), c(2,3,4,1), c(5,2,3,4),
                 c(1,2,3,4), c(4,0,0,4), c(3,0,2,3), c(2,1,0,3))) {
  p<-cfg[1];q<-cfg[2];r<-cfg[3];s<-cfg[4]; n<-p+q; l<-r+s
  x <- counts(p,q,r,s)
  A <- pf(p)*pf(s) + pf(q)*pf(r)
  cat(sprintf("%-14s %5d %5d %5d | A=%-3d disc? dec vs %-4d| char=(%d)(%d)=%d\n",
      paste(cfg,collapse=","), x$conc, x$disc, x$dec, A, A + x$disc, n-1, l-1, (n-1)*(l-1)))
}
