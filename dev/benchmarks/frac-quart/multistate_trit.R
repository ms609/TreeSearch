# Verify: for a multistate character, GF(2) rank of its decisive quartets
# equals Sum_{i<j}(n_i-1)(n_j-1), and rank of concordant quartets (vs a split)
# equals Sum_{i<j} A_ij.  Blocks are edge-disjoint => trits add across pairs.

gf2rank <- function(M) {                      # rank over GF(2)
  if (nrow(M) == 0 || ncol(M) == 0) return(0L)
  M <- M %% 2L; r <- 0L; nc <- ncol(M)
  for (col in seq_len(nc)) {
    piv <- which(M[, col] == 1L)
    piv <- piv[piv > r]
    if (!length(piv)) next
    p <- piv[1]; r <- r + 1L
    M[c(r, p), ] <- M[c(p, r), ]
    others <- which(M[, col] == 1L); others <- others[others != r]
    if (length(others)) M[others, ] <- sweep(M[others, , drop = FALSE], 2, M[r, ], `+`) %% 2L
    if (r == nrow(M)) break
  }
  r
}

# edge id for unordered pair (a,b) among t taxa
eid <- function(a, b, t) { if (a > b) { tmp <- a; a <- b; b <- tmp }; (a - 1) * t + b }

# All decisive quartets of `char` (state vector), each as its char-induced
# 4-cycle edge vector over C(t,2)... encoded on t*t (sparse-safe) columns.
quartet_rows <- function(char, split = NULL, concordant_only = FALSE) {
  t <- length(char)
  qs <- combn(t, 4)
  rows <- list()
  for (k in seq_len(ncol(qs))) {
    q <- qs[, k]; st <- char[q]
    tb <- table(st)
    if (length(tb) == 2 && all(tb == 2)) {          # "xx yy": decisive
      states <- as.integer(names(tb))
      g1 <- q[st == states[1]]; g2 <- q[st == states[2]]   # the two same-state pairs
      if (concordant_only) {
        s1 <- split[g1]; s2 <- split[g2]
        # concordant iff each same-state pair lies wholly on one split side,
        # and the two pairs are on opposite sides
        ok <- length(unique(s1)) == 1 && length(unique(s2)) == 1 && s1[1] != s2[1]
        if (!ok) next
      }
      v <- integer(t * t)
      for (a in g1) for (b in g2) v[eid(a, b, t)] <- 1L    # cross edges = 4-cycle
      rows[[length(rows) + 1]] <- v
    }
  }
  if (!length(rows)) return(matrix(0L, 0, t * t))
  do.call(rbind, rows)
}

Aij <- function(p, q, r, s) max(p - 1, 0) * max(s - 1, 0) + max(q - 1, 0) * max(r - 1, 0)

set.seed(1)
cat(sprintf("%-28s %6s %6s %6s\n", "case", "rankD", "SumW", "match"))
for (trial in 1:6) {
  t <- sample(6:9, 1)
  r <- sample(2:4, 1)
  char <- sample(seq_len(r), t, replace = TRUE)
  while (length(unique(char)) < 2 || any(tabulate(char) == 1 & tabulate(char) > 0 & FALSE)) char <- sample(seq_len(r), t, replace = TRUE)
  ns <- as.integer(tabulate(char)); ns <- ns[ns > 0]
  sumW <- sum(outer(seq_along(ns), seq_along(ns), Vectorize(function(i, j)
    if (i < j) (ns[i] - 1) * (ns[j] - 1) else 0)))
  rD <- gf2rank(quartet_rows(char))
  cat(sprintf("t=%d states=%-14s %6d %6d %6s\n",
              t, paste(ns, collapse = ","), rD, sumW, rD == sumW))
}

cat("\n-- concordant vs split --\n")
cat(sprintf("%-34s %6s %6s %6s\n", "case", "rankC", "SumAij", "match"))
for (trial in 1:6) {
  t <- sample(6:9, 1); r <- sample(2:4, 1)
  char <- sample(seq_len(r), t, replace = TRUE)
  while (length(unique(char)) < 2) char <- sample(seq_len(r), t, replace = TRUE)
  split <- sample(c(TRUE, FALSE), t, replace = TRUE)
  states <- sort(unique(char))
  sumA <- 0
  for (ii in seq_along(states)) for (jj in seq_along(states)) if (ii < jj) {
    i <- states[ii]; j <- states[jj]
    p <- sum(char == i & split); qc <- sum(char == i & !split)
    rr <- sum(char == j & split); ss <- sum(char == j & !split)
    sumA <- sumA + Aij(p, qc, rr, ss)
  }
  rC <- gf2rank(quartet_rows(char, split, concordant_only = TRUE))
  cat(sprintf("t=%d states=%-12s split=%-2d/%-2d %6d %6d %6s\n",
              t, paste(as.integer(tabulate(char))[tabulate(char) > 0], collapse = ","),
              sum(split), sum(!split), rC, sumA, rC == sumA))
}
