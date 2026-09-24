# Brute-force check of concordant/discordant/independent-quartet formulae
# for two binary bipartitions with 2x2 taxon table {p,q,r,s}.
choose2 <- function(x) x * (x - 1) / 2

# quartet resolution induced by a 0/1 labelling of 4 taxa: returns the
# unordered pair-of-pairs if 2-2, else NA
res <- function(lab) {
  if (sum(lab) != 2) return(NA_character_)   # not 2-2
  paste(sort(which(lab == lab[order(lab)][1]))[1:2], collapse = "")  # placeholder
}
# cleaner: canonical partition key for a 2-2 split of positions 1:4
key22 <- function(side) {
  if (length(unique(side)) != 2 || sum(side == side[1]) != 2) return(NA_character_)
  a <- sort(which(side == side[1]))
  paste(a, collapse = "")   # the pair on 'first' side identifies the 2|2 split up to complement
}

set.seed(1)
ok <- TRUE
for (rep in 1:2000) {
  t <- sample(6:11, 1)
  charState <- sample(0:1, t, replace = TRUE)   # I = state1
  splitSide <- sample(0:1, t, replace = TRUE)   # A = side1
  p <- sum(charState == 1 & splitSide == 1)
  q <- sum(charState == 1 & splitSide == 0)
  r <- sum(charState == 0 & splitSide == 1)
  s <- sum(charState == 0 & splitSide == 0)

  conc <- 0; disc <- 0; decisive <- 0
  combs <- combn(t, 4)
  for (j in seq_len(ncol(combs))) {
    idx <- combs[, j]
    ck <- key22(charState[idx]); sk <- key22(splitSide[idx])
    if (!is.na(ck)) decisive <- decisive + 1
    if (!is.na(ck) && !is.na(sk)) {
      # same 2|2 split iff the pair-key matches OR is complementary
      same <- (ck == sk) || (ck == paste(setdiff(1:4, as.integer(strsplit(ck,"")[[1]])), collapse=""))
      if (same) conc <- conc + 1 else disc <- disc + 1
    }
  }
  fConc <- choose2(p)*choose2(s) + choose2(q)*choose2(r)
  fDisc <- p*q*r*s
  n <- p + q
  fDecisive <- choose2(n) * choose2(t - n)
  if (conc != fConc || disc != fDisc || decisive != fDecisive) {
    cat("MISMATCH pqrs=", p,q,r,s,
        " conc", conc, fConc, " disc", disc, fDisc,
        " dec", decisive, fDecisive, "\n"); ok <- FALSE; break
  }
  # fractional (independent) agreement with floors
  pf <- function(x) pmax(x - 1, 0)
  indConc <- pf(p)*pf(s) + pf(q)*pf(r)
  indDec  <- pf(n)*pf(t - n)           # char's own independent quartets
  if (indDec > 0) {
    FQ <- indConc / indDec
    if (FQ < -1e-9 || FQ > 1 + 1e-9) { cat("FQ out of range", FQ, p,q,r,s, "\n"); ok <- FALSE; break }
  }
  # self-agreement: char == split  => q=r=0 => FQ==1
  if (q == 0 && r == 0 && indDec > 0 && abs(FQ - 1) > 1e-9) {
    cat("self-agree != 1", FQ, p,q,r,s, "\n"); ok <- FALSE; break
  }
}
cat(if (ok) "ALL CHECKS PASS\n" else "FAILED\n")

# ind_decisive - ind_concordant identity (advisor: = pr + qs - 1 without floors)
set.seed(2); id_ok <- TRUE
for (i in 1:5000) {
  p<-sample(0:6,1);q<-sample(0:6,1);r<-sample(0:6,1);s<-sample(0:6,1)
  n<-p+q; tt<-p+q+r+s
  lhs <- (n-1)*(tt-n-1) - ((p-1)*(s-1) + (q-1)*(r-1))
  if (lhs != p*r + q*s - 1) { cat("identity fail\n"); id_ok<-FALSE; break }
}
cat(if (id_ok) "IDENTITY (no-floor) ind_dec - ind_conc == pr+qs-1: CONFIRMED\n" else "IDENTITY FAIL\n")
