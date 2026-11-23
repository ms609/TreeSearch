# B2 min-length-0 CRITERION ORACLE (2026-06-22).
#
# Aggressive collapse-during-search (TNT `collapse 3`) collapses edges of minimum
# possible length 0 = "EXISTS a most-parsimonious reconstruction (MPR) with no
# state change along the edge".  Validate the candidate kernel criterion against
# brute-force ground truth on tiny trees before porting to C++ (advisor: don't
# bless a formula by derivation).
#
# Ground truth (per char, per non-root edge p-c): enumerate ALL internal
# labelings; M = global min total changes; edge is soft-zero-length for that char
# iff SOME labeling achieving M has label(p)==label(c).  Collapsible iff
# soft-zero-length for EVERY char.  Pure enumeration; hand-verified on ((A,B),(C,D)).
#
# Candidate criterion: F_final = final[p] & final[c] != 0 for all chars, where
# final = Fitch up-pass (MPR) set computed with proper intersection/union node
# typing (Swofford & Maddison).  This is what the C++ kernel's `final_` array
# holds, so a match here licenses a direct kernel port.

bit <- function(s) bitwShiftL(1L, s)

randTree <- function(nt) {
  active <- 1:nt; nextNode <- nt + 1L; parent <- integer(2L * nt - 1L)
  while (length(active) > 1L) {
    ij <- sample(length(active), 2L); a <- active[ij[1]]; b <- active[ij[2]]
    parent[a] <- nextNode; parent[b] <- nextNode
    active <- c(active[-ij], nextNode); nextNode <- nextNode + 1L
  }
  root <- 2L * nt - 1L; parent[root] <- 0L
  list(parent = parent, nt = nt, root = root, nNode = 2L * nt - 1L)
}
childrenOf <- function(tr) {
  ch <- vector("list", tr$nNode)
  for (v in seq_len(tr$nNode)) if (tr$parent[v] != 0L)
    ch[[tr$parent[v]]] <- c(ch[[tr$parent[v]]], v)
  ch
}
postVisit <- function(tr, ch) { po <- integer(0)
  rec <- function(v) { for (k in ch[[v]]) rec(k); po[[length(po) + 1L]] <<- v }
  rec(tr$root); unlist(po) }

# --- Brute-force ground truth -------------------------------------------------
oracleEdgeZero <- function(tr, ch, tipMat, nState) {
  nInt <- tr$nt - 1L; intNodes <- (tr$nt + 1L):tr$nNode
  grid <- as.matrix(expand.grid(rep(list(0:(nState - 1L)), nInt)))
  edges <- which(seq_len(tr$nNode) != tr$root)
  zero <- rep(TRUE, length(edges)); nChar <- ncol(tipMat)
  for (cc in seq_len(nChar)) {
    tipState <- tipMat[, cc]
    nc <- apply(grid, 1L, function(g) {
      st <- c(tipState, g)
      sum(vapply(seq_len(tr$nNode), function(v)
        if (tr$parent[v] != 0L && st[v] != st[tr$parent[v]]) 1L else 0L, integer(1)))
    })
    M <- min(nc); optRows <- which(nc == M)
    for (ei in seq_along(edges)) {
      c0 <- edges[ei]; p0 <- tr$parent[c0]; feasible <- FALSE
      for (r in optRows) {
        st <- c(tipState, grid[r, ]); if (st[c0] == st[p0]) { feasible <- TRUE; break }
      }
      if (!feasible) zero[ei] <- FALSE
    }
  }
  setNames(zero, edges)
}

# --- Correct Fitch down + up pass (bitsets), with int/union typing ------------
fitchFinal <- function(tr, ch, tipMat, nState) {
  nChar <- ncol(tipMat)
  prelim <- matrix(0L, tr$nNode, nChar); final <- matrix(0L, tr$nNode, nChar)
  isInt  <- matrix(FALSE, tr$nNode, nChar)        # node formed by intersection?
  po <- postVisit(tr, ch); pre <- rev(po)
  for (cc in seq_len(nChar)) {
    for (v in po) {                                # downpass
      kids <- ch[[v]]
      if (!length(kids)) { prelim[v, cc] <- bit(tipMat[v, cc]); next }
      a <- prelim[kids[1], cc]; b <- prelim[kids[2], cc]; inter <- bitwAnd(a, b)
      if (inter != 0L) { prelim[v, cc] <- inter; isInt[v, cc] <- TRUE }
      else prelim[v, cc] <- bitwOr(a, b)
    }
    for (v in pre) {                               # uppass (Swofford & Maddison)
      if (v == tr$root) { final[v, cc] <- prelim[v, cc]; next }
      kids <- ch[[v]]
      if (!length(kids)) { final[v, cc] <- prelim[v, cc]; next }   # tip: fixed
      p <- tr$parent[v]; fp <- final[p, cc]; sv <- prelim[v, cc]
      if (bitwAnd(fp, sv) == fp) {                 # prelim contains all parent final
        final[v, cc] <- fp
      } else if (isInt[v, cc]) {
        sl <- prelim[kids[1], cc]; sr <- prelim[kids[2], cc]
        final[v, cc] <- bitwOr(sv, bitwAnd(fp, bitwOr(sl, sr)))
      } else {
        final[v, cc] <- bitwOr(sv, fp)
      }
    }
  }
  final
}

# --- Compare F_final to oracle ------------------------------------------------
set.seed(1); nState <- 3L; nChar <- 6L
tot <- 0L; mismatch <- 0L; falsePos <- 0L; falseNeg <- 0L
for (trial in seq_len(400L)) {
  nt <- sample(4:6, 1L); tr <- randTree(nt); ch <- childrenOf(tr)
  tipMat <- matrix(sample(0:(nState - 1L), nt * nChar, replace = TRUE), nt, nChar)
  orc <- oracleEdgeZero(tr, ch, tipMat, nState)
  final <- fitchFinal(tr, ch, tipMat, nState)
  edges <- as.integer(names(orc))
  for (ei in seq_along(edges)) {
    c0 <- edges[ei]; p0 <- tr$parent[c0]
    fF <- all(vapply(seq_len(nChar), function(cc)
      bitwAnd(final[p0, cc], final[c0, cc]) != 0L, logical(1)))
    tot <- tot + 1L
    if (fF != orc[ei]) { mismatch <- mismatch + 1L
      if (fF && !orc[ei]) falsePos <- falsePos + 1L else falseNeg <- falseNeg + 1L }
  }
}
cat(sprintf("Tested %d edges over 400 trees (nt 4-6, %d states, %d chars)\n",
            tot, nState, nChar))
cat(sprintf("F_final (final[p] & final[c] != 0, all chars): %d mismatches (%.2f%%)\n",
            mismatch, 100 * mismatch / tot))
cat(sprintf("  false-positive (formula says collapse, oracle no): %d\n", falsePos))
cat(sprintf("  false-negative (formula no, oracle yes):           %d\n", falseNeg))
cat(if (mismatch == 0L) "=> F_final MATCHES oracle exactly -> port to kernel.\n"
    else if (falsePos == 0L) "=> F_final is CONSERVATIVE (no false collapses; safe, may under-collapse).\n"
    else "=> F_final has FALSE POSITIVES -> would over-collapse; not safe as-is.\n")
