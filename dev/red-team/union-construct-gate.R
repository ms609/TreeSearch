# Union-construct lower-bound gate (reproducible, pure R, no C++ engine needed).
# See dev/red-team/proofs/union-construct-lower-bound.md for the verdict.
# Run:  Rscript dev/red-team/union-construct-gate.R
#
# Tests Goloboff (1996) "union construct" screen:
#   UC(D) = union_set(D)  UNION  MP(A)          (A = parent of D)
# CLAIM (lower-bound / coverage): for EVERY descendant edge D' in subtree(D)
# (including D itself), the exact MPR edge set M(D') = combine(P(D'),Up(D'))
# satisfies  M(D') SUBSET-OF UC(D).  Antitone disjointness indicator then makes
#   added_UC(D) = |{c : r_c disjoint UC(D)_c}|  <=  added_exact(D') for all r.
#
# We test the CONTAINMENT directly (strongest form; holds for all clip states r):
#   variant TRUE : MP(A) = M(A)  = exact MPR finals  (fitchUp)
#   variant COLL : MP(A) = fin(A)= engine uppass_node collapse (fitchUpEngine)
# Single unordered (Fitch) character at a time; state s -> bit s.  Polymorphic
# tips = multi-bit masks.  Missing "?" = full mask.  (Brazeau inapplicable = NA
# is handled analytically in the proof; see its NA section.)

suppressMessages({ library(ape) })

fullMask <- function(k) bitwShiftL(1L, k) - 1L
combineMask <- function(a, b) { i <- bitwAnd(a,b); ifelse(i != 0L, i, bitwOr(a,b)) }

getRoot <- function(tr) setdiff(unique(tr$edge[,1]), unique(tr$edge[,2]))
childrenOf <- function(e, n) e[e[,1]==n, 2]
parentOf <- function(e, n) e[e[,2]==n, 1]
siblingOf <- function(e, n) { p <- parentOf(e,n); setdiff(childrenOf(e,p), n) }

preorderNodes <- function(tr) {
  root <- getRoot(tr); e <- tr$edge; out <- integer(0); st <- root
  while (length(st)) { nd <- st[length(st)]; st <- st[-length(st)]
    out <- c(out, nd); ch <- childrenOf(e, nd); if (length(ch)) st <- c(st, ch) }
  out
}
postorderInternal <- function(tr) { pre <- rev(preorderNodes(tr)); pre[pre > length(tr$tip.label)] }

# Down-pass: prelim P per node (single char, integer mask vector indexed by node).
fitchDown <- function(tr, tipMask) {
  nTip <- length(tr$tip.label); nNode <- nTip + tr$Nnode
  P <- integer(nNode); P[seq_len(nTip)] <- tipMask
  for (nd in postorderInternal(tr)) {
    ch <- childrenOf(tr$edge, nd); acc <- P[ch[1]]
    if (length(ch) >= 2) for (j in 2:length(ch)) acc <- combineMask(acc, P[ch[j]])
    P[nd] <- acc
  }
  P
}
# union_set per node: recursive UNION (not combine) of descendant tip masks.
unionSets <- function(tr, tipMask) {
  nTip <- length(tr$tip.label); nNode <- nTip + tr$Nnode
  U <- integer(nNode); U[seq_len(nTip)] <- tipMask
  for (nd in postorderInternal(tr)) {
    ch <- childrenOf(tr$edge, nd); U[nd] <- Reduce(bitwOr, U[ch], 0L)
  }
  U
}
# Exact directional up-pass: Up[nd] = message into nd; E[nd]=combine(P,Up)=M(nd).
fitchUp <- function(tr, P) {
  nTip <- length(tr$tip.label); nNode <- nTip + tr$Nnode; e <- tr$edge; root <- getRoot(tr)
  Up <- integer(nNode); E <- integer(nNode); E[root] <- P[root]
  for (nd in preorderNodes(tr)) {
    if (nd == root) next
    A <- parentOf(e, nd); sib <- siblingOf(e, nd)
    if (A == root) { u <- P[sib[1]]; if (length(sib) > 1) for (k in 2:length(sib)) u <- combineMask(u, P[sib[k]]) }
    else { u <- combineMask(Up[A], P[sib[1]]); if (length(sib) > 1) for (k in 2:length(sib)) u <- combineMask(u, P[sib[k]]) }
    Up[nd] <- u; E[nd] <- combineMask(P[nd], u)
  }
  list(Up = Up, E = E)
}
# Engine's COLLAPSED up-pass (uppass_node): fin(root)=P; fin(nd)=fin(A)&P[nd] if
# non-empty else P[nd].  Always a SUBSET of P[nd].
fitchUpEngine <- function(tr, P) {
  nTip <- length(tr$tip.label); nNode <- nTip + tr$Nnode; e <- tr$edge; root <- getRoot(tr)
  Fe <- integer(nNode); Fe[root] <- P[root]
  for (nd in preorderNodes(tr)) {
    if (nd == root) next
    A <- parentOf(e, nd); is <- bitwAnd(Fe[A], P[nd]); Fe[nd] <- ifelse(is != 0L, is, P[nd])
  }
  Fe
}
# Brute-force MPR finals for a single UNAMBIGUOUS char (validates reference).
bruteFinals <- function(tr, tipState, k) {
  nTip <- length(tr$tip.label); internal <- (nTip+1):(nTip+tr$Nnode); e <- tr$edge
  grid <- as.matrix(expand.grid(rep(list(0:(k-1)), length(internal))))
  full <- integer(nTip + tr$Nnode); full[seq_len(nTip)] <- tipState
  lens <- integer(nrow(grid))
  for (g in seq_len(nrow(grid))) { full[internal] <- grid[g,]; lens[g] <- sum(full[e[,1]] != full[e[,2]]) }
  opt <- which(lens == min(lens))
  fm <- integer(nTip + tr$Nnode); fm[seq_len(nTip)] <- bitwShiftL(1L, tipState)
  for (i in seq_along(internal)) fm[internal[i]] <- Reduce(bitwOr, bitwShiftL(1L, unique(grid[opt,i])), 0L)
  fm
}

# subtree node membership (nodes strictly below D, plus D) — the descendant edges.
subtreeNodes <- function(tr, D) {
  e <- tr$edge; out <- integer(0); st <- D
  while (length(st)) { nd <- st[length(st)]; st <- st[-length(st)]; out <- c(out, nd)
    ch <- childrenOf(e, nd); if (length(ch)) st <- c(st, ch) }
  out
}

# ---- Core check for one (tree, single-char tip mask) ----
checkOne <- function(tr, tipMask) {
  nTip <- length(tr$tip.label); root <- getRoot(tr); e <- tr$edge
  P <- fitchDown(tr, tipMask); U <- unionSets(tr, tipMask)
  up <- fitchUp(tr, P); E <- up$E; Up <- up$Up
  Fe <- fitchUpEngine(tr, P)
  vio_true <- 0L; vio_coll <- 0L; vio_true_baseedge <- 0L; vio_true_root <- 0L
  nodes <- setdiff(seq_len(nTip + tr$Nnode), root)
  for (D in nodes) {
    A <- parentOf(e, D)
    UCt <- bitwOr(U[D], E[A])   # MP(A) = exact MPR final M(A)  (E[root]=P[root])
    UCc <- bitwOr(U[D], Fe[A])  # MP(A) = collapsed fin(A)
    dl <- setdiff(subtreeNodes(tr, D), root)   # descendant edges incl. D
    for (Dp in dl) {
      miss_t <- bitwAnd(E[Dp], bitwNot(UCt)) != 0L
      miss_c <- bitwAnd(E[Dp], bitwNot(UCc)) != 0L
      if (miss_t) { vio_true <- vio_true + 1L
        if (Dp == D) vio_true_baseedge <- vio_true_baseedge + 1L
        if (A == root) vio_true_root <- vio_true_root + 1L }
      if (miss_c) vio_coll <- vio_coll + 1L
    }
  }
  c(vio_true = vio_true, vio_coll = vio_coll,
    vio_true_baseedge = vio_true_baseedge, vio_true_root = vio_true_root)
}

# ---- exhaustive single-state regime ----
runExhaustive <- function(nTip, k, nTrees, seed) {
  set.seed(seed)
  assigns <- as.matrix(expand.grid(rep(list(0:(k-1)), nTip)))
  tot <- c(vio_true=0L, vio_coll=0L, vio_true_baseedge=0L, vio_true_root=0L)
  nchecks <- 0L
  for (t in seq_len(nTrees)) {
    tr <- ape::rtree(nTip); tr$edge.length <- NULL
    tr$tip.label <- paste0("t", seq_len(nTip))
    for (g in seq_len(nrow(assigns))) {
      tm <- bitwShiftL(1L, as.integer(assigns[g,]))
      tot <- tot + checkOne(tr, tm); nchecks <- nchecks + 1L
    }
  }
  list(label=sprintf("exhaustive nTip=%d k=%d (%d trees x %d assigns)", nTip, k, nTrees, nrow(assigns)),
       tot=tot, ncharTrees=nchecks)
}

# ---- random polymorphic / missing regime ----
runPoly <- function(nTip, k, nCharTrees, seed, pPoly=0.4, pMiss=0.1) {
  set.seed(seed)
  tot <- c(vio_true=0L, vio_coll=0L, vio_true_baseedge=0L, vio_true_root=0L)
  for (i in seq_len(nCharTrees)) {
    tr <- ape::rtree(nTip); tr$edge.length <- NULL; tr$tip.label <- paste0("t", seq_len(nTip))
    tm <- integer(nTip)
    for (j in seq_len(nTip)) {
      u <- runif(1)
      if (u < pMiss) tm[j] <- fullMask(k)
      else if (u < pMiss + pPoly) { # polymorphic: random nonempty subset (>=1 bit)
        bits <- which(as.logical(rbinom(k, 1, 0.5))) - 1L
        if (!length(bits)) bits <- sample(0:(k-1), 1)
        tm[j] <- Reduce(bitwOr, bitwShiftL(1L, bits), 0L)
      } else tm[j] <- bitwShiftL(1L, sample(0:(k-1), 1))
    }
    tot <- tot + checkOne(tr, tm)
  }
  list(label=sprintf("polymorphic nTip=%d k=%d (%d char-trees, pPoly=%.2f pMiss=%.2f)",
                      nTip, k, nCharTrees, pPoly, pMiss), tot=tot, ncharTrees=nCharTrees)
}

# ---- validate reference Fitch: E[D]=combine(P,Up) is the EXACT per-edge added-
# length set.  Brute check: for a single-state clip r, physically attach a tip of
# mask r as sister of node D (new node W splits edge (A,D)); the Fitch length
# increase must equal [ r disjoint E[D] ].  This validates E as the object the
# proof uses, independent of the up-pass code. ----
fitchLen1 <- function(tr, tipMask) {  # total length, single char
  nTip <- length(tr$tip.label); P <- integer(nTip + tr$Nnode); P[seq_len(nTip)] <- tipMask
  len <- 0L
  for (nd in postorderInternal(tr)) {
    ch <- childrenOf(tr$edge, nd); acc <- P[ch[1]]
    if (length(ch) >= 2) for (j in 2:length(ch)) { b <- P[ch[j]]; i <- bitwAnd(acc,b)
      if (i != 0L) acc <- i else { acc <- bitwOr(acc,b); len <- len + 1L } }
    P[nd] <- acc
  }
  len
}
attachTipAbove <- function(tr, D, rMaskLabel) {  # new tip 'rtip' as sister of D
  nTip <- length(tr$tip.label); e <- tr$edge; A <- parentOf(e, D)
  newTip <- nTip + 1L                       # relabel: shift internals by +1
  shift <- function(x) ifelse(x <= nTip, x, x + 1L)
  e2 <- cbind(shift(e[,1]), shift(e[,2]))
  W <- (nTip + 1L) + tr$Nnode + 1L          # new internal node
  Ds <- shift(D); As <- shift(A)
  e2 <- e2[!(e2[,1]==As & e2[,2]==Ds), , drop=FALSE]
  e2 <- rbind(e2, c(As, W), c(W, Ds), c(W, newTip))
  structure(list(edge = e2, tip.label = c(tr$tip.label, "rtip"),
                 Nnode = tr$Nnode + 1L), class = "phylo")
}
validateReference <- function(nTip, k, nTrees, seed) {
  set.seed(seed); mism <- 0L; ntest <- 0L
  for (t in seq_len(nTrees)) {
    tr <- ape::rtree(nTip); tr$edge.length <- NULL; tr$tip.label <- paste0("t", seq_len(nTip))
    root <- getRoot(tr); e <- tr$edge
    for (g in seq_len(20)) {
      tm <- bitwShiftL(1L, sample(0:(k-1), nTip, replace=TRUE))
      L0 <- fitchLen1(tr, tm)
      P <- fitchDown(tr, tm); E <- fitchUp(tr, P)$E
      r <- bitwShiftL(1L, sample(0:(k-1), 1))
      for (D in setdiff(seq_len(nTip+tr$Nnode), root)) {
        aug <- attachTipAbove(tr, D, r)
        # tip masks for augmented tree: internal shift means original tip i keeps
        # index i (<=nTip), new tip is nTip+1
        augMask <- c(tm, r)
        addBrute <- fitchLen1(aug, augMask) - L0
        addPred  <- as.integer(bitwAnd(r, E[D]) == 0L)
        ntest <- ntest + 1L; if (addBrute != addPred) mism <- mism + 1L
      }
    }
  }
  cat(sprintf("Reference validation (E==brute per-edge added length): %d mismatches / %d edges\n", mism, ntest))
}

# ---- added-length cross-check: containment => inequality for all clip masks r.
# Sample random clip basal masks r and confirm, for every node D, that
#   added_UC_true(D) <= added_exact(D')  for every descendant edge D' in subtree(D).
# (added_exact(D') = disjointCount(r, E[D']) is the exact added Fitch length; see
# sibling proof S1.)  Reports max positive (over-reject) excess; expect 0. ----
addedCheck <- function(nTip, k, nCharTrees, seed) {
  set.seed(seed); worst <- 0L; nbad <- 0L; nchk <- 0L
  for (i in seq_len(nCharTrees)) {
    tr <- ape::rtree(nTip); tr$edge.length <- NULL; tr$tip.label <- paste0("t", seq_len(nTip))
    root <- getRoot(tr); e <- tr$edge
    tm <- bitwShiftL(1L, sample(0:(k-1), nTip, replace=TRUE))
    P <- fitchDown(tr, tm); U <- unionSets(tr, tm); E <- fitchUp(tr, P)$E
    r <- bitwShiftL(1L, sample(0:(k-1), 1))   # single-state clip basal set
    for (D in setdiff(seq_len(nTip+tr$Nnode), root)) {
      A <- parentOf(e, D); UCt <- bitwOr(U[D], E[A])
      addUC <- as.integer(bitwAnd(r, UCt) == 0L)
      for (Dp in setdiff(subtreeNodes(tr, D), root)) { nchk <- nchk + 1L
        addEx <- as.integer(bitwAnd(r, E[Dp]) == 0L)
        if (addUC > addEx) { nbad <- nbad + 1L; worst <- max(worst, addUC - addEx) }
      }
    }
  }
  cat(sprintf("Added-length check (UC_true <= exact, all D'): %d over-rejections / %d checks (worst excess %d)\n",
              nbad, nchk, worst))
}

cat("== Union-construct lower-bound gate ==\n\n")
validateReference(5, 3, 20, 11)
validateReference(6, 3, 12, 12)
addedCheck(6, 4, 4000, 21)
addedCheck(8, 5, 2000, 22)

regimes <- list(
  runExhaustive(4, 3, 40, 101),
  runExhaustive(5, 3, 40, 102),
  runExhaustive(5, 4, 20, 103),
  runExhaustive(6, 4, 12, 104),
  runExhaustive(7, 3, 8,  105),
  runPoly(6, 4, 4000, 201),
  runPoly(7, 4, 3000, 202),
  runPoly(8, 5, 2000, 203)
)
cat("\n")
for (r in regimes) {
  cat(sprintf("%-58s : vio_TRUE=%d  vio_COLL=%d  (of which base-edge D'=D: %d, A=root: %d)\n",
      r$label, r$tot["vio_true"], r$tot["vio_coll"], r$tot["vio_true_baseedge"], r$tot["vio_true_root"]))
}
cat("\nInterpretation: vio_TRUE=0  => UC with EXACT MPR finals M(A) is a sound coverage/lower bound.\n")
cat("                vio_COLL>0  => UC with COLLAPSED uppass_node fin(A) is NOT sound.\n")
