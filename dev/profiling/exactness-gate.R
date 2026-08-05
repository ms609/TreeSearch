# Exactness gate harness (reproducible). See exactness-gate.md for the verdict.
# Run:  Rscript dev/profiling/exactness-gate.R
# Requires: TreeSearch, phangorn, ape, TreeTools. Pure R; independent of the C++ engine.

# Helper library: from-scratch bitmask Fitch (down + up), oracle wrappers,
# brute-force finals, prune/regraft.  Equal-weights Fitch, unordered characters.
# States are 0-indexed; bitmask for state s is bit (s).  Missing/"?"/"-" = all bits.

suppressMessages({
  library(TreeSearch)
  library(phangorn)
  library(ape)
})

fullMask <- function(nStates) bitwShiftL(1L, nStates) - 1L

stateToMask <- function(x, nStates) {
  if (is.na(x) || x == "?" || x == "-") return(fullMask(nStates))
  bitwShiftL(1L, as.integer(x))
}

# charMat: character matrix, rows = tips (rownames = tip labels), cols = chars.
masksFromChars <- function(charMat, nStates) {
  m <- matrix(0L, nrow(charMat), ncol(charMat))
  for (i in seq_len(nrow(charMat)))
    for (j in seq_len(ncol(charMat)))
      m[i, j] <- stateToMask(charMat[i, j], nStates)
  m
}

# Build a tip-mask matrix whose ROW i corresponds to NODE i of `tree`
# (tree$tip.label[i]), aligning arbitrary charMat row order to the tree's tip
# numbering.  ape::rtree() permutes tip labels, so this alignment is essential.
alignedTipMask <- function(tree, charMat, nStates) {
  stopifnot(all(tree$tip.label %in% rownames(charMat)))
  masksFromChars(charMat[tree$tip.label, , drop = FALSE], nStates)
}
# Aligned integer states (0-indexed; NA for missing "?"/"-"), row i = node i.
alignedStates <- function(tree, charMat) {
  cm <- charMat[tree$tip.label, , drop = FALSE]
  out <- suppressWarnings(matrix(as.integer(cm), nrow(cm), ncol(cm)))
  out
}

getRoot <- function(tree) setdiff(unique(tree$edge[, 1]), unique(tree$edge[, 2]))
childrenOf <- function(edge, node) edge[edge[, 1] == node, 2]
parentOf <- function(edge, node) edge[edge[, 2] == node, 1]
siblingOf <- function(edge, node) {
  p <- parentOf(edge, node)
  setdiff(childrenOf(edge, p), node)
}

# Preorder node list (parents strictly before descendants).
preorderNodes <- function(tree) {
  root <- getRoot(tree); edge <- tree$edge
  out <- integer(0); stack <- root
  while (length(stack)) {
    nd <- stack[length(stack)]; stack <- stack[-length(stack)]
    out <- c(out, nd)
    ch <- childrenOf(edge, nd)
    if (length(ch)) stack <- c(stack, ch)
  }
  out
}
# Postorder internal nodes (all children before their parent).
postorderInternal <- function(tree) {
  pre <- preorderNodes(tree)
  r <- rev(pre)
  r[r > length(tree$tip.label)]
}

combineMask <- function(a, b) {
  inter <- bitwAnd(a, b)
  ifelse(inter != 0L, inter, bitwOr(a, b))
}

# Fitch down-pass: returns prelim (per node) and total length.  Handles >2
# children by sequential folding (exact for binary; folding for polytomies).
fitchDown <- function(tree, tipMask) {
  nTip <- length(tree$tip.label)
  nNode <- nTip + tree$Nnode
  nChar <- ncol(tipMask)
  prelim <- matrix(0L, nNode, nChar)
  prelim[seq_len(nTip), ] <- tipMask
  len <- integer(nChar)
  for (nd in postorderInternal(tree)) {
    ch <- childrenOf(tree$edge, nd)
    acc <- prelim[ch[1], ]
    if (length(ch) >= 2) for (j in 2:length(ch)) {
      b <- prelim[ch[j], ]
      inter <- bitwAnd(acc, b)
      hit <- inter != 0L
      acc <- ifelse(hit, inter, bitwOr(acc, b))
      len <- len + as.integer(!hit)
    }
    prelim[nd, ] <- acc
  }
  list(prelim = prelim, length = sum(len), perChar = len)
}

# Fitch up-pass (rooted, degree-2 root): up[nd] = message into nd from parent
# side; final[nd] = combine(prelim[nd], up[nd]).  root final = prelim (unused
# for edges).
fitchUp <- function(tree, prelim) {
  nTip <- length(tree$tip.label)
  nNode <- nTip + tree$Nnode
  nChar <- ncol(prelim)
  edge <- tree$edge
  root <- getRoot(tree)
  up <- matrix(0L, nNode, nChar)
  final <- matrix(0L, nNode, nChar)
  final[root, ] <- prelim[root, ]
  for (nd in preorderNodes(tree)) {
    if (nd == root) next
    A <- parentOf(edge, nd)
    sib <- siblingOf(edge, nd)
    if (A == root) {
      # degree-2 root: message is sibling's prelim
      upnd <- prelim[sib[1], ]
      if (length(sib) > 1) for (k in 2:length(sib)) upnd <- combineMask(upnd, prelim[sib[k], ])
    } else {
      upnd <- combineMask(up[A, ], prelim[sib[1], ])
      if (length(sib) > 1) for (k in 2:length(sib)) upnd <- combineMask(upnd, prelim[sib[k], ])
    }
    up[nd, ] <- upnd
    final[nd, ] <- combineMask(prelim[nd, ], upnd)
  }
  list(up = up, final = final)
}

# Engine's SIMPLIFIED up-pass (replicates uppass_node / fitch_uppass in
# src/ts_fitch.cpp).  Top-down: final(root) = prelim(root); for each other node
# final(node) = final(parent) & prelim(node) if non-empty, else prelim(node).
# Never unions -> engine_final(node) is always a SUBSET of prelim(node).
fitchUpEngine <- function(tree, prelim) {
  nTip <- length(tree$tip.label)
  nNode <- nTip + tree$Nnode
  nChar <- ncol(prelim)
  edge <- tree$edge
  root <- getRoot(tree)
  final <- matrix(0L, nNode, nChar)
  final[root, ] <- prelim[root, ]
  for (nd in preorderNodes(tree)) {
    if (nd == root) next
    A <- parentOf(edge, nd)
    isect <- bitwAnd(final[A, ], prelim[nd, ])
    final[nd, ] <- ifelse(isect != 0L, isect, prelim[nd, ])
  }
  list(final = final)
}

# Regraft subtree S onto the edge above node D of `base` (a new node w splits
# that edge; w's children = D and S's root).  Explicit edge-matrix surgery, then
# normalise numbering via write/read newick.  Branch lengths dropped (Fitch is
# length-invariant).  Returns a phylo.
regraft <- function(base, S, D) {
  nTb <- length(base$tip.label); nTs <- length(S$tip.label)
  Nb <- base$Nnode; Ns <- S$Nnode
  nT <- nTb + nTs
  mapBase <- function(x) ifelse(x <= nTb, x, nT + (x - nTb))
  mapS    <- function(x) ifelse(x <= nTs, nTb + x, nT + Nb + (x - nTs))
  w <- nT + Nb + Ns + 1L
  eb <- base$edge; es <- S$edge
  ebm <- cbind(mapBase(eb[, 1]), mapBase(eb[, 2]))
  esm <- cbind(mapS(es[, 1]), mapS(es[, 2]))
  rsMapped <- mapS(getRoot(S))
  Dm <- mapBase(D); Am <- mapBase(parentOf(eb, D))
  ebm <- ebm[!(ebm[, 1] == Am & ebm[, 2] == Dm), , drop = FALSE]
  newEdges <- rbind(ebm, esm, c(Am, w), c(w, Dm), c(w, rsMapped))
  tr <- structure(list(edge = newEdges,
                       tip.label = c(base$tip.label, S$tip.label),
                       Nnode = Nb + Ns + 1L), class = "phylo")
  ape::read.tree(text = ape::write.tree(tr))
}

# Oracle: build phyDat and score with TreeSearch::TreeLength (Fitch, EW).
makePhyDat <- function(charMat, nStates) {
  phangorn::phyDat(charMat, type = "USER", levels = as.character(0:(nStates - 1)))
}
treeLen <- function(tree, pd) as.numeric(TreeSearch::TreeLength(tree, pd))

# Brute-force finals for a SINGLE unambiguous character.  tipState: integer
# vector length nTip (single state per tip, 0-indexed).  Returns list:
#   optLen, finalMask (per node bitmask of states appearing in some optimum).
bruteFinals <- function(tree, tipState, nStates) {
  nTip <- length(tree$tip.label)
  internal <- (nTip + 1):(nTip + tree$Nnode)
  edge <- tree$edge
  nInt <- length(internal)
  # enumerate all assignments of states to internal nodes
  grid <- as.matrix(expand.grid(rep(list(0:(nStates - 1)), nInt)))
  fullNode <- integer(nTip + tree$Nnode)
  fullNode[seq_len(nTip)] <- tipState
  lens <- integer(nrow(grid))
  for (g in seq_len(nrow(grid))) {
    fullNode[internal] <- grid[g, ]
    lens[g] <- sum(fullNode[edge[, 1]] != fullNode[edge[, 2]])
  }
  optLen <- min(lens)
  optRows <- which(lens == optLen)
  finalMask <- integer(nTip + tree$Nnode)
  finalMask[seq_len(nTip)] <- bitwShiftL(1L, tipState)
  for (k in seq_along(internal)) {
    nd <- internal[k]
    states <- unique(grid[optRows, k])
    finalMask[nd] <- Reduce(bitwOr, bitwShiftL(1L, states), 0L)
  }
  list(optLen = optLen, finalMask = finalMask)
}

# ============================ EXPERIMENT ============================

disjointCount <- function(X, S) sum(bitwAnd(X, S) == 0L)
scoreSub <- function(tr, cm, nStates) treeLen(tr, makePhyDat(cm[tr$tip.label,,drop=FALSE], nStates))

# single-tip "subtree": root (Nnode=1) with the one tip as its only child
oneTipTree <- function(lbl) structure(list(edge = matrix(c(2L,1L),1,2),
  tip.label = lbl, Nnode = 1L), class = "phylo")

runClip <- function(tr, cm, nStates, s) {
  nTipT <- length(tr$tip.label)
  if (s <= nTipT) {                       # clip a single leaf
    Stips <- tr$tip.label[s]; S <- oneTipTree(Stips)
  } else {
    S <- ape::extract.clade(tr, s); Stips <- S$tip.label
  }
  if (length(Stips) > nTipT - 3) return(NULL)
  base <- ape::drop.tip(tr, Stips, collapse.singles = TRUE)
  if (length(base$tip.label) < 3) return(NULL)

  Smask <- alignedTipMask(S, cm, nStates)
  Sdn <- fitchDown(S, Smask)
  X <- Sdn$prelim[getRoot(S), ]
  Lwithin <- Sdn$length
  Lbase <- scoreSub(base, cm, nStates)

  bMask <- alignedTipMask(base, cm, nStates)
  bDn <- fitchDown(base, bMask)
  E  <- fitchUp(base, bDn$prelim)$final       # edge sets (engine EXACT path)
  Fe <- fitchUpEngine(base, bDn$prelim)$final  # engine SIMPLIFIED final_
  root <- getRoot(base)
  cands <- setdiff(seq_len(length(base$tip.label) + base$Nnode), root)

  rows <- vector("list", length(cands)); k <- 0
  for (D in cands) {
    A <- parentOf(base$edge, D)
    rec <- tryCatch(regraft(base, S, D), error = function(e) NULL)
    if (is.null(rec)) next
    added <- as.integer(round(scoreSub(rec, cm, nStates) - Lbase - Lwithin))
    k <- k + 1
    rows[[k]] <- data.frame(nStates = nStates, clipSize = length(Stips),
      baseTips = length(base$tip.label), added = added,
      pEdge = disjointCount(X, E[D, ]),
      pUnionSimp = disjointCount(X, bitwOr(Fe[A, ], Fe[D, ])),
      pUnionEdge = disjointCount(X, bitwOr(E[A, ], E[D, ])),
      pSingleSimp = disjointCount(X, Fe[D, ]))
  }
  if (!k) return(NULL)
  do.call(rbind, rows[seq_len(k)])
}

randChars <- function(nTip, nChar, nStates, pMissing = 0) {
  m <- matrix(as.character(sample(0:(nStates-1), nTip*nChar, replace = TRUE)), nTip, nChar)
  if (pMissing > 0) { miss <- matrix(runif(nTip*nChar) < pMissing, nTip, nChar); m[miss] <- "?" }
  rownames(m) <- paste0("t", seq_len(nTip)); m
}
cladeSize <- function(tr, node) if (node <= length(tr$tip.label)) 1L else length(ape::extract.clade(tr, node)$tip.label)

runRegime <- function(label, nStates, pMissing, nTrees, seed) {
  set.seed(seed); all <- list()
  for (r in seq_len(nTrees)) {
    nTip <- sample(8:13, 1); nChar <- sample(8:16, 1)
    cm <- randChars(nTip, nChar, nStates, pMissing)
    tr <- ape::rtree(nTip, tip.label = rownames(cm))
    internal <- (nTip+1):(nTip+tr$Nnode); root <- getRoot(tr)
    sizes <- sapply(internal, function(n) cladeSize(tr, n))
    goodInt <- internal[internal != root & sizes >= 2 & sizes <= nTip - 3]
    clips <- c(if (length(goodInt)) sample(goodInt, min(4, length(goodInt))),
               sample(seq_len(nTip), 2))          # + 2 leaf clips
    for (s in clips) {
      res <- tryCatch(runClip(tr, cm, nStates, s), error = function(e) NULL)
      if (!is.null(res)) all[[length(all)+1]] <- res
    }
  }
  df <- do.call(rbind, all); df$regime <- label; df
}

reportRegime <- function(df) {
  cat(sprintf("\n=== %s (nStates=%d, %d candidates, clipSizes %d..%d) ===\n",
      df$regime[1], df$nStates[1], nrow(df), min(df$clipSize), max(df$clipSize)))
  show <- function(tag, pred) {
    e <- df$added - pred
    cat(sprintf("  %-22s err(added-pred): min=%d max=%d | over(<0)=%d under(>0)=%d exact(0)=%d (%.1f%% exact)\n",
        tag, min(e), max(e), sum(e<0), sum(e>0), sum(e==0), 100*mean(e==0)))
  }
  show("P1 edge-set(exact)", df$pEdge)
  show("P2 union-simplified", df$pUnionSimp)   # DEPLOYED fitch_indirect_length
  show("P3 union-edgesets", df$pUnionEdge)
  show("P4 single-simplified", df$pSingleSimp)
}

cat("Running 4 regimes...\n")
res <- list(
  bin    = runRegime("binary-noNA",      2, 0.00, 30, 101),
  binNA  = runRegime("binary-15NA",      2, 0.15, 30, 102),
  mult   = runRegime("multistate4-noNA", 4, 0.00, 30, 103),
  multNA = runRegime("multistate4-20NA", 4, 0.20, 30, 104))
for (df in res) reportRegime(df)
cat("\nSaved results.rds\n")
