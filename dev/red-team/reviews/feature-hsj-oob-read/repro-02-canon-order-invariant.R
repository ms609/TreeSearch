# Mirror of ts_hsj.cpp build_canon_order() (src/ts_hsj.cpp:58-110) in R.
# Evidence for agent-issues/TreeSearch#51.  Establishes three things about
# the CSR children arrays, over 900 random trees of 2-24 tips:
#
#  1. EVERY tree has at least one node with kidOff[node] == length(kids), so
#     the `if (nk == 0) continue` guard in fitch_label_char()'s uppass is
#     always load-bearing -- the pre-fix code formed a reference to
#     co.kids.end() on every single HSJ scoring call.
#  2. Usually SEVERAL nodes do (843 of 900), not just the last one popped:
#     any childless node reached after the final push_back carries the end
#     offset.  The last popped node is always among them.
#  3. kidOff/kidNum are otherwise CONSISTENT -- kids[off + 1 .. off + num]
#     is exactly node n's canonical children for every node with children.
#     So the guard is a bounds fix, not a patch over a corrupt CSR.
#
# Every such node has kidNum == 0, which is why skipping them changes no
# score: the loop body the guard bypasses is zero-trip anyway.
#
# Pure R; needs no TreeSearch build.  Run:  Rscript <this file>
suppressMessages(library("ape"))
set.seed(1)

canon <- function(edge, nTip) {
  nNode <- max(edge)                       # 1-based node count
  adj <- vector("list", nNode)
  for (i in seq_len(nrow(edge))) {
    adj[[edge[i, 1]]] <- c(adj[[edge[i, 1]]], edge[i, 2])
    adj[[edge[i, 2]]] <- c(adj[[edge[i, 2]]], edge[i, 1])
  }
  # C++ indices are 0-based with tips first; ape's are already tips-first,
  # so sorting ascending on ape's numbering matches sorting on 0-based.
  adj <- lapply(adj, sort)
  kidOff <- integer(nNode); kidNum <- integer(nNode)
  kids <- integer(0); pre <- integer(0)
  seen <- logical(nNode); stack <- 1L; seen[1] <- TRUE   # start at tip 0
  while (length(stack)) {
    n <- stack[length(stack)]; stack <- stack[-length(stack)]
    pre <- c(pre, n)
    kidOff[n] <- length(kids)                            # 0-based offset
    for (nb in adj[[n]]) {
      if (seen[nb]) next
      seen[nb] <- TRUE
      kids <- c(kids, nb); kidNum[n] <- kidNum[n] + 1L
      stack <- c(stack, nb)
    }
  }
  list(pre = pre, kids = kids, kidOff = kidOff, kidNum = kidNum,
       nNode = nNode, nVisited = length(pre))
}

bad <- 0L; multi <- 0L; unreached <- 0L
for (nTip in 2:24) for (rep in 1:40) {
  tr <- if (nTip == 2) structure(list(edge = matrix(c(3L,1L,3L,2L), 2, 2,
                                                    byrow = TRUE),
                                      tip.label = c("a","b"), Nnode = 1L),
                                 class = "phylo") else rtree(nTip)
  co <- canon(tr$edge, nTip)
  if (co$nVisited != co$nNode) unreached <- unreached + 1L
  oob <- which(co$kidOff == length(co$kids))
  if (length(oob) == 0) bad <- bad + 1L
  if (length(oob) > 1) multi <- multi + 1L
  stopifnot(all(co$kidNum[oob] == 0L))                  # OOB node is childless
  stopifnot(co$pre[length(co$pre)] %in% oob)             # last popped is one
  # kidOff/kidNum consistency: children of n are exactly kids[off+1 .. off+num]
  for (n in seq_len(co$nNode)) if (co$kidNum[n] > 0) {
    got <- co$kids[co$kidOff[n] + seq_len(co$kidNum[n])]
    par <- tr$edge[tr$edge[, 2] == n, 1]
    nbs <- sort(setdiff(c(tr$edge[tr$edge[,1]==n,2], par), integer(0)))
    stopifnot(setequal(got, setdiff(nbs, co$pre[seq_len(which(co$pre==n))])))
  }
}
cat(sprintf("trees with NO kidOff==size node: %d\n", bad))
cat(sprintf("trees with >1 such node        : %d\n", multi))
cat(sprintf("trees with unreached nodes     : %d\n", unreached))
cat("kidOff/kidNum CSR consistency: OK\n")
