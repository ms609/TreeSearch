# Independent TBR-neighbourhood oracle for agent-issues/TreeSearch#147.
#
# Deliberately shares no code with src/rearrange.cpp, R/TBR.R or src/ts_tbr.cpp:
# the tree is held as a bare undirected edge list, splits are recovered by
# deleting an edge and flood-filling one side, and TBR is performed by the
# textbook definition (delete an edge, suppress the two exposed degree-2
# vertices, subdivide one edge in each fragment, join the two new vertices).
#
# Not part of the package; not sourced by any test.

# ---- undirected representation -------------------------------------------

# Convert a `phylo` (rooted or not) to an undirected edge list.
AsUndirected <- function(tree) {
  edge <- tree[["edge"]]
  list(nTip = length(tree[["tip.label"]]),
       from = as.integer(edge[, 1]),
       to = as.integer(edge[, 2]))
}

# Vertices reachable from `start`, optionally refusing to cross edge `blocked`.
Reach <- function(g, start, blocked = 0L) {
  nV <- max(c(g[["from"]], g[["to"]], start))
  seen <- logical(nV)
  seen[[start]] <- TRUE
  keep <- seq_along(g[["from"]]) != blocked
  from <- g[["from"]][keep]
  to <- g[["to"]][keep]
  repeat {
    grew <- FALSE
    hit <- seen[from] & !seen[to]
    if (any(hit)) { seen[to[hit]] <- TRUE; grew <- TRUE }
    hit <- seen[to] & !seen[from]
    if (any(hit)) { seen[from[hit]] <- TRUE; grew <- TRUE }
    if (!grew) break
  }
  which(seen)
}

# Canonical key for the *unrooted* topology of `g`.
#
# Each edge contributes the set of leaves on one side, complemented if needed so
# that leaf 1 is always absent.  Trivial splits are dropped, so a degree-2 root
# and any suppressed vertex are invisible to the key.
SplitKey <- function(g) {
  nTip <- g[["nTip"]]
  full <- seq_len(nTip)
  masks <- character(0)
  for (i in seq_along(g[["from"]])) {
    side <- intersect(Reach(g, start = g[["to"]][[i]], blocked = i), full)
    if (1L %in% side) {
      side <- setdiff(full, side)
    }
    if (length(side) < 2L || length(side) > nTip - 2L) next
    masks <- c(masks, paste0(sort(side), collapse = "."))
  }
  paste0(sort(unique(masks)), collapse = "|")
}

# Delete a degree-2 vertex, joining its two neighbours directly.
SuppressVertex <- function(g, v) {
  inc <- which(g[["from"]] == v | g[["to"]] == v)
  stopifnot(length(inc) == 2L)
  nbr <- c(g[["from"]][inc], g[["to"]][inc])
  nbr <- nbr[nbr != v]
  stopifnot(length(nbr) == 2L)
  g[["from"]] <- c(g[["from"]][-inc], nbr[[1]])
  g[["to"]] <- c(g[["to"]][-inc], nbr[[2]])
  g
}

# Drop every degree-2 internal vertex, e.g. the root of a rooted `phylo`, so
# that enumeration runs over the 2n - 3 edges of the unrooted tree rather than
# the 2n - 2 edges of its rooted representation.
Unroot <- function(g) {
  repeat {
    deg <- tabulate(c(g[["from"]], g[["to"]]))
    twos <- setdiff(which(deg == 2L), seq_len(g[["nTip"]]))
    if (!length(twos)) break
    g <- SuppressVertex(g, twos[[1]])
  }
  g
}

# Subdivide edge `i` of `g` with a new vertex `v`; returns the modified graph.
Subdivide <- function(g, i, v) {
  far <- g[["to"]][[i]]
  g[["to"]][[i]] <- v
  g[["from"]] <- c(g[["from"]], v)
  g[["to"]] <- c(g[["to"]], far)
  g
}

# Every unrooted topology one TBR move from `tree`, as canonical split keys.
# `includeSelf = FALSE` drops the starting topology, which TBR reaches whenever
# it rejoins where it cut.
TbrOracle <- function(tree, includeSelf = FALSE) {
  g0 <- Unroot(AsUndirected(tree))
  self <- SplitKey(g0)
  nTip <- g0[["nTip"]]
  spare <- max(c(g0[["from"]], g0[["to"]])) + 1L

  out <- character(0)
  for (cut in seq_along(g0[["from"]])) {
    u <- g0[["from"]][[cut]]
    v <- g0[["to"]][[cut]]
    # A leaf on each side of the cut, identified before any suppression.
    anchorU <- intersect(Reach(g0, start = u, blocked = cut), seq_len(nTip))[[1]]
    anchorV <- intersect(Reach(g0, start = v, blocked = cut), seq_len(nTip))[[1]]

    g <- g0
    g[["from"]] <- g[["from"]][-cut]
    g[["to"]] <- g[["to"]][-cut]
    for (end in c(u, v)) {
      if (end > nTip) g <- SuppressVertex(g, end)
    }

    compU <- Reach(g, anchorU)
    compV <- Reach(g, anchorV)
    stopifnot(!length(intersect(compU, compV)))

    edgesU <- which(g[["from"]] %in% compU)
    edgesV <- which(g[["from"]] %in% compV)
    # A single-leaf fragment offers no edge to subdivide: attach at the leaf.
    attachU <- if (length(edgesU)) edgesU else 0L
    attachV <- if (length(edgesV)) edgesV else 0L

    for (a in attachU) {
      for (b in attachV) {
        h <- g
        if (a == 0L) {
          newU <- compU
        } else {
          newU <- spare
          h <- Subdivide(h, a, newU)
        }
        if (b == 0L) {
          newV <- compV
        } else {
          newV <- spare + 1L
          h <- Subdivide(h, b, newV)
        }
        h[["from"]] <- c(h[["from"]], newU)
        h[["to"]] <- c(h[["to"]], newV)
        out <- c(out, SplitKey(h))
      }
    }
  }
  out <- unique(out)
  if (!includeSelf) out <- setdiff(out, self)
  sort(out)
}

# Split key of a `phylo`, for comparing package output against the oracle.
TreeKey <- function(tree) SplitKey(AsUndirected(tree))
