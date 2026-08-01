# Soft-Sankoff reference implementation (pure R, exploration only).
#
# See dev/plans/2026-08-01-soft-sankoff-temperature-dial.md
#
# Sankoff's DP with `min` replaced by a soft-min at temperature T:
#
#   softmin_T(x) = -T * log(sum_j exp(-x_j / T))
#   S_v(i)       = sum_{c in children(v)} softmin_T_j [ cost_c(i, j) + S_c(j) ]
#
# T -> 0 recovers hard weighted parsimony.  T = 1 with cost = -log(P) recovers
# Felsenstein's pruning algorithm exactly.  The identity is easiest to see by
# writing l_v(i) = exp(-S_v(i) / T) and w_c(i, j) = exp(-cost_c(i, j) / T):
#
#   l_v(i) = prod_c sum_j w_c(i, j) l_c(j)
#
# i.e. soft-Sankoff at temperature T *is* pruning with transition weights
# raised elementwise to the power 1/T.  Everything below works in cost (log)
# space regardless, because w underflows for small T.
#
# Nothing here is exported or optimised; it exists to make the claim falsifiable
# without a compiled build.

# ---------------------------------------------------------------------------
# Primitives
# ---------------------------------------------------------------------------

# Numerically stable soft-minimum.  `temperature = 0` returns the hard minimum.
SoftMin <- function(x, temperature) {
  finite <- is.finite(x)
  if (!any(finite)) {
    return(Inf)
  }
  m <- min(x[finite])
  if (temperature <= 0) {
    return(m)
  }
  # softmin = m - T * log(sum exp(-(x - m) / T)); the shift keeps every
  # exponent in [-Inf, 0], so the sum is in [1, length(x)] and cannot overflow.
  m - temperature * log(sum(exp(-(x[finite] - m) / temperature)))
}

# Postorder sequence of internal nodes (children before parents) for a rooted
# binary tree in `ape` edge-matrix form.  Deliberately does not use
# ape::reorder.phylo, so a bug there cannot mask a bug here.
PostorderNodes <- function(edge, nTip) {
  nNode <- max(edge) - nTip
  kids <- vector("list", nTip + nNode)
  for (e in seq_len(nrow(edge))) {
    parent <- edge[e, 1]
    kids[[parent]] <- c(kids[[parent]], edge[e, 2])
  }
  done <- c(rep(TRUE, nTip), rep(FALSE, nNode))
  order <- integer(0)
  repeat {
    ready <- which(!done & vapply(seq_along(done), function(v) {
      length(kids[[v]]) > 0 && all(done[kids[[v]]])
    }, logical(1)))
    if (!length(ready)) break
    order <- c(order, ready)
    done[ready] <- TRUE
  }
  if (!all(done)) {
    stop("edge matrix does not describe a rooted tree")
  }
  list(order = order, kids = kids, root = setdiff(seq_along(done), edge[, 2]))
}

# Tip cost matrix: 0 for an observed state, Inf otherwise.  `states` is a list
# of integer vectors (allowing polymorphism / ambiguity) of length nTip.
TipCosts <- function(states, nStates) {
  out <- matrix(Inf, length(states), nStates)
  for (i in seq_along(states)) {
    out[i, states[[i]]] <- 0
  }
  out
}

# ---------------------------------------------------------------------------
# Down-pass
# ---------------------------------------------------------------------------

# `cost` is either a single nStates x nStates matrix (cost[from, to]) applied to
# every edge, or a list indexed by *child node* giving a per-branch matrix.
.EdgeCost <- function(cost, child) {
  if (is.list(cost)) cost[[child]] else cost
}

# Returns the nNode x nStates matrix of conditional costs S_v(i), indexed by
# node number (tip rows carry the tip costs).
SoftSankoffDownpass <- function(edge, tipCosts, cost, temperature) {
  nTip <- nrow(tipCosts)
  nStates <- ncol(tipCosts)
  po <- PostorderNodes(edge, nTip)
  S <- matrix(NA_real_, max(edge), nStates)
  S[seq_len(nTip), ] <- tipCosts
  for (v in po[["order"]]) {
    acc <- numeric(nStates)
    for (child in po[["kids"]][[v]]) {
      cc <- .EdgeCost(cost, child)
      acc <- acc + vapply(seq_len(nStates), function(i) {
        SoftMin(cc[i, ] + S[child, ], temperature)
      }, numeric(1))
    }
    S[v, ] <- acc
  }
  attr(S, "root") <- po[["root"]]
  attr(S, "kids") <- po[["kids"]]
  attr(S, "order") <- po[["order"]]
  S
}

# Total score.  `rootCost` is an additive per-state cost at the root; supply
# -log(pi) to reproduce a likelihood with root frequencies pi, or leave at 0
# for the parsimony convention.
SoftSankoffScore <- function(edge, tipCosts, cost, temperature,
                             rootCost = 0) {
  S <- SoftSankoffDownpass(edge, tipCosts, cost, temperature)
  root <- attr(S, "root")
  SoftMin(S[root, ] + rootCost, temperature)
}

# ---------------------------------------------------------------------------
# Up-pass: marginal state probabilities at T > 0
# ---------------------------------------------------------------------------

# Returns an nNode x nStates matrix of marginal posterior probabilities of each
# state at each node.  At T -> 0 this degenerates toward an indicator of the MPR
# set, which is why it is only meaningful for T > 0.
SoftSankoffMarginals <- function(edge, tipCosts, cost, temperature,
                                 rootCost = 0) {
  stopifnot(temperature > 0)
  S <- SoftSankoffDownpass(edge, tipCosts, cost, temperature)
  nStates <- ncol(tipCosts)
  root <- attr(S, "root")
  kids <- attr(S, "kids")
  O <- matrix(NA_real_, nrow(S), nStates)
  O[root, ] <- rootCost
  # Preorder: reverse of the postorder over internal nodes.
  for (v in rev(attr(S, "order"))) {
    for (child in kids[[v]]) {
      siblings <- setdiff(kids[[v]], child)
      # Outside cost of `v`, plus every sibling subtree's contribution.
      A <- O[v, ]
      for (s in siblings) {
        sc <- .EdgeCost(cost, s)
        A <- A + vapply(seq_len(nStates), function(i) {
          SoftMin(sc[i, ] + S[s, ], temperature)
        }, numeric(1))
      }
      cc <- .EdgeCost(cost, child)
      O[child, ] <- vapply(seq_len(nStates), function(j) {
        SoftMin(A + cc[, j], temperature)
      }, numeric(1))
    }
  }
  total <- S + O
  t(apply(total, 1, function(row) {
    if (!any(is.finite(row))) return(rep(NA_real_, length(row)))
    w <- exp(-(row - min(row[is.finite(row)])) / temperature)
    w[!is.finite(row)] <- 0
    w / sum(w)
  }))
}

# ---------------------------------------------------------------------------
# Independent oracles
# ---------------------------------------------------------------------------

# Hard Sankoff, written with explicit min-loops rather than SoftMin(), so a bug
# in SoftMin cannot hide inside its own T -> 0 limit.
HardSankoffScore <- function(edge, tipCosts, cost, rootCost = 0) {
  nTip <- nrow(tipCosts)
  nStates <- ncol(tipCosts)
  po <- PostorderNodes(edge, nTip)
  S <- matrix(NA_real_, max(edge), nStates)
  S[seq_len(nTip), ] <- tipCosts
  for (v in po[["order"]]) {
    acc <- numeric(nStates)
    for (child in po[["kids"]][[v]]) {
      cc <- .EdgeCost(cost, child)
      best <- numeric(nStates)
      for (i in seq_len(nStates)) {
        best[i] <- min(cc[i, ] + S[child, ])
      }
      acc <- acc + best
    }
    S[v, ] <- acc
  }
  min(S[po[["root"]], ] + rootCost)
}

# Felsenstein's pruning algorithm, in probability space, written independently
# of everything above.  `P` is a single nStates x nStates transition matrix or a
# list indexed by child node.  Returns the log likelihood.
FelsensteinLogLik <- function(edge, tipLik, P, rootFreq) {
  nTip <- nrow(tipLik)
  nStates <- ncol(tipLik)
  po <- PostorderNodes(edge, nTip)
  L <- matrix(NA_real_, max(edge), nStates)
  L[seq_len(nTip), ] <- tipLik
  for (v in po[["order"]]) {
    acc <- rep(1, nStates)
    for (child in po[["kids"]][[v]]) {
      pp <- if (is.list(P)) P[[child]] else P
      acc <- acc * as.vector(pp %*% L[child, ])
    }
    L[v, ] <- acc
  }
  log(sum(rootFreq * L[po[["root"]], ]))
}

# Fitch parsimony (equal costs, unordered), independent of the Sankoff path.
FitchScore <- function(edge, states, nStates) {
  nTip <- length(states)
  po <- PostorderNodes(edge, nTip)
  sets <- vector("list", max(edge))
  sets[seq_len(nTip)] <- states
  steps <- 0
  for (v in po[["order"]]) {
    kids <- po[["kids"]][[v]]
    inter <- Reduce(intersect, sets[kids])
    if (length(inter)) {
      sets[[v]] <- inter
    } else {
      sets[[v]] <- Reduce(union, sets[kids])
      steps <- steps + length(kids) - 1
    }
  }
  steps
}

# ---------------------------------------------------------------------------
# Convenience: an Mk / Jukes-Cantor transition matrix
# ---------------------------------------------------------------------------

MkTransition <- function(nStates, branchLength) {
  # Mk: equal rates, equal frequencies.  P = exp(Q t) has the closed form below.
  same <- 1 / nStates + (1 - 1 / nStates) * exp(-branchLength * nStates /
                                                  (nStates - 1))
  diff <- (1 - same) / (nStates - 1)
  out <- matrix(diff, nStates, nStates)
  diag(out) <- same
  out
}
