# Shared helpers for the T-364/T-370 three-arm constrained-search battery.
#
# Loaded by t364_probe.R (orientation-rate gate) and t364_cell.R (wall battery).
# Deliberately free of TreeSearch calls so it can be sourced against any arm.
#
# ROOTING VOCABULARY, used throughout.  A constraint split is an *unrooted*
# bipartition.  build_constraint() (ts_constraint.cpp) canonicalizes every mask
# so that dataset tip 0 -- R tip names(dataset)[1] -- is OUTSIDE it, and the
# pre-T-384 map_constraint_nodes() then demanded a node whose descendant tip set
# EQUALS that canonical mask.  So a tree can display the constraint yet still be
# unmappable, if only the tip-0 side is a rooted clade.  Hence three distinct
# quantities per tree, never conflated:
#   clade_canonical  -- the tip-0-EXCLUDING side is a rooted clade  => maps
#   clade_complement -- the tip-0-INCLUDING side is a rooted clade  => displays
#                       the split but mapped to -1 before the T-384 fix
#   compliant        -- either of the above (the honest unrooted test)
#
# All set comparisons are done on tip LABELS via direct edge-matrix descendant
# accumulation.  NOT as.Splits() + %in%: S4 %in% on Splits silently falls
# through to base::%in% when TreeTools is unattached and answers FALSE
# (memory loadall-is-not-rcmdcheck, splits-in-operator-testthat-trap).  The one
# sound figure in the T-384 record came from exactly this direct computation.

MBANK_FIXED_SAMPLE <- c(
  # Small (20-30 taxa)
  "project532", "project2346", "project2451", "project4501",
  "project944", "project971_(1)", "project2762",
  # Medium (31-60 taxa)
  "project826", "project561", "project571", "project4146_(3)",
  "project3688", "project4049", "project423",
  # Large (61-120 taxa)
  "project4286", "project4359", "project4397", "project2084_(1)",
  "project2771", "project2184", "project3938",
  # XLarge (121+ taxa)
  "syab07201", "project4133", "project804", "project4284"
)

find_first_dir <- function(cands) {
  for (d in cands) if (nzchar(d) && dir.exists(d)) {
    return(normalizePath(d, winslash = "/"))
  }
  NA_character_
}

neotrans_dir <- function() {
  d <- find_first_dir(c(
    Sys.getenv("NEOTRANS_DIR", ""),
    file.path("/nobackup", Sys.getenv("USER"), "neotrans", "inst", "matrices"),
    "C:/Users/pjjg18/GitHub/neotrans/inst/matrices"
  ))
  if (is.na(d)) stop("neotrans matrices dir not found; set NEOTRANS_DIR")
  d
}

catalogue_path <- function() {
  p <- Sys.getenv("CAT_CSV", "")
  if (nzchar(p) && file.exists(p)) return(p)
  for (q in c("dev/benchmarks/mbank_catalogue.csv", "mbank_catalogue.csv",
              file.path(Sys.getenv("REPO", "."), "dev/benchmarks/mbank_catalogue.csv"))) {
    if (file.exists(q)) return(q)
  }
  stop("mbank_catalogue.csv not found; set CAT_CSV")
}

load_catalogue <- function() {
  ct <- read.csv(catalogue_path(), stringsAsFactors = FALSE)
  rownames(ct) <- ct$key
  ct
}

# EW Fitch, gaps -> missing: the regime kick_anytime.R uses, kept identical so
# these numbers sit alongside the rest of the corpus work.  Constraint handling
# is orthogonal to the scoring kernel.
load_matrix <- function(key, catalogue) {
  if (!key %in% catalogue$key) stop("key not in catalogue: ", key)
  row <- catalogue[key, ]
  # SEQUESTER the validation split (memory validation-set-sequestered).
  if (!identical(row$split, "training")) {
    stop(sprintf("key %s is split='%s' -- validation is SEQUESTERED", key, row$split))
  }
  f <- file.path(neotrans_dir(), row$filename)
  if (!file.exists(f)) stop("matrix file not found: ", f)
  pd <- suppressWarnings(TreeTools::ReadAsPhyDat(f))
  m <- TreeTools::PhyDatToMatrix(pd, ambigNA = FALSE)
  m[m == "-"] <- "?"
  TreeTools::MatrixToPhyDat(m)
}

# ---------------------------------------------------------------------------
# Descendant tip sets, by direct edge-matrix accumulation.
# ---------------------------------------------------------------------------
# Returns a list over ALL nodes (1..max(edge)) of sorted tip-label vectors, plus
# the root id.  Postorder accumulation: children before parents.
desc_tip_labels <- function(tree) {
  edge <- tree$edge
  n_tip <- length(tree$tip.label)
  n_all <- max(edge)
  parent <- edge[, 1]
  child <- edge[, 2]
  root <- setdiff(parent, child)
  if (length(root) != 1L) stop("tree has ", length(root), " roots; expected 1")

  # Order edges so that every child is processed before its parent: repeatedly
  # peel nodes all of whose children are done.  ape's postorder reorder is the
  # fast path; fall back to an explicit peel if reorder is unavailable.
  tr <- ape::reorder.phylo(tree, "postorder")
  pe <- tr$edge[, 1]
  ce <- tr$edge[, 2]

  sets <- vector("list", n_all)
  for (i in seq_len(n_tip)) sets[[i]] <- tree$tip.label[i]
  for (k in seq_along(ce)) {
    p <- pe[k]
    c_ <- ce[k]
    sets[[p]] <- c(sets[[p]], sets[[c_]])
  }
  sets <- lapply(sets, function(s) if (is.null(s)) character(0) else sort(s))
  list(sets = sets, root = root, n_tip = n_tip)
}

# Edge-count depth from a focal node to every tip, over the tree as an
# UNDIRECTED graph.  O(n) breadth-first; used to pick a topologically distant
# rogue without an O(n^2) distance matrix (the corpus reaches 4062 tips).
tip_depths_from <- function(tree, focal) {
  edge <- tree$edge
  n_all <- max(edge)
  adj_from <- c(edge[, 1], edge[, 2])
  adj_to <- c(edge[, 2], edge[, 1])
  ord <- order(adj_from)
  adj_from <- adj_from[ord]
  adj_to <- adj_to[ord]
  starts <- c(1L, 1L + cumsum(tabulate(adj_from, nbins = n_all)))

  depth <- rep(NA_integer_, n_all)
  depth[focal] <- 0L
  frontier <- focal
  while (length(frontier)) {
    nxt <- integer(0)
    for (v in frontier) {
      lo <- starts[v]
      hi <- starts[v + 1L] - 1L
      if (hi >= lo) {
        nb <- adj_to[lo:hi]
        nb <- nb[is.na(depth[nb])]
        if (length(nb)) {
          depth[nb] <- depth[v] + 1L
          nxt <- c(nxt, nb)
        }
      }
    }
    frontier <- nxt
  }
  n_tip <- length(tree$tip.label)
  stats::setNames(depth[seq_len(n_tip)], tree$tip.label)
}

# ---------------------------------------------------------------------------
# Constraint generation: clade + distant rogue.
# ---------------------------------------------------------------------------
# WHY THIS SHAPE.  An easily-satisfied constraint measures nothing, and a group
# the data already hands you is easily satisfied.  So the group is built to
# CONFLICT with the unconstrained MP tree by construction: take a real clade C of
# that tree and add one tip r that sits maximally far from C, forcing the search
# to drag a distant taxon in.  That is a moderate, realistic conflict -- the kind
# a worker imposes from external evidence -- unlike a maximally-interleaved
# every-nth group, which on a 173-tip matrix forces an enormous penalty and a
# search that may not converge in budget.
#
# want_tip0 STRATIFIES on which SIDE of the bipartition holds dataset tip 0, and
# that is the axis that decides whether arm 2 can express its pathology at all.
# The reasoning, which is worth not re-deriving:
#   * build_constraint() canonicalizes the mask so tip 0 is OUTSIDE it, so the
#     canonical mask is ALWAYS the tip-0-excluding side -- which of the two sides
#     I hand in as the "1" group is immaterial, it gets flipped anyway.
#   * The pre-T-384 mapping needs that canonical side to be a rooted CLADE.  It
#     is not a clade exactly when the tree's root position lies inside it.
#   * So putting tip 0 in the SMALL side makes the canonical mask the LARGE side,
#     which the construction root very probably falls inside => complement-rooted.
#     Putting tip 0 in the large side makes the canonical mask small => the root
#     rarely falls inside => mapping usually succeeds.
# Hence: want_tip0 = TRUE is the stratum that stresses arm 2.
#
# Groups are drawn from tree_u's UNROOTED splits, both sides of each, not from its
# rooted clades: MaximizeParsimony() returns trees rooted at tip 0 (it re-roots so
# the root's first child is a tip, and Preorder() can only put tip 1 there), so no
# non-root clade ever contains tip 0 and a rooted-clade search cannot reach the
# want_tip0 = TRUE stratum at all.
#
# Returns NULL when no split side falls in the requested size band.
gen_constraint <- function(tree_u, tip0, want_tip0, size_lo, size_hi, rng) {
  dt <- desc_tip_labels(tree_u)
  n_tip <- dt$n_tip
  all_tips <- tree_u$tip.label
  non_root <- setdiff(seq_along(dt$sets), dt$root)

  # Candidate = (defining node, chosen side).  Both sides of every split.
  cand <- list()
  for (v in non_root) {
    a <- dt$sets[[v]]
    b <- sort(setdiff(all_tips, a))
    for (side in list(a, b)) {
      if (length(side) >= size_lo && length(side) <= size_hi &&
            (tip0 %in% side) == want_tip0) {
        cand[[length(cand) + 1L]] <- list(v = v, side = side)
      }
    }
  }
  if (!length(cand)) return(NULL)

  # Prefer the largest admissible side: a bigger enforced group is a stronger,
  # less trivially-satisfiable constraint.  Deterministic tie-break by rng.
  sizes <- vapply(cand, function(x) length(x$side), integer(1))
  best <- which(sizes == max(sizes))
  pick <- cand[[if (length(best) == 1L) best else best[rng(length(best))]]]
  mrca <- pick$v
  clade <- pick$side

  depths <- tip_depths_from(tree_u, mrca)
  outside <- setdiff(tree_u$tip.label, clade)
  if (!length(outside)) return(NULL)
  od <- depths[outside]
  rogue <- outside[which.max(od)]
  # A rogue that is tip 0 would silently flip the stratum.
  if (rogue == tip0) {
    od2 <- od[names(od) != tip0]
    if (!length(od2)) return(NULL)
    rogue <- names(od2)[which.max(od2)]
  }
  group <- sort(c(clade, rogue))
  if (length(group) >= n_tip - 1L) return(NULL)

  # Structural non-triviality: the group must NOT already be a clade of the
  # unconstrained tree, on either side (complement-aware).
  comp <- sort(setdiff(tree_u$tip.label, group))
  displayed <- any(vapply(dt$sets, function(s) {
    identical(s, group) || identical(s, comp)
  }, logical(1)))
  if (displayed) return(NULL)

  list(group = group, clade = clade, rogue = rogue,
       has_tip0 = tip0 %in% group, mrca_size = length(clade))
}

# One binary character: "1" = enforced group, "0" = the rest.  .PrepareConstraint
# reads the "1" set as the split to enforce (MaximizeParsimony.R:128-146).
constraint_phydat <- function(group, all_tips) {
  m <- matrix(ifelse(all_tips %in% group, "1", "0"),
              nrow = length(all_tips), dimnames = list(all_tips, "cons"))
  TreeTools::MatrixToPhyDat(m)
}

# ---------------------------------------------------------------------------
# Rooting-explicit compliance test (see vocabulary at the top of this file).
# ---------------------------------------------------------------------------
classify_tree <- function(tree, group, tip0) {
  dt <- desc_tip_labels(tree)
  grp <- sort(group)
  comp <- sort(setdiff(tree$tip.label, grp))
  # Which side is canonical: build_constraint() puts tip 0 OUTSIDE the mask.
  canon <- if (tip0 %in% grp) comp else grp
  other <- if (tip0 %in% grp) grp else comp
  non_root <- setdiff(seq_along(dt$sets), dt$root)
  hit <- function(target) {
    any(vapply(non_root, function(v) identical(dt$sets[[v]], target), logical(1)))
  }
  cc <- hit(canon)
  cx <- hit(other)
  list(clade_canonical = cc, clade_complement = cx, compliant = cc || cx)
}
