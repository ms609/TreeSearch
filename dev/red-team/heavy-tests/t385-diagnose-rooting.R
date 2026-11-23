# T-385 diagnostic: what rooting are the returned trees actually at, and do the
# pool members agree under a genuinely common one?
#
# Motivated by an apparent contradiction in t385-xform-report-agreement.R:
# all 32 returned trees scored 183 AS RETURNED, but scored 183/182 after
# RootTree(t, t$tip.label[1]).  Both were described as "a common rooting", so at
# most one of those descriptions is right.  Settle it before designing the fix.

lib <- if (dir.exists(".agent-t385")) ".agent-t385" else NULL
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = lib)
  library("TreeTools")
})

if (!file.exists("dev/red-team/heavy-tests/t385-make-xform-data.R")) {
  stop("Run this from the package root: the source() path below is relative.")
}
source("dev/red-team/heavy-tests/t385-make-xform-data.R")

dat <- MakeXformData()
ds <- dat$dataset
h <- dat$hierarchy

set.seed(11)
res <- suppressWarnings(
  MaximizeParsimony(ds, inapplicable = "xform", hierarchy = h,
                    maxReplicates = 4L, verbosity = 0L)
)
cat("reported:", unique(attr(res, "score")), "| n trees:", length(res), "\n\n")

Score <- function(tr) TreeLength(tr, ds, inapplicable = "xform", hierarchy = h)

# Is tip 1 (the dataset's first taxon == the kernel's tip 0) a child of the root?
RootedAtTip1 <- function(tr) {
  root <- length(tr$tip.label) + 1L
  tip1 <- match(names(ds)[1], tr$tip.label)
  tip1 %in% tr$edge[tr$edge[, 1] == root, 2]
}

# Does the root have exactly 2 children (rooted) or 3+ (unrooted-as-stored)?
RootDegree <- function(tr) {
  root <- length(tr$tip.label) + 1L
  sum(tr$edge[, 1] == root)
}

info <- data.frame(
  i           = seq_along(res),
  asReturned  = vapply(res, Score, numeric(1)),
  rootDeg     = vapply(res, RootDegree, integer(1)),
  atTip1      = vapply(res, RootedAtTip1, logical(1)),
  viaRootTree = vapply(res, function(tr) Score(RootTree(tr, names(ds)[1])),
                       numeric(1))
)
info$rtRootDeg <- vapply(res, function(tr) RootDegree(RootTree(tr, names(ds)[1])),
                         integer(1))

cat("--- per-tree ---\n")
print(head(info, 12))
cat("...\n\n")

cat("as-returned scores        :", paste(sort(unique(info$asReturned)),
                                        collapse = " "), "\n")
cat("via RootTree(., taxon 1)  :", paste(sort(unique(info$viaRootTree)),
                                        collapse = " "), "\n")
cat("root degree as returned   :", paste(sort(unique(info$rootDeg)),
                                        collapse = " "), "\n")
cat("root degree via RootTree  :", paste(sort(unique(info$rtRootDeg)),
                                        collapse = " "), "\n")
cat("tip1 is a root child      :", sum(info$atTip1), "/", nrow(info), "\n\n")

# The crux: are the trees whose score MOVES under RootTree the same ones that
# were not already rooted at taxon 1?
moved <- info$asReturned != info$viaRootTree
cat("--- trees whose score moves under RootTree ---\n")
cat("n moved:", sum(moved), "\n")
if (any(moved)) {
  cat("of those, already rooted at tip1:", sum(info$atTip1[moved]), "\n")
  cat("root degree of movers (as returned):",
      paste(sort(unique(info$rootDeg[moved])), collapse = " "), "\n")
  cat("deltas:", paste(unique(info$viaRootTree[moved] - info$asReturned[moved]),
                       collapse = " "), "\n")
}

# Are the MOVED trees topologically distinct from the others, or the same
# topology stored differently?  If RootTree changes the score, and the score is
# supposed to depend only on the unrooted topology, then either the topologies
# differ or the scorer is rooting-sensitive (which is the finding).
cat("\n--- topology check ---\n")
# Canonical key: root every tree at the same taxon, sort, serialise.  (Avoids
# as.Splits comparisons -- see the Splits/%in% dispatch trap.)
TopoKey <- function(tr) {
  ape::write.tree(SortTree(RootTree(tr, names(ds)[1])))
}
keys <- vapply(res, TopoKey, character(1))
cat("distinct topologies among returned trees:", length(unique(keys)),
    "of", length(res), "\n")
# Do trees sharing one topology share a score?
byTopo <- tapply(info$viaRootTree, keys, function(x) length(unique(x)))
cat("topologies whose duplicate copies disagree on score:",
    sum(byTopo > 1L), "\n")
cat("score range within the pool at a common rooting:",
    paste(range(info$viaRootTree), collapse = " - "), "\n")
