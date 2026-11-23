# T-385 (P1) — does MaximizeParsimony's reported XFORM score agree with
# TreeLength() on the tree it just returned?
#
# The finding (2026-07-31, found while building the T-378 benchmark panel):
# attr(res, "score") disagreed with TreeLength(res[[1]], ...) by 10 steps on a
# 36-tip synthetic XFORM dataset, and rerooting the *same* returned topology
# gave a third value again (386 reported / 396 as-returned / 394 rerooted).
#
# This script reproduces that against whatever tip it is run on, BEFORE any fix
# is written, because the tip has moved since the finding was recorded.  It
# reports rather than asserts: it prints a verdict line and exits 0 when the
# three quantities agree, 1 when they do not.
#
# Run:  Rscript dev/red-team/heavy-tests/t385-xform-report-agreement.R
# (expects the package in .agent-t385, i.e. the AGENTS.md tarball build)

lib <- if (dir.exists(".agent-t385")) ".agent-t385" else NULL
suppressPackageStartupMessages({
  library("TreeSearch", lib.loc = lib)
  library("TreeTools")
})

cat("TreeSearch:", as.character(packageVersion("TreeSearch")), "\n")

# ---- A synthetic XFORM dataset -----------------------------------------------
# Generator shared with t385-diagnose-rooting.R so the two cannot drift apart.
if (!file.exists("dev/red-team/heavy-tests/t385-make-xform-data.R")) {
  stop("Run this from the package root: the source() path below is relative.")
}
source("dev/red-team/heavy-tests/t385-make-xform-data.R")

dat <- MakeXformData()
ds <- dat$dataset
h <- dat$hierarchy
cat(sprintf("dataset: %d tips, %d chars, %d blocks x %d secondaries\n",
            length(ds), attr(ds, "nr"), dat$nBlock, dat$nSec))
cat(sprintf("bound on rooting overstatement (sum nSec per block): %d\n",
            dat$nBlock * dat$nSec))

# ---- Search ------------------------------------------------------------------
set.seed(11)
res <- suppressWarnings(
  MaximizeParsimony(ds, inapplicable = "xform", hierarchy = h,
                    maxReplicates = 4L, verbosity = 0L)
)

reported <- attr(res, "score")
cat("\n--- reported ---\n")
cat("attr(res, 'score'):", paste(unique(reported), collapse = " "), "\n")
cat("n trees returned:", length(res), "\n")

# ---- Symptom 1: reported score vs TreeLength() of the returned tree ---------
asReturned <- vapply(seq_along(res), function(i) {
  TreeLength(res[[i]], ds, inapplicable = "xform", hierarchy = h)
}, numeric(1))

cat("\n--- symptom 1: reported vs TreeLength(returned tree) ---\n")
cat("TreeLength as-returned:", paste(round(asReturned, 4), collapse = " "), "\n")
gap1 <- max(abs(asReturned - min(reported)))
cat(sprintf("max |TreeLength - reported| = %g\n", gap1))

# ---- Symptom 2: rerooting the SAME topology changes the score ---------------
# The engine's rooting is an artefact (RandomTree/ts_collapse_pool tip-0), so a
# score that moves under rerooting is the un-pinned reporting gap.
tree1 <- res[[1]]
rootings <- vapply(tree1$tip.label[seq_len(min(8L, length(tree1$tip.label)))],
                   function(tip) {
                     TreeLength(RootTree(tree1, tip), ds,
                                inapplicable = "xform", hierarchy = h)
                   }, numeric(1))

cat("\n--- symptom 2: same topology, different rootings ---\n")
print(round(rootings, 4))
spread <- diff(range(rootings))
cat(sprintf("spread across %d rootings = %g (documented bound = %d)\n",
            length(rootings), spread, dat$nBlock * dat$nSec))

# ---- Symptom 3: is the returned MPT set internally consistent? -------------
cat("\n--- symptom 3: MPT-set self-consistency ---\n")
if (length(res) > 1L) {
  # Root every tree at the SAME NAMED taxon.  Do NOT use res[[i]]$tip.label[1]:
  # after Renumber() that is a different taxon for different trees, so it scores
  # each tree at a different rooting and manufactures a spurious disagreement.
  # (That error produced a retracted "2 distinct scores" claim on 2026-07-31.)
  rootTaxon <- names(ds)[1]
  canon <- vapply(seq_along(res), function(i) {
    TreeLength(RootTree(res[[i]], rootTaxon), ds,
               inapplicable = "xform", hierarchy = h)
  }, numeric(1))
  cat("scores at a common (tip-1) rooting:",
      paste(round(canon, 4), collapse = " "), "\n")
  nDistinct <- length(unique(round(canon, 8)))
  cat(sprintf("distinct scores among %d 'equally parsimonious' trees: %d\n",
              length(canon), nDistinct))
} else {
  canon <- asReturned
  nDistinct <- 1L
  cat("single tree returned; not informative\n")
}

# ---- Verdict -----------------------------------------------------------------
tol <- 1e-8
agree1 <- gap1 < tol
agree2 <- spread < tol
agree3 <- nDistinct == 1L

cat("\n===== VERDICT =====\n")
cat(sprintf("reported == TreeLength(returned) : %s\n",
            if (agree1) "AGREE" else sprintf("DISAGREE by %g", gap1)))
cat(sprintf("score invariant under rerooting  : %s\n",
            if (agree2) "INVARIANT" else sprintf("VARIES, spread %g", spread)))
cat(sprintf("MPT set shares one score         : %s\n",
            if (agree3) "CONSISTENT"
            else sprintf("INCONSISTENT, %d distinct", nDistinct)))

if (agree1 && agree2 && agree3) {
  cat("\nT-385 does NOT reproduce on this tip.\n")
  quit(status = 0)
} else {
  cat("\nT-385 REPRODUCES on this tip.\n")
  quit(status = 1)
}
