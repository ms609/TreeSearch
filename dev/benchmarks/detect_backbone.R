#!/usr/bin/env Rscript
# =====================================================================
# Mission A (cold-start basin capture) -- Phase 1: backbone DETECTABILITY cut
# =====================================================================
# THE QUESTION (advisor-framed, 2026-07-16):
#   restore238 already proved SUFFICIENCY -- get the deep backbone right and
#   TS's own TBR descends into the basin. So a cold-start constructor reduces
#   to a backbone DETECTOR: can the correct deep split be identified from a
#   COMPUTABLE, DATA-DERIVED signal WITHOUT knowing the answer?
#
#   Concretely: do the TRUE deep-backbone skeleton splits (present in the
#   best-known tree, ABSENT from every near-optimal TS tree -- incl. the
#   238-tip split) score HIGHER under a data-derived per-split signal than the
#   WRONG deep splits (present in TS trees, absent from the best-known tree)?
#
# DECISION RULE:
#   * If some signal ranks true-missing-deep ABOVE wrong-deep (AUC >> 0.5,
#     and the 238-split sits high among comparable-size splits)
#         -> a constructor engine EXISTS: impose top-signal deep splits as a
#            constraint skeleton, RAS-complete, let TS descend. Proceed to
#            Phase 2 (kick-scaffold + constrained-start descent probes).
#   * If no signal separates them (AUC ~ 0.5 or worse; 238-split not
#     distinguished) -> no cheap detector can build what it cannot see
#         -> PIVOT to the reweighting-kick SCHEDULE (angle #2).
#
# WHY DATA-DERIVED ONLY:
#   Search-derived signals are already REFUTED. transient-autoconstraint-5432:
#   a consensus of converged searches encodes the WRONG 1950-basin backbone
#   ("on 5432 the pool converges to 1950 within a few reps, so ANY pool
#   consensus is a 1950 consensus"). QuartetConcordance / PhylogeneticConcordance
#   are computed from the CHARACTER MATRIX, independent of any trajectory.
#
# SIGNALS (all exported, tested, data-derived; R/Concordance.R):
#   * QuartetConcordance(tree, dataset, return="edge")
#       proportion of DECISIVE quartets concordant with each split
#       (Minh2020 sCF analogue; built to sidestep Goloboff2024 criticisms;
#        ambiguous/inapplicable tokens treated as "?").
#   * PhylogeneticConcordance(tree, dataset)
#       proportion of INFORMATIVE characters compatible with each split
#       (the character-compatibility signal).
#   Split size (min side) is recorded as a covariate -- compare within size bins
#   because concordance can correlate with clade size.
#
# USAGE (run on Hamilton where data + current TreeSearch live):
#   Rscript detect_backbone.R <matrix> <best_tree> <ts_trees> [out_prefix]
#     <matrix>     project5432 matrix (.nex -> ReadAsPhyDat, or .tnt)
#     <best_tree>  best-known tree (1943 floor now; 1939 once regenerated)
#     <ts_trees>   multi-Newick file OR directory of near-optimal TS .tre trees
#     [out_prefix] output prefix (default: detect_backbone)
#
# NOTE: concordance is a GROUPING measure -> it is immune to the 5432
# ordering/polymorphism SCORING subtleties (diversity-generation-gates); it
# only reads which taxa share tokens, treating ambiguity as "?".
# =====================================================================

suppressWarnings(suppressMessages({
  library(TreeSearch)
  library(TreeTools)
  library(ape)
}))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: detect_backbone.R <matrix> <best_tree> <ts_trees> [out_prefix]")
}
matPath  <- args[[1]]
bestPath <- args[[2]]
tsPath   <- args[[3]]
outPref  <- if (length(args) >= 4L) args[[4]] else "detect_backbone"

msg <- function(...) cat(sprintf(...), "\n", sep = "")

## -- 1. Load matrix ---------------------------------------------------------
## project5432 regime = EW Fitch, gaps=missing, inapplicable=missing.
## For concordance, "-"/ambiguous are treated as "?" internally, so the signal
## is regime-robust; we still read faithfully.
readMatrix <- function(p) {
  if (grepl("\\.tnt$", p, ignore.case = TRUE)) {
    ReadTntCharacters(p)
  } else {
    ReadAsPhyDat(p)
  }
}
dataset <- readMatrix(matPath)
msg("Loaded matrix: %d taxa x %d characters", length(dataset),
    attr(dataset, "nr"))

## -- 2. Load trees ----------------------------------------------------------
## best tree may be TNT format (`tread` header, space-separated, 0-based numeric
## tips) -> relabel indices to dataset tip names (nexus/m01.tnt share order).
## TS trees are taxon-labelled Newick (R MaximizeParsimony output).
readBest <- function(p, tipNames) {
  first <- tryCatch(readLines(p, n = 1L, warn = FALSE), error = function(e) "")
  isTnt <- grepl("tread", first, ignore.case = TRUE) ||
    grepl("^\\(\\s*[0-9]", first)
  if (isTnt) {
    t <- tryCatch(ReadTntTree(p, tipLabels = tipNames), error = function(e) NULL)
    if (is.null(t)) t <- ReadTntTree(p)          # fall back, relabel manually
    if (!inherits(t, "phylo")) t <- t[[1]]
    if (all(grepl("^[0-9]+$", t$tip.label))) {   # 0-based TNT indices
      t$tip.label <- tipNames[as.integer(t$tip.label) + 1L]
    }
    t
  } else {
    t <- ape::read.tree(p); if (!inherits(t, "phylo")) t <- t[[1]]; t
  }
}
## tsPath: a glob ("...seed*.tre"), a directory, or a single file.
readNewickSet <- function(p) {
  files <- if (grepl("[*?]", p)) Sys.glob(p) else
    if (dir.exists(p)) list.files(p, pattern = "\\.(tre|nwk|tree|txt)$",
                                  full.names = TRUE, ignore.case = TRUE) else p
  out <- list()
  for (f in files) {
    t <- tryCatch(ape::read.tree(f), error = function(e) NULL)
    if (is.null(t)) next
    if (inherits(t, "phylo")) out[[length(out) + 1L]] <- t else
      for (tt in t) out[[length(out) + 1L]] <- tt
  }
  out
}
bestTree <- readBest(bestPath, names(dataset))
tsTrees  <- readNewickSet(tsPath)
stopifnot(length(tsTrees) >= 1L)
msg("Loaded best-known tree (%d tips) + %d near-optimal TS trees",
    length(bestTree$tip.label), length(tsTrees))

## Restrict everything to the common tip set (defensive).
commonTips <- Reduce(intersect,
                     c(list(names(dataset), bestTree$tip.label),
                       lapply(tsTrees, function(t) t$tip.label)))
stopifnot(length(commonTips) >= 4L)
dataset  <- dataset[commonTips, drop = FALSE]
tipOrder <- names(dataset)                     # canonical tip order
bestTree <- keep.tip(bestTree, commonTips)
tsTrees  <- lapply(tsTrees, keep.tip, commonTips)
nTip <- length(tipOrder)
msg("Common tip set: %d taxa", nTip)

## -- 3. Split extraction + canonical keys ----------------------------------
## Membership over the FIXED dataset tip order; canonicalize so tip 1 is FALSE.
splitLogivec <- function(tree) {
  sp <- as.Splits(tree, dataset)              # aligned to `dataset` tip order
  lapply(seq_along(sp), function(i) as.logical(sp[[i]]))
}
canonKey <- function(v) {
  if (isTRUE(v[[1]])) v <- !v
  paste(which(v), collapse = ",")
}
cladeSize <- function(v) min(sum(v), sum(!v))

bestVecs <- splitLogivec(bestTree)
bestKeys <- vapply(bestVecs, canonKey, character(1))
bestSize <- vapply(bestVecs, cladeSize, integer(1))

tsKeySet <- unique(unlist(lapply(tsTrees, function(t)
  vapply(splitLogivec(t), canonKey, character(1)))))

## -- 4. Per-split data-derived signals -------------------------------------
## QuartetConcordance / PhylogeneticConcordance return per-EDGE vectors aligned
## to as.Splits(tree, dataset) order, i.e. to bestVecs / each ts tree's splits.
msg("Computing QuartetConcordance (best tree)... [may take minutes at 482t]")
qcBest <- as.numeric(QuartetConcordance(bestTree, dataset, return = "edge"))
msg("Computing PhylogeneticConcordance (best tree)...")
pcBest <- as.numeric(PhylogeneticConcordance(bestTree, dataset))

## For the WRONG deep splits we need signals on TS-tree splits. Compute per TS
## tree and collect the splits ABSENT from the best-known tree.
bestKeySet <- unique(bestKeys)
tsRows <- list()
for (k in seq_along(tsTrees)) {
  t <- tsTrees[[k]]
  vecs <- splitLogivec(t)
  keys <- vapply(vecs, canonKey, character(1))
  qc <- as.numeric(QuartetConcordance(t, dataset, return = "edge"))
  pc <- as.numeric(PhylogeneticConcordance(t, dataset))
  sz <- vapply(vecs, cladeSize, integer(1))
  absent <- !(keys %in% bestKeySet)            # present in TS, absent from best
  if (any(absent)) {
    tsRows[[length(tsRows) + 1L]] <- data.frame(
      tree = paste0("ts", k), key = keys[absent],
      size = sz[absent], qc = qc[absent], pc = pc[absent],
      class = "wrong", stringsAsFactors = FALSE)
  }
}
tsDf <- if (length(tsRows)) do.call(rbind, tsRows) else
  data.frame(tree = character(), key = character(), size = integer(),
             qc = numeric(), pc = numeric(), class = character())
## dedupe wrong splits (a wrong split can recur across TS trees); keep max signal
if (nrow(tsDf)) {
  tsDf <- do.call(rbind, by(tsDf, tsDf$key, function(d) d[which.max(d$qc), ]))
}

## Best-tree splits: classify true-missing (absent from ALL TS trees) vs shared.
bestDf <- data.frame(
  tree = "best", key = bestKeys, size = bestSize, qc = qcBest, pc = pcBest,
  class = ifelse(bestKeys %in% tsKeySet, "shared", "true_missing"),
  stringsAsFactors = FALSE)

allDf <- rbind(bestDf, tsDf)
write.csv(allDf, paste0(outPref, "_splits.csv"), row.names = FALSE)

## -- 5. The cut -------------------------------------------------------------
## Compare true_missing vs wrong, focusing on DEEP splits. Report overall and
## within size bands; AUC = P(signal(true) > signal(wrong)).
aucStat <- function(pos, neg) {
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}
report <- function(sig, label) {
  tm <- allDf[allDf$class == "true_missing", ]
  wr <- allDf[allDf$class == "wrong", ]
  msg("\n=== SIGNAL: %s ===", label)
  msg("  true_missing splits: %d   wrong splits: %d", nrow(tm), nrow(wr))
  bands <- list(all = 0L, deep30 = 30L, deep60 = 60L, deep120 = 120L)
  for (bn in names(bands)) {
    thr <- bands[[bn]]
    tmv <- tm[tm$size >= thr, sig]; wrv <- wr[wr$size >= thr, sig]
    msg("  [size>=%3d] true n=%2d med=%.3f | wrong n=%2d med=%.3f | AUC(true>wrong)=%.3f",
        thr, length(tmv), stats::median(tmv, na.rm = TRUE),
        length(wrv), stats::median(wrv, na.rm = TRUE), aucStat(tmv, wrv))
  }
}
report("qc", "QuartetConcordance")
report("pc", "PhylogeneticConcordance")

## -- 6. The headline: the deepest true split (the ~238-split) --------------
tm <- allDf[allDf$class == "true_missing", ]
if (nrow(tm)) {
  deepest <- tm[which.max(tm$size), ]
  msg("\n=== DEEPEST TRUE-MISSING SPLIT (the backbone split TS never generates) ===")
  msg("  size=%d  QC=%.3f  PC=%.3f", deepest$size, deepest$qc, deepest$pc)
  ## its percentile among ALL deep splits (true+wrong) of comparable size
  comp <- allDf[allDf$class %in% c("true_missing", "wrong") &
                  allDf$size >= 0.7 * deepest$size, ]
  for (sig in c("qc", "pc")) {
    pct <- mean(comp[[sig]] <= deepest[[sig]], na.rm = TRUE)
    msg("  %s percentile among %d comparable-size deep splits: %.1f%%",
        toupper(sig), nrow(comp), 100 * pct)
  }
}

msg("\nWrote: %s_splits.csv", outPref)
msg("INTERPRET: AUC >> 0.5 on deep bands AND deepest-split high percentile ->",
    " detector exists (Phase 2). AUC ~ 0.5 -> pivot to kick-schedule (angle #2).")
