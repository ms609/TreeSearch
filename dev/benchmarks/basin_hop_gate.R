#!/usr/bin/env Rscript
# Mission A, route #2 (basin-HOPPING schedule) -- GO/NO-GO GATE (advisor steps 1+2).
#
# Reuses the PROVEN project5432 idioms from reeval/geom_full.R + basin_5432.R +
# recon5432b.R (ordered-char-aware scorer; TNT integer->name relabel; self-checked
# against the meta's 1943). Computes the decisive predictor:
#   TBRDist( TS best-kick floor tree  ->  1943 TNT floor )
# Short bounded uphill (a few TBR) -> accept-worse (drift schedule) can plausibly bridge
# it, probe/build worth it. 15+ TBR -> accept-worse from the floor is a random walk and
# route #2 is likely dead. TS HOLDING its floor already proves it is outside 1943's
# ~8-12 TBR capture radius, so the distance is >~12; this measures HOW MUCH.
# BUILDLESS: reads existing artifacts, re-scores + TBRDist. No search.
#
# Env: NEXUS (matrix), FLOOR (1943 TNT tree, TNT format), REEVAL (dir of TS trees).

.libPaths(c("/nobackup/pjjg18/prlib", "/nobackup/pjjg18/TreeSearch/lib", .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

NEXUS  <- Sys.getenv("NEXUS",  "/nobackup/pjjg18/lgsweep/matrices/project5432.nex")
FLOOR  <- Sys.getenv("FLOOR",  "/nobackup/pjjg18/floors/project5432_tnt_floor_one.tre")
REEVAL <- Sys.getenv("REEVAL", "/nobackup/pjjg18/reeval")
have_tbr <- requireNamespace("TBRDist", quietly = TRUE)

pd <- ReadAsPhyDat(NEXUS); tips <- names(pd); nt <- length(tips)
score <- function(tr) tryCatch(
  as.numeric(TreeLength(Preorder(tr), pd, concavity = Inf, inapplicable = "missing")),
  error = function(e) NA_real_)
ren <- function(t) RenumberTips(Preorder(t), tips)

cat(sprintf("MATRIX %s : %d taxa (ordered-char EW; TS objective)\n", basename(NEXUS), nt))

# ---- 1943 target (TNT format, integer labels -> relabel by matrix taxon order) ----
tnt <- tryCatch(ReadTntTree(FLOOR, tipLabels = tips), error = function(e) NULL)
if (inherits(tnt, "multiPhylo")) tnt <- tnt[[1]]
if (is.null(tnt)) stop("could not read FLOOR: ", FLOOR)
s1943 <- score(tnt)
cat(sprintf("TARGET 1943 floor: score=%.0f  [self-check ==1943: %s]\n",
            s1943, isTRUE(s1943 == 1943)))
if (!isTRUE(s1943 == 1943))
  cat("  [WARN] self-check FAILED -- relabel/ordering off; TBRDist below is untrustworthy.\n")

# ---- candidate TS trees: kick sweep (mo_*), cold seeds (ts5432_*), intensity, sd_* ----
cand_files <- unique(c(
  Sys.glob(file.path(REEVAL, "mo_5432_*.tre")),
  Sys.glob(file.path(REEVAL, "ts5432_seed*.tre")),
  Sys.glob(file.path(REEVAL, "tsint5432_*.tre")),
  Sys.glob(file.path(REEVAL, "sd_5432_*.tre"))))
cand <- lapply(cand_files, function(f) tryCatch(Preorder(ape::read.tree(f)), error = function(e) NULL))
ok <- !vapply(cand, is.null, logical(1))
cand <- cand[ok]; cand_files <- cand_files[ok]
csc <- vapply(cand, score, numeric(1))
# keep only trees with the full taxon set (TBRDist/relabel need identical leaves)
full <- vapply(cand, function(t) setequal(t$tip.label, tips), logical(1))
cat(sprintf("\nloaded %d candidate TS trees (%d with full %d-taxon set)\n",
            length(cand), sum(full), nt))
o <- order(csc)
cat("=== candidate TS-tree scores (sorted; gap vs 1943) ===\n")
for (i in o) cat(sprintf("  %-46s %s (gap +%s)%s\n", basename(cand_files[i]),
    ifelse(is.na(csc[i]), "NA", format(csc[i])),
    ifelse(is.na(csc[i]), "?", format(csc[i] - 1943)),
    ifelse(full[i], "", "  [partial taxa - skipped for TBRDist]")))

# ---- TBRDist(1943 -> best kick tree + cold contrast) ----
if (have_tbr && isTRUE(s1943 == 1943)) {
  usable <- o[full[o] & is.finite(csc[o])]
  best_kick <- usable[seq_len(min(3, length(usable)))]        # the 3 lowest-scoring (~1945-1946)
  cold_cands <- usable[grepl("ts5432_seed", basename(cand_files[usable]))]
  cold <- if (length(cold_cands)) cold_cands[which.max(csc[cold_cands])] else tail(usable, 1)  # a ~1950 cold tree
  # method-diverse near-optimal trees: closes the reweight-schedule loophole -- if every
  # generation method sits ~100 TBR from 1943, no accept-worse schedule (reweight incl.) approaches it.
  pick_lowest <- function(pat) { ix <- usable[grepl(pat, basename(cand_files[usable]))]
    if (length(ix)) ix[which.min(csc[ix])] else integer(0) }
  diverse <- c(pick_lowest("random"), pick_lowest("largefirst"), pick_lowest("invweight"))
  probe <- unique(c(best_kick, cold, diverse))
  tntR <- ren(tnt)
  cat("\n=== TBRDist(1943 -> tree)   [tbr_min .. tbr_max | CID] ===\n")
  rows <- list()
  for (i in probe) {
    tR <- ren(cand[[i]])
    d <- tryCatch(TBRDist::TBRDist(tntR, tR, exact = FALSE),
                  error = function(e) list(tbr_min = NA, tbr_max = NA))
    cc <- tryCatch(as.numeric(TreeDist::ClusteringInfoDist(tntR, tR, normalize = TRUE)),
                   error = function(e) NA_real_)
    cat(sprintf("  %-46s score=%s  TBR=%s..%s  CID=%.3f\n", basename(cand_files[i]),
        format(csc[i]), d$tbr_min, d$tbr_max, cc))
    rows[[length(rows) + 1L]] <- data.frame(tree = basename(cand_files[i]), score = csc[i],
        tbr_min = d$tbr_min, tbr_max = d$tbr_max, cid = cc, stringsAsFactors = FALSE)
  }
  write.csv(do.call(rbind, rows), file.path(REEVAL, "basin_hop_gate_tbr.csv"), row.names = FALSE)
  cat("\nREAD: the lowest-score row is the TS best-kick floor. Its tbr_min is the gate.\n")
} else if (!have_tbr) {
  cat("\n[WARN] TBRDist package not on .libPaths -- cannot compute the predictor.\n")
}
cat("GATE_DONE\n")
