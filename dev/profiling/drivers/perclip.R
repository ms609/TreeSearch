# Per-clip economy: localize the stable ~1.8x TS-vs-TNT throughput factor.
#
# framing.R showed efficiency (candidates) is at parity-or-better, so the
# residual wall gap is THROUGHPUT (~1.3-1.84x same-machine).  Two possible
# sources (advisor's b-vs-c):
#   (b) per-candidate scan cost  -> ns/candidate
#   (c) per-clip overhead        -> amortized over candidates-per-clip
# The recently-added compute_insertion_edge_sets (Wagner fix) does O(n*chars)
# + 3 heap allocs PER CLIP.  If candidates-per-clip is large, that overhead is
# well amortized (not the lever); if small, it matters.  This measures it on
# the SHIPPING kernel via ts_tbr_diagnostics (one TBR-to-convergence per start).
#
# Run single-threaded with nothing else competing (timing-sensitive).

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-p0"),
                                              winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")

fitch_convert <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m)
}
bundle <- function(phy) list(
  contrast = attr(phy, "contrast"),
  tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
  weight   = attr(phy, "weight"), levels = attr(phy, "levels"),
  nTip = length(phy), phy = phy)

dsNames <- strsplit(trimws(Sys.getenv("TS_DATASETS",
            "Wortley2006 Giles2015 Zhu2013")), "\\s+")[[1]]
seeds <- 1:3

rows <- list()
for (nm in dsNames) {
  phy <- fitch_convert(inapplicable.phyData[[nm]]); d <- bundle(phy)
  for (sd in seeds) {
    set.seed(sd); start <- RandomTree(phy, root = TRUE)
    edge <- Preorder(RenumberTips(start, names(phy)))[["edge"]]
    set.seed(sd)
    tt <- system.time(res <- TreeSearch:::ts_tbr_diagnostics(
      edge, d$contrast, d$tip_data, d$weight, d$levels,
      maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, unrooted = TRUE))
    p <- res$passes
    totClips <- sum(p$n_clips_tried)
    totCand  <- sum(p$n_candidates_evaluated)
    wall     <- as.double(tt["elapsed"])
    rows[[length(rows) + 1]] <- data.frame(
      dataset = nm, tips = d$nTip, seed = sd, score = res$score,
      n_passes = nrow(p), clips = totClips, cand = totCand,
      cand_per_clip = round(totCand / totClips, 1),
      ns_per_cand = round(wall / totCand * 1e9, 1),
      us_per_clip = round(wall / totClips * 1e6, 1),
      wall = round(wall, 2), stringsAsFactors = FALSE)
  }
}
res <- do.call(rbind, rows)
agg <- do.call(rbind, lapply(split(res, res$dataset), function(x) data.frame(
  dataset = x$dataset[1], tips = x$tips[1],
  cand_per_clip = round(median(x$cand_per_clip), 1),
  ns_per_cand = round(median(x$ns_per_cand), 1),
  us_per_clip = round(median(x$us_per_clip), 1),
  clips = round(median(x$clips)), cand_M = round(median(x$cand) / 1e6, 1),
  stringsAsFactors = FALSE)))
agg <- agg[order(agg$tips), ]

cat("\n=== Per-clip economy (shipping kernel, RandomTree start, TBR to convergence) ===\n")
cat("cand_per_clip: candidates examined per clip (amortization of per-clip overhead).\n")
cat("ns_per_cand:   wall per candidate (incl. amortized per-clip overhead).\n")
cat("us_per_clip:   wall per clip (per-clip overhead + its candidate scans).\n\n")
print(agg, row.names = FALSE)
write.csv(res, "dev/profiling/perclip_latest.csv", row.names = FALSE)
