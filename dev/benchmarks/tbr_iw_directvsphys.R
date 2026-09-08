# tbr_iw_directvsphys.R -- for ONE failing (nTip, idx), compare the direct
# unrooted path vs the physical-reroot path (TS_PHYS_REROOT=1, exact multi-
# rooting) on the SAME start.  If physical reaches a lower IW than direct, the
# residual is a DIRECT-PATH gap (single-rooting enumeration or incremental
# scoring); if both miss the oracle's improver, it is a deeper/plateau issue.
#
# Usage: Rscript dev/benchmarks/tbr_iw_directvsphys.R [nTip] [idx] [concavity]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-tbr"),
            winslash = "/"))
  library(TreeTools)
})
args <- commandArgs(trailingOnly = TRUE)
nTip <- if (length(args) >= 1) as.integer(args[[1]]) else 16L
idx  <- if (length(args) >= 2) as.integer(args[[2]]) else 19L
conc <- if (length(args) >= 3) as.numeric(args[[3]]) else 10
nChar <- 60L; nState <- 3L

randomData <- function(seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  m <- matrix(sample(0:(nState - 1L), nTip * nChar, replace = TRUE),
              nrow = nTip, dimnames = list(tips, NULL))
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, nTip = length(phy), labels = names(phy))
}
scoreTree <- function(tree, d) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  TreeSearch:::ts_fitch_score(edge, d$contrast, d$tip_data, d$weight, d$levels, concavity = conc)
}
kernelTbr <- function(tree, d, phys) {
  edge <- Preorder(RenumberTips(tree, d$labels))[["edge"]]
  if (phys) Sys.setenv(TS_PHYS_REROOT = "1") else Sys.unsetenv("TS_PHYS_REROOT")
  res <- TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L, concavity = conc, unrooted = TRUE)
  Sys.unsetenv("TS_PHYS_REROOT")
  list(tree = structure(list(edge = res$edge, Nnode = d$nTip - 1L, tip.label = d$labels),
                        class = "phylo"), score = res$score)
}

d <- randomData(1000L + idx)
set.seed(7000L + idx); start <- RandomTree(d$phy, root = TRUE)
set.seed(idx); rD <- kernelTbr(start, d, FALSE)
set.seed(idx); rP <- kernelTbr(start, d, TRUE)
cat(sprintf("=== direct vs physical, nTip=%d tree#%d conc=%g ===\n", nTip, idx, conc))
cat(sprintf("DIRECT   : reported=%.5f  ts_fitch=%.5f\n", rD$score, scoreTree(rD$tree, d)))
cat(sprintf("PHYSICAL : reported=%.5f  ts_fitch=%.5f\n", rP$score, scoreTree(rP$tree, d)))
dlo <- scoreTree(rD$tree, d); plo <- scoreTree(rP$tree, d)
cat(sprintf("=> %s\n", if (plo < dlo - 1e-6) "PHYSICAL reaches lower => DIRECT-PATH gap (single-rooting)"
            else if (dlo < plo - 1e-6) "DIRECT lower (physical incomplete?!)"
            else "EQUAL => both miss it (deeper / plateau / common enumeration gap)"))
