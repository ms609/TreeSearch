# IW element isolation, NEAR-OPTIMAL regime (faithfulness companion to
# bench_iw_element_isolation.R, which scans a RANDOM tree = worst case for IW).
#
# Real search spends most candidate-scans near a local optimum, where needs_step
# is SPARSE and IW's per-set-bit gather loop is short. Here each regime scans the
# neighbourhood of ITS OWN TBR-converged tree (EW tree for EW, IW tree for IW) --
# the realistic production regime. If the random-tree 3x IW per-candidate penalty
# collapses toward ~1x here, the gather is a worst-case artifact, not a lever.
#
# Usage: Rscript dev/benchmarks/bench_iw_element_nearopt.R [lib] [concavity] [reps]
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  CONC <- if (length(args) >= 2) as.numeric(args[[2]]) else 10
  REPS <- if (length(args) >= 3) as.integer(args[[3]]) else 7L
  library(TreeSearch, lib.loc = LIB); library(TreeTools)
}))
data("inapplicable.phyData", package = "TreeSearch")

recode_na_off <- function(phy) {
  at <- attributes(phy); ct <- at$contrast; lv <- at$levels; al <- at$allLevels
  inapp_col <- which(lv == "-"); dash_row <- which(al == "-")
  if (length(dash_row)) ct[dash_row, ] <- 1
  if (length(inapp_col)) { ct <- ct[, -inapp_col, drop = FALSE]; lv <- lv[-inapp_col] }
  attr(phy, "contrast") <- ct; attr(phy, "levels") <- lv; phy
}
prep <- function(phy) {
  at <- attributes(phy)
  list(contrast = at$contrast,
       tip = matrix(unlist(phy, use.names = FALSE), nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, n = length(phy), ns = ncol(at$contrast))
}
bench_one <- function(p, edge, conc) {
  fr <- ci <- ind <- numeric(REPS); nc <- ncand <- NA_integer_
  for (r in seq_len(REPS)) {
    set.seed(1000 + r)
    res <- TreeSearch:::ts_bench_tbr_phases(edge, p$contrast, p$tip, p$weight,
                                            p$levels, integer(0), conc)
    fr[r] <- res$time_full_rescore_us; ci[r] <- res$time_clip_incr_us
    ind[r] <- res$time_indirect_us; nc <- res$n_clips; ncand <- res$n_candidates
  }
  list(full = median(fr), clip = median(ci), indirect = median(ind),
       n_clips = nc, n_cand = ncand)
}

DATASETS <- c("Vinther2008", "Wortley2006", "Zanol2014", "Giles2015", "Dikow2009")
cat(sprintf("IW element isolation NEAR-OPTIMAL | concavity=%g reps=%d\n", CONC, REPS))
cat(sprintf("%-13s %3s %2s | %-26s | %-22s | %-14s\n", "dataset", "n", "ns",
            "per-cand ns (EW/IW/Δ ratio)", "per-clip us (EW/IW)", "full us EW/IW"))

for (nm in DATASETS) {
  phy <- recode_na_off(inapplicable.phyData[[nm]]); p <- prep(phy)
  set.seed(7294)
  redge <- RandomTree(names(inapplicable.phyData[[nm]]), root = TRUE)$edge
  # TBR-converge each regime to its own near-optimal tree
  ew_opt <- TreeSearch:::ts_tbr_search(redge, p$contrast, p$tip, p$weight,
              p$levels, maxHits = 1L, concavity = -1.0)$edge
  iw_opt <- TreeSearch:::ts_tbr_search(redge, p$contrast, p$tip, p$weight,
              p$levels, maxHits = 1L, concavity = CONC)$edge
  ew <- bench_one(p, ew_opt, -1.0)
  iw <- bench_one(p, iw_opt, CONC)
  pc_ew <- ew$indirect * 1000 / ew$n_cand
  pc_iw <- iw$indirect * 1000 / iw$n_cand
  cat(sprintf("%-13s %3d %2d | %6.1f/%6.1f/%+5.1f (%.2fx) | %6.1f/%6.1f | %5.0f/%5.0f\n",
              nm, p$n, p$ns, pc_ew, pc_iw, pc_iw - pc_ew, pc_iw / pc_ew,
              ew$clip / ew$n_clips, iw$clip / iw$n_clips, ew$full, iw$full))
}
