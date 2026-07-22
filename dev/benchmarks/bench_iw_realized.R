# IW realized-cost isolation: the PRODUCTION tbr_search loop (real bail logic +
# EW x4-flat vs IW scalar), instrumented via env TS_IW_TIMING. Prints per regime:
#   per_clip_us  = bail-INDEPENDENT IW precompute (extract_char_steps + iw_delta)
#   per_cand_ns  = bail-DEPENDENT realized per-candidate scan
# EW per_clip_us ~ 0 (no IW block); the EW/IW per_cand_ns ratio is the realized
# (not unbounded) per-candidate penalty that gates whether the gather/batching
# is a mission lever or a worst-case artifact.
#
# Usage: Rscript dev/benchmarks/bench_iw_realized.R [lib] [concavity]
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  CONC <- if (length(args) >= 2) as.numeric(args[[2]]) else 10
  # mode: "pure" (recode -> ?, has_na=FALSE, EW-comparable) or "native" (raw,
  # has_na=TRUE => NA+IW path, the datasets' real default).
  MODE <- if (length(args) >= 3) args[[3]] else "pure"
  library(TreeSearch, lib.loc = LIB); library(TreeTools)
}))
data("inapplicable.phyData", package = "TreeSearch")
Sys.setenv(TS_IW_TIMING = "1")

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

DATASETS <- c("Vinther2008", "Wortley2006", "Zanol2014", "Giles2015", "Dikow2009")
for (nm in DATASETS) {
  raw <- inapplicable.phyData[[nm]]
  phy <- if (identical(MODE, "native")) raw else recode_na_off(raw)
  p <- prep(phy)
  set.seed(7294)
  redge <- RandomTree(names(inapplicable.phyData[[nm]]), root = TRUE)$edge
  for (regime in c("EW", "IW")) {
    conc <- if (regime == "EW") -1.0 else CONC
    message(sprintf("=== %s n=%d ns=%d %s ===", nm, p$n, p$ns, regime))
    invisible(TreeSearch:::ts_tbr_search(redge, p$contrast, p$tip, p$weight,
                p$levels, maxHits = 1L, concavity = conc))
  }
}
