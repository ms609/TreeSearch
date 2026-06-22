# Lever 1 validation: IW reroot x4 batching.
#   byte-identity: score AND n_evaluated must match scalar (TS_IW_NOX4=1).
#   speedup:       REROOT per_cand_ns (TS_IW_TIMING) x4-on vs x4-off.
# Pure-IW (recode -> ?). One ts_tbr_search to convergence per (dataset, mode).
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  CONC <- if (length(args) >= 2) as.numeric(args[[2]]) else 10
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
       weight = at$weight, levels = at$levels, n = length(phy))
}
# set.seed immediately before each call: ts_tbr_search seeds its internal
# main_edges partial-shuffle from R's RNG (make_rng -> unif_rand), so x4-on and
# x4-off must share the SAME R seed or they take different (valid) trajectories.
run <- function(p, edge) {
  set.seed(424242L)
  TreeSearch:::ts_tbr_search(edge, p$contrast, p$tip, p$weight, p$levels,
                             maxHits = 1L, concavity = CONC)
}

DATASETS <- c("Vinther2008", "Wortley2006", "Zanol2014", "Giles2015", "Dikow2009")
for (nm in DATASETS) {
  phy <- recode_na_off(inapplicable.phyData[[nm]]); p <- prep(phy)
  set.seed(7294)
  edge <- RandomTree(names(inapplicable.phyData[[nm]]), root = TRUE)$edge
  Sys.setenv(TS_IW_NOX4 = "1"); message(sprintf("=== %s SCALAR ===", nm)); off <- run(p, edge)
  Sys.unsetenv("TS_IW_NOX4");    message(sprintf("=== %s X4 ===", nm));     on  <- run(p, edge)
  ok <- isTRUE(all.equal(off$score, on$score)) && off$n_evaluated == on$n_evaluated
  message(sprintf("RESULT %-12s score %.4f/%.4f neval %d/%d => %s",
                  nm, off$score, on$score, off$n_evaluated, on$n_evaluated,
                  if (ok) "IDENTICAL" else "*** DIFF ***"))
}
