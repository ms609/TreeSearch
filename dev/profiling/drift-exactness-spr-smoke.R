#!/usr/bin/env Rscript
# Validate the spr_search exact path (TS_SPR_EXACT): the returned score must
# equal an INDEPENDENT full Fitch re-score of the returned tree (wiring), and the
# hill-climb must improve on the start.  Also contrasts union vs exact reach --
# spr_search is a pure hill-climb with no compensating exact-TBR phase, so the
# union over-count (which trips the `dominated` gate at ts_search.cpp:379 and
# skips improving moves) causes premature convergence the exact path avoids.
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("GATE_LIB", ".agent-hj"))
  library(TreeTools)
})
fitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m)
}
phy <- fitchPhy(ReadAsPhyDat("data-raw/Zhu2013.nex"))
at <- attributes(phy)
ds <- list(contrast = at$contrast,
           tip_data = matrix(unlist(phy, use.names = FALSE),
                             nrow = length(phy), byrow = TRUE),
           weight = at$weight, levels = at$levels)
indep <- function(edge) TreeSearch:::ts_fitch_score(edge, ds$contrast, ds$tip_data,
                          ds$weight, ds$levels, min_steps = integer(0), concavity = -1)
for (s in 1:5) {
  set.seed(s)
  st <- RandomTree(length(phy), root = TRUE); st$tip.label <- names(phy)
  start_len <- indep(st$edge)
  for (ex in c(FALSE, TRUE)) {
    if (ex) Sys.setenv(TS_SPR_EXACT = "1") else Sys.unsetenv("TS_SPR_EXACT")
    set.seed(100 + s)
    r <- TreeSearch:::ts_spr_search(st$edge, ds$contrast, ds$tip_data,
                                    ds$weight, ds$levels, maxHits = 5L)
    Sys.unsetenv("TS_SPR_EXACT")
    chk <- indep(r$edge)
    cat(sprintf("seed=%d %-5s start=%.0f returned=%.0f indep=%.0f MATCH=%s improved=%s\n",
                s, if (ex) "exact" else "union", start_len, r$score, chk,
                abs(r$score - chk) < 1e-6, r$score <= start_len))
  }
}
