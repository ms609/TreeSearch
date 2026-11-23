#!/usr/bin/env Rscript
# Drift EW decision-flip census: quantify how far the deployed union-of-finals
# drift scorer diverges from the exact directional path, on the three gate
# datasets recoded to pure Fitch (EW).  Prints the DRIFT-SCANCHK line for the
# union (default) and exact paths.  Local, bounded -- severity probe only.
suppressMessages({
  library(TreeSearch, lib.loc = ".agent-hj")
  library(TreeTools)
})

fitchPhy <- function(p) {
  m <- PhyDatToMatrix(p, ambigNA = FALSE)
  m[m == "-"] <- "?"                       # inapplicable -> missing => standard Fitch
  MatrixToPhyDat(m)
}

makeDs <- function(dataset) {
  at <- attributes(dataset)
  list(contrast = at$contrast,
       tip_data = matrix(unlist(dataset, use.names = FALSE),
                         nrow = length(dataset), byrow = TRUE),
       weight = at$weight,
       levels = at$levels)
}

datasets <- c("Zhu2013", "Zanol2014", "Dikow2009")

for (nm in datasets) {
  phy <- fitchPhy(ReadAsPhyDat(file.path("data-raw", paste0(nm, ".nex"))))
  ds <- makeDs(phy)
  nTip <- length(phy)
  hasGap <- "-" %in% ds$levels
  cat(sprintf("\n=== %s : %d tips, %d chars(weighted), levels={%s} hasGap=%s ===\n",
              nm, nTip, sum(ds$weight), paste(ds$levels, collapse = ""), hasGap))

  # Production-representative start: random tree -> TBR to a local optimum, then
  # census drift from there (drift runs AFTER a hill-climb in real searches).
  set.seed(1L)
  start <- RandomTree(nTip, root = TRUE)
  # Align edge tip index i with tip_data row i (= names(phy)[i]).
  start$tip.label <- names(phy)
  tbr <- TreeSearch:::ts_tbr_search(start$edge, ds$contrast, ds$tip_data,
                                    ds$weight, ds$levels, maxHits = 5L)
  startTree <- list(edge = tbr$edge)
  cat(sprintf("  local-opt score after TBR = %.0f\n", tbr$score))

  runDrift <- function(exact) {
    if (exact) Sys.setenv(TS_DRIFT_EXACT = "1") else Sys.unsetenv("TS_DRIFT_EXACT")
    Sys.setenv(TS_DRIFT_SCANCHK = "1")
    set.seed(42L)
    r <- TreeSearch:::ts_drift_search(startTree$edge, ds$contrast, ds$tip_data,
                                      ds$weight, ds$levels,
                                      nCycles = 12L, afdLimit = 3L,
                                      rfdLimit = 0.1, maxHits = 1L)
    Sys.unsetenv("TS_DRIFT_SCANCHK"); Sys.unsetenv("TS_DRIFT_EXACT")
    cat(sprintf("  [%s] final drift score = %.0f (drift_moves=%d tbr_moves=%d)\n",
                if (exact) "EXACT" else "UNION", r$score,
                r$total_drift_moves, r$total_tbr_moves))
    invisible(r)
  }
  runDrift(FALSE)  # union (deployed default) -- severity
  runDrift(TRUE)   # exact -- wiring validation (applied_mismatch must be 0)
}
cat("\nDONE\n")
