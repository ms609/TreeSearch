# Re-gate Lever 1 on the MAIN checkout: reproduce the subagent's bit-identical
# table.  If score AND candidates_evaluated match the documented gate values
# exactly, this build is byte-identical to the validated mod build (which the
# subagent showed == the unmodified base over 12 runs).
suppressMessages({
  library(TreeSearch, lib.loc = Sys.getenv("TS_LIB", ".agent-l1"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

# Documented gate values (subagent report, base == mod, reps=3):
gate <- list(
  "Wortley2006:1" = c(479, 3020328),
  "Wortley2006:2" = c(479, 3294030),
  "Zhu2013:1"     = c(625, 38297730),
  "Zhu2013:2"     = c(624, 40310443),
  "Zanol2014:1"   = c(1261, 40629194),
  "Zanol2014:2"   = c(1261, 36348409))

allok <- TRUE
for (key in names(gate)) {
  parts <- strsplit(key, ":")[[1]]
  nm <- parts[1]; sd <- as.integer(parts[2])
  phy <- fc(inapplicable.phyData[[nm]])
  set.seed(sd)
  r <- suppressWarnings(MaximizeParsimony(
    phy, maxReplicates = 3L, nThreads = 1L, strategy = "auto", verbosity = 0L))
  sc <- as.numeric(attr(r, "score"))
  ca <- as.numeric(attr(r, "candidates_evaluated"))
  exp_sc <- gate[[key]][1]; exp_ca <- gate[[key]][2]
  ok <- isTRUE(sc == exp_sc) && isTRUE(ca == exp_ca)
  allok <- allok && ok
  cat(sprintf("%-14s score %g (exp %g) %s | cand %.0f (exp %.0f) %s => %s\n",
              key, sc, exp_sc, ifelse(sc == exp_sc, "ok", "DIFF"),
              ca, exp_ca, ifelse(ca == exp_ca, "ok", "DIFF"),
              ifelse(ok, "MATCH", "*** MISMATCH ***")))
}
cat(if (allok) "\nALL 6 BIT-IDENTICAL TO VALIDATED BUILD\n" else "\n*** GATE FAILED ***\n")
