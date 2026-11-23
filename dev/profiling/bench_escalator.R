# Smoke validation for stallEscalateFactor (the online stall-escalator).
# Does escalating ratchet perturbation on stall help close the TNT gap on a
# hard dataset? Standard-Fitch path ('-' -> '?'), nThreads=1, several seeds.
# Compares factor=1.0 (off, baseline) against escalation factors.
#   env: TREESEARCH_VTUNE_LIB (lib path), TS_DATASET (Wortley2006),
#        TS_REPS (30), TS_SECONDS (0 = rep-limited), TS_SEEDS ("1 2 3")
LIBDIR <- normalizePath(Sys.getenv("TREESEARCH_VTUNE_LIB"), winslash = "/")
suppressMessages(library(TreeSearch, lib.loc = LIBDIR))
suppressMessages(library(TreeTools))

dsName <- Sys.getenv("TS_DATASET", unset = "Wortley2006")
nReps  <- as.integer(Sys.getenv("TS_REPS", unset = "30"))
maxSec <- as.double(Sys.getenv("TS_SECONDS", unset = "0"))
seeds  <- as.integer(strsplit(Sys.getenv("TS_SEEDS", unset = "1 2 3"), " ")[[1]])
factors <- c(1.0, 1.5, 2.0)               # 1.0 = escalator off (baseline)

raw <- inapplicable.phyData[[dsName]]
m <- PhyDatToMatrix(raw, ambigNA = FALSE)
m[m == "-"] <- "?"                         # TNT-parity standard Fitch
dataset <- MatrixToPhyDat(m)
cat(sprintf("%s | %d tips | factors {%s} | seeds {%s} | reps %d | %ss\n\n",
            dsName, length(dataset), paste(factors, collapse = ", "),
            paste(seeds, collapse = ","), nReps,
            if (maxSec > 0) maxSec else "rep-limited"))

for (fac in factors) {
  scores <- integer(0)
  for (sd in seeds) {
    set.seed(sd)
    r <- suppressWarnings(MaximizeParsimony(
      dataset, maxReplicates = nReps, nThreads = 1L, strategy = "auto",
      maxSeconds = maxSec, verbosity = 0L,
      control = SearchControl(stallEscalateFactor = fac)))
    scores <- c(scores, attr(r, "score"))
  }
  cat(sprintf("factor %.1f : best=%d  median=%.1f  scores={%s}\n",
              fac, min(scores), median(scores),
              paste(scores, collapse = ",")))
}
