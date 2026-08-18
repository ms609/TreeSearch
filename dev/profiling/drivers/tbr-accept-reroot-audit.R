# Score-identity battery for the TBR-rerooting accept path — issue #38 (T-300)
#
# Extending the dirty-set incremental rescore from SPR accepts to TBR-rerooting
# accepts claims to be SCORE-IDENTICAL: the incremental score must equal what
# full_rescore would have returned, so every accept/reject decision — and hence
# the whole search trajectory — must be unchanged.  Run this against the
# pre-patch and post-patch libraries and diff the CSVs.
#
# Pair it with TS_TBR_ACCEPTCHK=1, which makes the DLL cross-check every
# incremental accept against full_rescore in-flight and abort on drift.  That is
# the guard the reverted first attempt (b7303ee5, systematic delta = -3) lacked.
#
# Usage:
#   TS_AUDIT_OUT=base.csv  Rscript dev/profiling/drivers/tbr-accept-reroot-audit.R .agent-i38
#   TS_TBR_ACCEPTCHK=1 TS_AUDIT_OUT=patched.csv \
#     Rscript dev/profiling/drivers/tbr-accept-reroot-audit.R .agent-i38b

args <- commandArgs(trailingOnly = TRUE)
libDir <- if (length(args) >= 1L) args[[1]] else ".agent-i38b"
library(TreeSearch, lib.loc = libDir)
library(TreeTools, quietly = TRUE)

MakeData <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tipData = matrix(unlist(dataset, use.names = FALSE),
                     nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    nTip = length(dataset)
  )
}

# Weak-signal random matrices accept long chains of moves, which is what drives
# reroot accepts; the real matrices add NA blocks and realistic state counts.
cases <- list()
for (nTip in c(12L, 18L, 25L)) {
  set.seed(1000 + nTip)
  mat <- matrix(sample(0:3, nTip * 8L, replace = TRUE), nrow = nTip,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  cases[[paste0("rand", nTip)]] <- MatrixToPhyDat(mat)
}
data("inapplicable.phyData", package = "TreeSearch")
for (nm in c("Longrich2010", "Vinther2008", "Sansom2010", "DeAssis2011")) {
  cases[[nm]] <- inapplicable.phyData[[nm]]
}

rows <- list()
for (nm in names(cases)) {
  dataset <- cases[[nm]]
  d <- MakeData(dataset)
  minSteps <- as.integer(MinimumLength(dataset, compress = TRUE))
  for (mode in c("EW", "IW")) {
    searchConcavity <- if (mode == "EW") -1 else 10
    scoreConcavity <- if (mode == "EW") Inf else 10
    ms <- if (mode == "EW") integer(0) else minSteps
    for (start in c(1, 17, 88, 256, 777)) {
      tree <- as.phylo(start, d$nTip)
      set.seed(start)
      res <- TreeSearch:::ts_tbr_search(
        tree$edge, d$contrast, d$tipData, d$weight, d$levels,
        maxHits = 50L, min_steps = ms, concavity = searchConcavity)
      independent <- TreeSearch:::ts_fitch_score(
        res$edge, d$contrast, d$tipData, d$weight, d$levels,
        min_steps = ms, concavity = scoreConcavity)
      rows[[length(rows) + 1L]] <- data.frame(
        case = nm, mode = mode, start = start,
        score = res$score, independent = independent,
        nAccepted = res$n_accepted, nEvaluated = res$n_evaluated,
        stringsAsFactors = FALSE)
    }
  }
}

tab <- do.call(rbind, rows)
tab$drift <- tab$score - tab$independent
bad <- tab[abs(tab$drift) > 1e-9, , drop = FALSE]

cat(sprintf("cells: %d | total accepts: %d | in-flight audit: %s\n",
            nrow(tab), sum(tab$nAccepted),
            if (nzchar(Sys.getenv("TS_TBR_ACCEPTCHK"))) "ON" else "off"))
if (nrow(bad)) {
  cat("SCORE DRIFT vs independent recomputation:\n")
  print(bad, row.names = FALSE)
} else {
  cat("all reported scores match an independent full recomputation\n")
}

outFile <- Sys.getenv("TS_AUDIT_OUT", unset = "")
if (nzchar(outFile)) {
  write.csv(tab[, c("case", "mode", "start", "score", "nAccepted", "nEvaluated")],
            outFile, row.names = FALSE)
  cat("wrote", outFile, "\n")
}
if (nrow(bad)) quit(status = 1L)
