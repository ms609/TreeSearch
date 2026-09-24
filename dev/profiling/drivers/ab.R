# A/B gate for the compute_insertion_edge_sets buffer-hoist (Step 1).
#
# Behaviour-neutrality (correctness, paramount — this is the quality-fix fn):
# a pure buffer hoist must NOT change the search trajectory, so BOTH the final
# score AND candidates_evaluated must be BIT-IDENTICAL between baseline and the
# modified build, on the EW-fitch path (changed via tbr_search) AND the raw NA
# path (which still hits the fn via Wagner construction + sectorial).
#
# Run twice (AB_TAG=base TS_LIB=.agent-p0 ; AB_TAG=step1 TS_LIB=.bench-step1),
# appending to OUT_CSV, then diff.  nThreads=1L (candidates_evaluated only
# populated serially) + set.seed for determinism.

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB"), winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")

tag <- Sys.getenv("AB_TAG", "base")
out <- Sys.getenv("OUT_CSV", "dev/profiling/ab_latest.csv")
fc <- function(phy) { m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS", "Wortley2006 Zhu2013 Zanol2014")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2")), "\\s+")[[1]])
# NA path is owned by another agent; this gate covers the EW (fitch) path only.
objsN <- strsplit(trimws(Sys.getenv("TS_OBJS", "fitch")), "\\s+")[[1]]

rows <- list()
for (nm in dsN) {
  raw <- inapplicable.phyData[[nm]]; fit <- fc(raw)
  for (obj in objsN) {
    phy <- if (obj == "fitch") fit else raw
    for (sd in seeds) {
      set.seed(sd); t0 <- Sys.time()
      r <- suppressWarnings(MaximizeParsimony(
        phy, maxReplicates = 3L, nThreads = 1L, strategy = "auto", verbosity = 0L))
      w <- as.double(difftime(Sys.time(), t0, units = "secs"))
      rows[[length(rows) + 1]] <- data.frame(
        tag = tag, dataset = nm, obj = obj, seed = sd,
        score = attr(r, "score"), cand = attr(r, "candidates_evaluated"),
        wall = round(w, 2), stringsAsFactors = FALSE)
    }
  }
}
res <- do.call(rbind, rows)
write.table(res, out, sep = ",", row.names = FALSE,
            col.names = !file.exists(out), append = file.exists(out))
cat(sprintf("%s: %d runs, total wall %.1f s\n", tag, nrow(res), sum(res$wall)))
