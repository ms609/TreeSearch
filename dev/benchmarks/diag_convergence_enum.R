# Does cheap MPT enumeration recover consensus completeness when the convergence
# stop is ON?  Enumeration (TBR plateau walk in finish:) is SKIPPED while the
# pool is at its cap (100), so the early-stop's clustered 100-tree pool stays
# over-resolved.  Test: give the pool room (poolMaxSize up) so enumeration runs
# and injects within-island diversity; measure whether the strict-consensus node
# count collapses back toward the full run's, and at what wall cost.
#
# Reference for CID = leave-one-out union of the FULL-run MPTs (unbiased).
# "Truth" is the MOST collapsed consensus (fewest internal nodes); lower armNode
# = closer to truth = better.  Win if cs6+enum: score 0-loss, armNode <= fullNode
# (or close), wall still << full.  Else it's a genuine speed/consensus tradeoff.
#
# Env: TS_LIB (default .agent-stop), NSEED (default 3), POOL (default 400).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

nseed <- as.integer(Sys.getenv("NSEED", "3"))
pool  <- as.integer(Sys.getenv("POOL", "400"))
datasets <- c("Zanol2014", "Zhu2013")   # the two over-resolved cases
data("inapplicable.phyData", package = "TreeSearch")

.phy <- function(nm) {
  m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE); m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
.strict <- function(trees) {
  if (inherits(trees, "phylo")) return(trees)
  if (length(trees) == 1L) return(trees[[1]])
  ape::consensus(trees, p = 1)
}

rows <- list()
for (nm in datasets) {
  phy <- .phy(nm)
  fullTrees <- vector("list", nseed); armList <- vector("list", nseed)
  fullWall <- armWall <- numeric(nseed); fullSc <- armSc <- numeric(nseed)
  for (s in seq_len(nseed)) {
    set.seed(s)
    tf <- system.time(rf <- suppressWarnings(MaximizeParsimony(phy,
            strategy = "thorough", maxSeconds = 600, nThreads = 1L,
            verbosity = 0L)))                                    # full, pool 100
    set.seed(s)
    ta <- system.time(ra <- suppressWarnings(MaximizeParsimony(phy,
            strategy = "thorough", maxSeconds = 600, nThreads = 1L,
            verbosity = 0L, consensusStableReps = 6L,
            poolMaxSize = pool, enumTimeFraction = 0.3)))        # stop + enum room
    fullTrees[[s]] <- rf; armList[[s]] <- ra
    fullWall[s] <- tf["elapsed"]; armWall[s] <- ta["elapsed"]
    fullSc[s] <- min(as.double(attr(rf, "score")))
    armSc[s]  <- min(as.double(attr(ra, "score")))
  }
  for (s in seq_len(nseed)) {
    others <- setdiff(seq_len(nseed), s)
    refPool <- do.call(c, lapply(others, function(j) {
      tj <- fullTrees[[j]]; if (inherits(tj, "phylo")) list(tj) else tj }))
    class(refPool) <- "multiPhylo"
    refCons <- .strict(refPool)
    fullCons <- .strict(fullTrees[[s]]); armCons <- .strict(armList[[s]])
    cidFull <- as.double(ClusteringInfoDist(fullCons, refCons, normalize = TRUE))
    cidArm  <- as.double(ClusteringInfoDist(armCons,  refCons, normalize = TRUE))
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = nm, seed = s, scoreLoss = armSc[s] - fullSc[s],
      fullWall = round(fullWall[s], 1), armWall = round(armWall[s], 1),
      wallFrac = round(armWall[s] / fullWall[s], 2),
      fullMPT = length(fullTrees[[s]]), armMPT = length(armList[[s]]),
      fullNode = fullCons$Nnode, armNode = armCons$Nnode,
      nodeDelta = armCons$Nnode - fullCons$Nnode,
      cidFull2ref = round(cidFull, 4), cidArm2ref = round(cidArm, 4))
    cat(sprintf(paste0("%-12s s%d: loss %+.0f | wall %.0f%% (%.1f->%.1fs) | ",
                "MPT %d->%d | nodes full=%d arm=%d (%+d) | cid full=%.3f arm=%.3f\n"),
                nm, s, armSc[s] - fullSc[s], 100 * armWall[s] / fullWall[s],
                fullWall[s], armWall[s], length(fullTrees[[s]]), length(armList[[s]]),
                fullCons$Nnode, armCons$Nnode, armCons$Nnode - fullCons$Nnode,
                cidFull, cidArm))
  }
}
df <- do.call(rbind, rows)
write.csv(df, file.path(Sys.getenv("OUTDIR", "dev/benchmarks"),
                        "convergence_enum.csv"), row.names = FALSE)
cat("\n=== median by dataset (cs6 + poolMaxSize=", pool, " + enumFrac 0.3) ===\n", sep = "")
agg <- aggregate(cbind(scoreLoss, wallFrac, armMPT, nodeDelta, cidArm2ref, cidFull2ref) ~
                   dataset, df, median)
print(agg, row.names = FALSE)
cat("\nRecovery WIN if nodeDelta ~<=0 and wallFrac << 1. Else: genuine tradeoff -> ask user.\n")
