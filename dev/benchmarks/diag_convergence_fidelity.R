# Consensus-fidelity gate for the xmult-style convergence stop (consensusStableReps).
# Decides whether to DEFAULT the stop, using an UNBIASED reference:
#
#  (1) Leave-one-out union consensus: ref_s = strict consensus of the union of
#      the OTHER seeds' full-run MPTs (so a full run is never scored against a
#      set it belongs to).  Ship-clear if mean consCID(cs->ref) <= mean(full->ref):
#      early-stopping then costs nothing the seed lottery wasn't already costing.
#  (2) Resolution direction: internal-node count of each consensus.  If the
#      early-stop consensus is MORE resolved than the full-run consensus it
#      overstates support (harmful); equal-or-less is conservative (harmless).
#  (3) Zhu stress: extra seeds at the ship K, score-loss must stay 0 (Zhu is the
#      high-plateauSafeK case where cs3 already lost +1).
#
# Scope: `thorough` only (all tuning was on thorough).
# Env: TS_LIB (default .agent-stop), NSEED (default 3), SHIPK (default 6),
#      ZHU_EXTRA (default 9 extra Zhu seeds for the stress test).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

nseed   <- as.integer(Sys.getenv("NSEED", "3"))
shipK   <- as.integer(Sys.getenv("SHIPK", "6"))
zhuExtra<- as.integer(Sys.getenv("ZHU_EXTRA", "9"))
armsK   <- c(5L, shipK)          # cs5 (continuity) + ship candidate
datasets <- c("Wortley2006", "Zanol2014", "Zhu2013", "Giles2015")
target   <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
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
.run <- function(phy, seed, csReps) {
  set.seed(seed)
  t <- system.time(
    r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
           maxSeconds = 600, nThreads = 1L, verbosity = 0L,
           consensusStableReps = csReps)))
  attr(r, "wall") <- as.double(t["elapsed"]); r
}

# --- Phases 1+2: fidelity + resolution on 4 datasets x nseed -------------------
fid <- list()
for (nm in datasets) {
  phy <- .phy(nm)
  fullTrees <- vector("list", nseed); fullScore <- numeric(nseed)
  armTrees  <- list()                # [[K]][[seed]]
  for (K in armsK) armTrees[[as.character(K)]] <- vector("list", nseed)
  for (s in seq_len(nseed)) {
    rf <- .run(phy, s, 0L)
    fullTrees[[s]] <- rf; fullScore[s] <- min(as.double(attr(rf, "score")))
    for (K in armsK) armTrees[[as.character(K)]][[s]] <- .run(phy, s, K)
  }
  # Leave-one-out union reference + comparisons
  for (s in seq_len(nseed)) {
    others <- setdiff(seq_len(nseed), s)
    pool <- do.call(c, lapply(others, function(j) {
      tj <- fullTrees[[j]]; if (inherits(tj, "phylo")) list(tj) else tj
    }))
    class(pool) <- "multiPhylo"
    refCons <- .strict(pool)
    refNode <- refCons$Nnode
    fullCons <- .strict(fullTrees[[s]])
    cidFull  <- as.double(ClusteringInfoDist(fullCons, refCons, normalize = TRUE))
    for (K in armsK) {
      ar <- armTrees[[as.character(K)]][[s]]
      arCons <- .strict(ar)
      cidArm <- as.double(ClusteringInfoDist(arCons, refCons, normalize = TRUE))
      fid[[length(fid) + 1L]] <- data.frame(
        dataset = nm, seed = s, csReps = K,
        fullScore = fullScore[s], armScore = min(as.double(attr(ar, "score"))),
        scoreLoss = min(as.double(attr(ar, "score"))) - fullScore[s],
        fullWall = round(attr(fullTrees[[s]], "wall"), 1),
        armWall  = round(attr(ar, "wall"), 1),
        wallFrac = round(attr(ar, "wall") / attr(fullTrees[[s]], "wall"), 2),
        cidFull2ref = round(cidFull, 4), cidArm2ref = round(cidArm, 4),
        refNode = refNode, fullNode = fullCons$Nnode, armNode = arCons$Nnode,
        # >0 => arm MORE resolved than full (overstates support = harmful)
        nodeDelta = arCons$Nnode - fullCons$Nnode)
      cat(sprintf(paste0("%-12s s%d cs%d: loss %+.0f | wall %.0f%% | ",
                  "cid(full->ref)=%.3f cid(cs->ref)=%.3f | nodes full=%d cs=%d (%+d)\n"),
                  nm, s, K, min(as.double(attr(ar,"score"))) - fullScore[s],
                  100 * attr(ar,"wall")/attr(fullTrees[[s]],"wall"),
                  cidFull, cidArm, fullCons$Nnode, arCons$Nnode,
                  arCons$Nnode - fullCons$Nnode))
    }
  }
}
fdf <- do.call(rbind, fid)
write.csv(fdf, file.path(Sys.getenv("OUTDIR", "dev/benchmarks"),
                         "convergence_fidelity.csv"), row.names = FALSE)
cat("\n=== fidelity median by dataset x csReps ===\n")
agg <- aggregate(cbind(scoreLoss, wallFrac, cidFull2ref, cidArm2ref, nodeDelta) ~
                   dataset + csReps, fdf, median)
print(agg[order(agg$dataset, agg$csReps), ], row.names = FALSE)

# --- Phase 3: Zhu stress at ship K (score-loss only) ---------------------------
cat(sprintf("\n=== Zhu stress: seeds %d..%d at cs%d ===\n", nseed + 1L,
            nseed + zhuExtra, shipK))
phyZ <- .phy("Zhu2013"); zloss <- integer(0)
for (s in (nseed + 1L):(nseed + zhuExtra)) {
  rz <- .run(phyZ, s, shipK)
  sc <- min(as.double(attr(rz, "score")))
  zloss <- c(zloss, sc - target[["Zhu2013"]])
  cat(sprintf("Zhu s%d cs%d: score %.0f (%+.0f) | reps %d wall %.1fs\n",
              s, shipK, sc, sc - target[["Zhu2013"]],
              attr(rz, "replicates"), attr(rz, "wall")))
}
cat(sprintf("\nZhu stress at cs%d: max loss over %d extra seeds = %+d (want 0)\n",
            shipK, zhuExtra, max(zloss)))
