# How well does each engine SAMPLE tree space (not just "what's the best score")?
# Reframes the consensus-fidelity question: our early-stop consensus is more
# resolved than our OWN exhaustive full run -- but the engine we are matching is
# TNT, whose xmult also self-terminates.  So the fair benchmark is TNT's
# sampling, not our gold-plated full run.
#
# Tree-space sampling is BITNESS-INDEPENDENT (only wall-clock needs Hamilton), so
# local 32-bit TNT is valid here.  For each dataset x seed we collect the MPT set
# from three methods and compare the STRICT-CONSENSUS RESOLUTION (internal-node
# count) -- fewer nodes = more conservative = more thoroughly sampled plateau:
#   TNT  : xmult=level 10, hold 10000, best  (representative thorough user)
#   TSf  : TreeSearch thorough, full run (current default; no early stop)
#   TScs : TreeSearch thorough + consensusStableReps=6 (the proposed early stop)
#
# A "best available" reference consensus is built from the UNION of all three
# methods' MPTs (all seeds); each method's distance to it (ClusteringInfoDist)
# and its node count are reported.  Verdict logic in the trailer.
#
# Env: TS_LIB (default .agent-stop), NSEED (default 3),
#      TNT_EXE (default local 32-bit 1.6).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

nseed   <- as.integer(Sys.getenv("NSEED", "3"))
TNT_EXE <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
datasets <- c("Wortley2006", "Zanol2014", "Zhu2013", "Giles2015")
target   <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
data("inapplicable.phyData", package = "TreeSearch")
wd <- file.path(tempdir(), "tntsamp"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

.phy <- function(nm) {
  m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE); m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
.strict <- function(trees) {
  if (is.null(trees) || length(trees) == 0L) return(NULL)
  if (inherits(trees, "phylo")) return(trees)
  if (length(trees) == 1L) return(trees[[1]])
  ape::consensus(trees, p = 1)
}
.asMP <- function(tr) {        # normalise to multiPhylo list of phylo
  if (is.null(tr)) return(list())
  if (inherits(tr, "phylo")) return(list(tr))
  unclass(tr)
}

# --- TNT: xmult=level 10, retain MPTs, save all trees, read back --------------
runTNT <- function(phy, seed) {
  datafile <- file.path(wd, "d.tnt"); out <- file.path(wd, "tntout.tre")
  if (file.exists(out)) file.remove(out)
  WriteTntCharacters(phy, datafile)
  cmds <- c("mxram 1024;", sprintf("proc %s;", basename(datafile)),
            "hold 10000;", sprintf("rseed %d;", seed),
            "xmult=level 10;", "best;",
            "tsave *tntout.tre;", "save;", "tsave/;", "quit;")
  old <- setwd(wd); on.exit(setwd(old))
  system2(TNT_EXE, input = cmds, stdout = FALSE, stderr = FALSE)
  tr <- tryCatch(ReadTntTree("tntout.tre"), error = function(e) NULL)
  tr
}

allMPT <- list()           # method -> dataset -> seed -> multiPhylo (for union ref)
rows <- list()
for (nm in datasets) {
  phy <- .phy(nm); tgt <- target[[nm]]
  for (s in seq_len(nseed)) {
    # TNT
    tnt <- runTNT(phy, s)
    tntMP <- .asMP(tnt)
    tntSc <- if (length(tntMP)) min(vapply(tntMP, function(t) TreeLength(t, phy), 0)) else NA
    # TreeSearch full
    set.seed(s)
    tsf <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
             maxSeconds = 600, nThreads = 1L, verbosity = 0L))
    # TreeSearch + early stop
    set.seed(s)
    tscs <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
              maxSeconds = 600, nThreads = 1L, verbosity = 0L,
              consensusStableReps = 6L))
    allMPT[["TNT"]][[nm]][[s]]  <- tntMP
    allMPT[["TSf"]][[nm]][[s]]  <- .asMP(tsf)
    allMPT[["TScs"]][[nm]][[s]] <- .asMP(tscs)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = nm, seed = s, target = tgt,
      tntScore = round(tntSc), tntMPT = length(tntMP),
      tntNode = { c <- .strict(tnt); if (is.null(c)) NA else c$Nnode },
      tsfScore = min(as.double(attr(tsf, "score"))), tsfMPT = length(tsf),
      tsfNode = .strict(tsf)$Nnode,
      tscsScore = min(as.double(attr(tscs, "score"))), tscsMPT = length(tscs),
      tscsNode = .strict(tscs)$Nnode)
    cat(sprintf(paste0("%-12s s%d: TNT %.0f (n=%d, nodes=%s) | ",
                "TSfull %.0f (n=%d, nodes=%d) | TScs6 %.0f (n=%d, nodes=%d)\n"),
                nm, s, tntSc, length(tntMP),
                { c <- .strict(tnt); if (is.null(c)) "NA" else c$Nnode },
                min(as.double(attr(tsf,"score"))), length(tsf), .strict(tsf)$Nnode,
                min(as.double(attr(tscs,"score"))), length(tscs), .strict(tscs)$Nnode))
  }
}
df <- do.call(rbind, rows)

# --- Union reference per dataset + CID of each method's per-dataset consensus --
cidRows <- list()
for (nm in datasets) {
  uni <- do.call(c, lapply(c("TNT", "TSf", "TScs"), function(meth)
    do.call(c, lapply(seq_len(nseed), function(s) allMPT[[meth]][[nm]][[s]]))))
  uni <- uni[!vapply(uni, is.null, TRUE)]
  class(uni) <- "multiPhylo"
  refCons <- .strict(uni); refNode <- refCons$Nnode
  for (meth in c("TNT", "TSf", "TScs")) {
    methAll <- do.call(c, lapply(seq_len(nseed), function(s) allMPT[[meth]][[nm]][[s]]))
    methAll <- methAll[!vapply(methAll, is.null, TRUE)]; class(methAll) <- "multiPhylo"
    mc <- .strict(methAll)
    cid <- tryCatch(as.double(ClusteringInfoDist(mc, refCons, normalize = TRUE)),
                    error = function(e) NA_real_)
    cidRows[[length(cidRows) + 1L]] <- data.frame(
      dataset = nm, method = meth, node = mc$Nnode, refNode = refNode,
      cid2union = round(cid, 4))
  }
}
cdf <- do.call(rbind, cidRows)
write.csv(df,  file.path(Sys.getenv("OUTDIR", "dev/benchmarks"), "treespace_sampling.csv"), row.names = FALSE)
write.csv(cdf, file.path(Sys.getenv("OUTDIR", "dev/benchmarks"), "treespace_cid.csv"), row.names = FALSE)

cat("\n=== resolution (median internal nodes; lower = more conservative sampling) ===\n")
print(aggregate(cbind(tntNode, tsfNode, tscsNode) ~ dataset, df, median), row.names = FALSE)
cat("\n=== pooled-by-dataset consensus vs union-of-all-methods reference ===\n")
print(cdf[order(cdf$dataset, cdf$method), ], row.names = FALSE)
cat("\nVERDICT: if tscsNode is between TNT and TSfull (i.e. TScs <= TNT), our early\n",
    "stop samples tree space at least as conservatively as TNT -> no regression vs\n",
    "the engine we are matching -> ship the stop. If TNT ~ TSfull << TScs, TNT truly\n",
    "samples better and the full run is worth keeping as the non-early-stop path.\n", sep = "")
