# Decisive tree-space sampling comparison, isolating THREE confounded effects on
# the strict-consensus resolution (internal-node count; fewer = more conservative
# = more thoroughly sampled plateau):
#   (A) pool-cap under-sampling : TS full pool 100  vs  TS full pool 10000
#   (B) early-stop island deficit: TS full  vs  TS cs6   (at the SAME pool)
#   (C) the benchmark           : TNT xmult=level 10, hold 10000
#
# Tree-space sampling is bitness-independent, so local 32-bit TNT is valid.
# Reference = strict consensus of the UNION of all methods' MPTs (per dataset);
# each method's consensus node count + ClusteringInfoDist-to-union reported.
#
# Verdict hinges on where TNT sits:
#  - TNT ~ TS-full-pool100 (both over-resolved) => TNT is ALSO island-limited in
#    fast mode; our early stop is no worse than the engine we match => SHIP stop.
#  - TNT ~ TS-full-pool10000 (well collapsed) << TS-cs6 => TNT samples better;
#    keep the full/large-pool path for conservative consensus (stop = opt-in).
#  - TS-full-pool100 >> TS-full-pool10000 => our DEFAULT pool of 100 under-samples
#    vs TNT regardless of the stop => raise poolMaxSize (separate, important fix).
#
# Env: TS_LIB (default .agent-stop), NSEED (default 3),
#      TNT_EXE (local 32-bit), BIGPOOL (default 10000), MAXSEC (default 300).
.libPaths(c(Sys.getenv("TS_LIB", ".agent-stop"), .libPaths()))
suppressMessages({ library(TreeSearch); library(TreeTools); library(TreeDist) })

nseed   <- as.integer(Sys.getenv("NSEED", "3"))
bigPool <- as.integer(Sys.getenv("BIGPOOL", "10000"))
maxSec  <- as.integer(Sys.getenv("MAXSEC", "300"))
TNT_EXE <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
datasets <- c("Zanol2014", "Zhu2013", "Wortley2006", "Giles2015")
target   <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)
data("inapplicable.phyData", package = "TreeSearch")
wd <- file.path(tempdir(), "tspool"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

.phy <- function(nm) {
  m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE); m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}
.strict <- function(tr) {
  if (is.null(tr) || length(tr) == 0L) return(NULL)
  if (inherits(tr, "phylo")) return(tr)
  if (length(tr) == 1L) return(tr[[1]])
  ape::consensus(tr, p = 1)
}
.asMP <- function(tr) { if (is.null(tr)) list() else if (inherits(tr, "phylo")) list(tr) else unclass(tr) }

runTNT <- function(phy, seed) {
  datafile <- file.path(wd, "d.tnt"); out <- file.path(wd, "o.tre")
  if (file.exists(out)) file.remove(out)
  WriteTntCharacters(phy, datafile)
  cmds <- c("mxram 1024;", sprintf("proc %s;", basename(datafile)), "hold 10000;",
            sprintf("rseed %d;", seed), "xmult=level 10;", "best;",
            "tsave *o.tre;", "save;", "tsave/;", "quit;")
  old <- setwd(wd); on.exit(setwd(old))
  t <- system.time(system2(TNT_EXE, input = cmds, stdout = FALSE, stderr = FALSE))
  tr <- tryCatch(ReadTntTree("o.tre"), error = function(e) NULL)
  list(tr = tr, wall = as.double(t["elapsed"]))
}
runTS <- function(phy, seed, csReps, poolSize) {
  set.seed(seed)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(phy, strategy = "thorough",
         maxSeconds = maxSec, nThreads = 1L, verbosity = 0L,
         consensusStableReps = csReps, poolMaxSize = poolSize)))
  list(tr = r, wall = as.double(t["elapsed"]))
}

methods <- list(
  TNT          = function(phy, s) { x <- runTNT(phy, s); list(mp = .asMP(x$tr), wall = x$wall) },
  TSf_p100     = function(phy, s) { x <- runTS(phy, s, 0L, 100L);     list(mp = .asMP(x$tr), wall = x$wall) },
  TSf_pBig     = function(phy, s) { x <- runTS(phy, s, 0L, bigPool);  list(mp = .asMP(x$tr), wall = x$wall) },
  TScs6_pBig   = function(phy, s) { x <- runTS(phy, s, 6L, bigPool);  list(mp = .asMP(x$tr), wall = x$wall) }
)

allMP <- list(); rows <- list()
for (nm in datasets) {
  phy <- .phy(nm); tgt <- target[[nm]]
  for (s in seq_len(nseed)) {
    rec <- list(dataset = nm, seed = s, target = tgt)
    for (mn in names(methods)) {
      res <- methods[[mn]](phy, s)
      mp <- res$mp; allMP[[mn]][[nm]][[s]] <- mp
      sc <- if (length(mp)) min(vapply(mp, function(t) TreeLength(t, phy), 0)) else NA_real_
      cons <- .strict(mp)
      rec[[paste0(mn, "_sc")]]   <- round(sc)
      rec[[paste0(mn, "_n")]]    <- length(mp)
      rec[[paste0(mn, "_node")]] <- if (is.null(cons)) NA_integer_ else cons$Nnode
      rec[[paste0(mn, "_wall")]] <- round(res$wall, 1)
    }
    rows[[length(rows) + 1L]] <- as.data.frame(rec)
    cat(sprintf("%-12s s%d | nodes: TNT=%s TSf100=%s TSfBig=%s cs6Big=%s | n: %s/%s/%s/%s | sc TNT=%s\n",
        nm, s, rec$TNT_node, rec$TSf_p100_node, rec$TSf_pBig_node, rec$TScs6_pBig_node,
        rec$TNT_n, rec$TSf_p100_n, rec$TSf_pBig_n, rec$TScs6_pBig_n, rec$TNT_sc))
  }
}
df <- do.call(rbind, rows)
write.csv(df, file.path(Sys.getenv("OUTDIR","dev/benchmarks"), "treespace_pool.csv"), row.names = FALSE)

# Union reference per dataset + per-method CID
cidRows <- list()
for (nm in datasets) {
  uni <- do.call(c, lapply(names(methods), function(mn)
    do.call(c, lapply(seq_len(nseed), function(s) allMP[[mn]][[nm]][[s]]))))
  uni <- uni[!vapply(uni, is.null, TRUE)]; class(uni) <- "multiPhylo"
  refCons <- .strict(uni)
  for (mn in names(methods)) {
    ma <- do.call(c, lapply(seq_len(nseed), function(s) allMP[[mn]][[nm]][[s]]))
    ma <- ma[!vapply(ma, is.null, TRUE)]; class(ma) <- "multiPhylo"
    mc <- .strict(ma)
    cid <- tryCatch(as.double(ClusteringInfoDist(mc, refCons, normalize = TRUE)),
                    error = function(e) NA_real_)
    cidRows[[length(cidRows)+1L]] <- data.frame(dataset = nm, method = mn,
      node = mc$Nnode, refNode = refCons$Nnode, cid2union = round(cid, 4))
  }
}
cdf <- do.call(rbind, cidRows)
write.csv(cdf, file.path(Sys.getenv("OUTDIR","dev/benchmarks"), "treespace_pool_cid.csv"), row.names = FALSE)

cat("\n=== median consensus internal nodes (lower = more conservative sampling) ===\n")
print(aggregate(cbind(TNT_node, TSf_p100_node, TSf_pBig_node, TScs6_pBig_node) ~ dataset,
                df, median), row.names = FALSE)
cat("\n=== median MPTs retained ===\n")
print(aggregate(cbind(TNT_n, TSf_p100_n, TSf_pBig_n, TScs6_pBig_n) ~ dataset, df, median), row.names = FALSE)
cat("\n=== median wall (s) ===\n")
print(aggregate(cbind(TNT_wall, TSf_p100_wall, TSf_pBig_wall, TScs6_pBig_wall) ~ dataset, df, median), row.names = FALSE)
cat("\n=== pooled consensus vs union-of-all-methods reference ===\n")
print(cdf[order(cdf$dataset, cdf$method), ], row.names = FALSE)
