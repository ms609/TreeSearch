# Hamilton wall-clock comparison: TreeSearch vs TNT 1.6 (64-bit), ONE dataset.
# Representative of what a real (not highly sophisticated) user runs in each engine,
# timed on identical 64-bit hardware. Scores re-computed in R via TreeLength
# (bitness-independent, authoritative); wall-clock is the comparison of interest.
# TNT results are static — cache once.
#
# Env: TS_LIB, TS_DATASET, TNT_EXE, OUTDIR, NSEED (default 3).
#   Requires LD_LIBRARY_PATH=<tnt>/TNT-bin and TERM=xterm in the job env.
.libPaths(c(Sys.getenv("TS_LIB", .libPaths()[1]), .libPaths()))
suppressMessages({
  library(TreeSearch)
  library(TreeTools)
})
nm    <- Sys.getenv("TS_DATASET", "Zanol2014")
TNT   <- Sys.getenv("TNT_EXE")
nseed <- as.integer(Sys.getenv("NSEED", "3"))
outdir<- Sys.getenv("OUTDIR", ".")
target <- c(Wortley2006 = 480, Zanol2014 = 1261, Zhu2013 = 624, Giles2015 = 670)

data("inapplicable.phyData", package = "TreeSearch")
m <- PhyDatToMatrix(inapplicable.phyData[[nm]], ambigNA = FALSE); m[m == "-"] <- "?"
phy <- MatrixToPhyDat(m); tgt <- target[[nm]]
wd <- file.path(tempdir(), "tnt"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)
WriteTntCharacters(phy, file.path(wd, "data.tnt"))

# --- TreeSearch: realistic preset runs, timed to completion ---
run_ts <- function(strat, seed) {
  set.seed(seed)
  t <- system.time(r <- suppressWarnings(MaximizeParsimony(phy, strategy = strat,
         maxSeconds = 600, nThreads = 1L, verbosity = 0L)))
  c(score = min(as.double(attr(r, "score"))), wall = as.double(t["elapsed"]))
}
# --- TNT: representative user configs (verified headless), timed; re-score tree ---
run_tnt <- function(cfg, seed) {
  out <- file.path(wd, "out.tre")
  if (file.exists(out)) file.remove(out)
  cmds <- c("mxram 2048;", "proc data.tnt;", sprintf("rseed %d;", seed),
            paste0(cfg, ";"), "tsave *out.tre;", "save;", "tsave/;", "quit;")
  old <- setwd(wd); on.exit(setwd(old))
  t <- system.time(system2(TNT, input = cmds, stdout = FALSE, stderr = FALSE))
  tr <- tryCatch(ReadTntTree("out.tre"), error = function(e) NULL)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  sc <- if (is.null(tr)) NA_real_ else TreeLength(tr, phy)
  c(score = sc, wall = as.double(t["elapsed"]))
}

configs <- list(
  list(engine = "TreeSearch", config = "default",      fn = function(s) run_ts("default", s)),
  list(engine = "TreeSearch", config = "thorough",     fn = function(s) run_ts("thorough", s)),
  list(engine = "TNT",        config = "mult-basic",   fn = function(s) run_tnt("mult=replic 10", s)),
  list(engine = "TNT",        config = "xmult-default",fn = function(s) run_tnt("xmult", s)),
  list(engine = "TNT",        config = "xmult-level10",fn = function(s) run_tnt("xmult=level 10", s))
)

rows <- list()
for (cf in configs) for (s in seq_len(nseed)) {
  v <- cf$fn(s)
  rows[[length(rows) + 1L]] <- data.frame(dataset = nm, target = tgt,
    engine = cf$engine, config = cf$config, seed = s,
    score = unname(v["score"]), over = unname(v["score"]) - tgt,
    wall_s = round(unname(v["wall"]), 1))
  cat(sprintf("%-12s %-13s s%d -> %.0f (%+.0f)  [%.1fs]\n",
              cf$engine, cf$config, s, v["score"], v["score"] - tgt, v["wall"]))
}
df <- do.call(rbind, rows)
write.csv(df, file.path(outdir, paste0("timing_", nm, ".csv")), row.names = FALSE)
cat("\n=== median by engine/config ===\n")
agg <- aggregate(cbind(over, wall_s) ~ engine + config, df, median)
print(agg[order(agg$wall_s), ], row.names = FALSE)
