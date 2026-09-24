# BIT-IDENTITY GATE for the early-abandonment seed change (ts_tbr.cpp:900).
# Runs a fixed, deterministic search (set.seed + fixed replicates + non-binding
# timeout, nThreads=1) across EW / NA / IW x strict / accept_equal, capturing
# score + MPT count + a topology checksum.  Run against the baseline lib then
# the edited lib; the change is behaviour-preserving iff every row matches.
#   TS_LIB=.agent-aband   Rscript dev/benchmarks/gate_abandon.R   # baseline
#   TS_LIB=.agent-aband2  Rscript dev/benchmarks/gate_abandon.R   # edited
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-aband"),
            winslash = "/"))
  library(TreeTools)
})
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }

lib <- Sys.getenv("TS_LIB", ".agent-aband")

# Deterministic topology checksum (no external pkg; identical trees -> identical).
tdig <- function(r) {
  tr <- if (inherits(r, "multiPhylo")) r else structure(list(r), class = "multiPhylo")
  nw <- vapply(tr, function(t) paste(ape::write.tree(t)), character(1))
  s  <- paste(sort(nw), collapse = "|")
  sprintf("%d:%d", length(tr), sum(as.integer(charToRaw(s))))
}

# config: dataset, scoremode (ew/na/iw), extra args (incl. accept_equal route)
cfgs <- list(
  list(id = "Zanol_ew_strict",   ds = "Zanol2014",  mode = "ew", extra = list(strategy = "thorough")),
  list(id = "Zanol_ew_acceq",    ds = "Zanol2014",  mode = "ew", extra = list(strategy = "thorough", sectorAcceptEqual = TRUE, rssRounds = 4L)),
  list(id = "Zanol_na_strict",   ds = "Zanol2014",  mode = "na", extra = list(strategy = "thorough")),
  list(id = "Wortley_ew_strict", ds = "Wortley2006", mode = "ew", extra = list(strategy = "thorough")),
  list(id = "Wortley_iw_k3",     ds = "Wortley2006", mode = "iw", extra = list(strategy = "thorough", concavity = 3))
)
seeds <- c(1L, 2L)

rows <- list()
for (cf in cfgs) {
  raw <- inapplicable.phyData[[cf$ds]]
  phy <- if (cf$mode == "na") raw else fitch(raw)
  for (sd in seeds) {
    set.seed(sd)
    args <- c(list(dataset = phy, maxReplicates = 1L, maxSeconds = 600,
                   nThreads = 1L, verbosity = 0L), cf$extra)
    t0 <- Sys.time()
    r <- suppressWarnings(do.call(MaximizeParsimony, args))
    wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
    rows[[length(rows) + 1]] <- data.frame(
      lib = basename(lib), config = cf$id, seed = sd,
      score = min(as.double(attr(r, "score"))), dig = tdig(r),
      wall = round(wall, 2), stringsAsFactors = FALSE)
    cat(sprintf("%-18s s%d | score=%.0f dig=%s wall=%.1fs\n",
                cf$id, sd, min(as.double(attr(r, "score"))), tdig(r), wall))
  }
}
S <- do.call(rbind, rows)
out <- sprintf("dev/benchmarks/gate_abandon_%s.csv", basename(lib))
write.csv(S, out, row.names = FALSE)
cat(sprintf("\nWrote %s (sum wall %.1fs)\n", out, sum(S$wall)))
