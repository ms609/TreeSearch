# 2x2 shared-start probe: rasStarts {1,3} x sectorAcceptEqual {F,T}.
#
# Tests the advisor's hypothesis: re-solve looked inert because every lateral
# sector move was reverted (revert-unless-strictly-better). accept_equal is the
# minimal relaxation (Goloboff 2014 plateau traversal). Same shared-start design
# as bench_ras_verify.R: our sectorial runs from TNT's OWN T0 (ratchet/drift off).
# gap = TS_sect - TNT_sect, both from identical T0. Lower = closer to TNT.
#
# Env: TS_LIB (default .agent-ratchet), TNT_EXE, TS_DATASETS, TS_SEEDS, TS_KPASS, OUT_CSV.
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-ratchet"),
            winslash = "/"))
  library(TreeTools)
})
TNT <- Sys.getenv("TNT_EXE", "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe")
seeds <- as.integer(strsplit(trimws(Sys.getenv("TS_SEEDS", "1 2 3")), "\\s+")[[1]])
K   <- as.integer(Sys.getenv("TS_KPASS", "8"))
dsN <- strsplit(trimws(Sys.getenv("TS_DATASETS",
         "Wortley2006 Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
out_csv <- Sys.getenv("OUT_CSV", "dev/benchmarks/sect_2x2.csv")
data("inapplicable.phyData", package = "TreeSearch")
fitch <- function(p) { m <- PhyDatToMatrix(p, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m) }
num <- function(x) as.double(gsub(",", "", x))
wd <- file.path(tempdir(), "sect2x2"); dir.create(wd, showWarnings = FALSE, recursive = TRUE)

run_tnt <- function(phy, seed, kpass) {
  WriteTntCharacters(phy, file.path(wd, "data.tnt"))
  script <- c("mxram 1024;", "proc data.tnt;", "hold 1;", sprintf("rseed %d;", seed),
              "taxname=;", "mult=replic 1;", "tsave *t0.tre;", "save;", "tsave/;",
              rep("sectsch=rss;", kpass), "quit;")
  writeLines(script, file.path(wd, "ss.run"))
  old <- setwd(wd); on.exit(setwd(old))
  out <- suppressWarnings(system2(TNT, args = "ss.run;", stdout = TRUE, stderr = TRUE))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  s_sect <- num(sub(".*best score:\\s*([0-9.]+).*", "\\1",
                    grep("Sectorial search \\(RSS\\), best score:", out, value = TRUE)))
  t0 <- tryCatch(ReadTntTree(file.path(wd, "t0.tre")), error = function(e) NULL)
  if (inherits(t0, "multiPhylo")) t0 <- t0[[1]]
  list(t0 = t0, s_sect = if (length(s_sect)) s_sect[length(s_sect)] else NA)
}
run_ts <- function(d, tree, rss, ras, aeq) {
  set.seed(1)
  nt <- length(d)
  smin <- as.integer(round(nt * 0.35)); smax <- as.integer(round(nt * 0.65))
  r <- suppressWarnings(MaximizeParsimony(d, tree = tree, maxReplicates = 1L,
    nThreads = 1L, strategy = "auto", maxSeconds = 0, verbosity = 0L,
    ratchetCycles = 0L, driftCycles = 0L, xssRounds = 0L, cssRounds = 0L,
    wagnerStarts = 1L, fuseInterval = 9999L, sectorMinSize = smin, sectorMaxSize = smax,
    rssRounds = as.integer(rss), rasStarts = as.integer(ras), sectorAcceptEqual = aeq))
  as.double(attr(r, "score"))
}

rows <- list()
for (nm in dsN) {
  phy <- fitch(inapplicable.phyData[[nm]])
  for (sd in seeds) {
    tn <- run_tnt(phy, sd, K)
    if (is.null(tn$t0)) { cat(sprintf("WARN %s s%d: no T0\n", nm, sd)); next }
    start <- TreeLength(tn$t0, phy)
    base <- run_ts(phy, tn$t0, K, 1L, FALSE)  # polish, strict
    aeqo <- run_ts(phy, tn$t0, K, 1L, TRUE)   # polish + accept_equal
    ras  <- run_ts(phy, tn$t0, K, 3L, FALSE)  # re-solve, strict
    both <- run_ts(phy, tn$t0, K, 3L, TRUE)   # re-solve + accept_equal
    rows[[length(rows) + 1]] <- data.frame(dataset = nm, seed = sd, start = start,
      tnt = tn$s_sect, base = base, aeq = aeqo, ras = ras, both = both,
      g_base = base - tn$s_sect, g_aeq = aeqo - tn$s_sect,
      g_ras = ras - tn$s_sect, g_both = both - tn$s_sect, stringsAsFactors = FALSE)
    cat(sprintf("%-11s s%d | TNT=%.0f | base=%.0f aeq=%.0f ras=%.0f both=%.0f | gaps %+.0f/%+.0f/%+.0f/%+.0f\n",
                nm, sd, tn$s_sect, base, aeqo, ras, both,
                base - tn$s_sect, aeqo - tn$s_sect, ras - tn$s_sect, both - tn$s_sect))
  }
}
S <- do.call(rbind, rows)
cat("\n== medians (gap = TS_sect - TNT_sect from identical T0; lower = closer to TNT) ==\n")
agg <- do.call(rbind, lapply(split(S, S$dataset), function(d) data.frame(
  dataset = d$dataset[1], TNT = median(d$tnt),
  g_base = median(d$g_base), g_aeq = median(d$g_aeq),
  g_ras = median(d$g_ras), g_both = median(d$g_both))))
print(agg, row.names = FALSE)
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
write.csv(S, out_csv, row.names = FALSE)
cat(sprintf("\nWrote %s\n", out_csv))
