# ratchet_race.R — #39 gate-2 ratchet isolation race: TS ratchet vs TNT ratchet
# from an IDENTICAL Wagner start. Answers "does TNT reach the optimum in fewer
# reweight cycles?" Unit = rearrangements (TS total_tbr_moves <-> TNT "Total
# rearrangements examined"); score-parity = validity gate; seed distributions.
# Both optimise the same Fitch objective (inapplicable -> '?').
#
# Env (tbr_shared_start_lib.R reads TS_LIB / TNT_EXE / T0_DIR): SHARED_LIB (path
#   to tbr_shared_start_lib.R), RACE_DATASETS, RACE_SEEDS, RACE_WAGSEED,
#   RAT_ITER, RACE_OUT.
source(Sys.getenv("SHARED_LIB", "dev/benchmarks/tbr_shared_start_lib.R"))

dsN     <- strsplit(trimws(Sys.getenv("RACE_DATASETS",
             "Wortley2006 Giles2015 Zhu2013 Zanol2014")), "\\s+")[[1]]
seeds   <- as.integer(strsplit(trimws(Sys.getenv("RACE_SEEDS", "1 2 3 4 5")), "\\s+")[[1]])
wagSeed <- as.integer(Sys.getenv("RACE_WAGSEED", "11"))
nIter   <- as.integer(Sys.getenv("RAT_ITER", "30"))

# As GrepNum but takes the LAST match (the ratchet's rearrangement line, after
# any tread/check bookkeeping).
GrepNumLast <- function(out, pat) {
  hit <- grep(pat, out, value = TRUE)
  if (!length(hit)) return(NA_real_)
  suppressWarnings(as.numeric(gsub(",", "", sub(pat, "\\1", hit[length(hit)]))))
}

# TNT ratchet from a tread'd start: `ratchet=iter N` (pinned syntax). Reports
# best score (R-scored from saved trees), total rearrangements examined, wall.
TntRatchet <- function(d, startTree, seed, nIter, hold = 1000) {
  script <- c("mxram 1024;", "taxname=;", "proc data.tnt;",
              paste0("rseed ", seed, ";"),
              paste0("hold ", hold, ";"),
              paste0("tread ", ToTntTree(startTree), ";"),
              paste0("ratchet=iter ", nIter, ";"),
              "tsave *out.tre;", "save;", "tsave/;", "quit;")
  wd <- file.path(tempdir(), paste0("tntrat", Sys.getpid()))
  unlink(wd, recursive = TRUE); dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  WriteTntCharacters(d$phy, file.path(wd, "data.tnt"))
  old <- setwd(wd); on.exit(setwd(old), add = TRUE)
  .t0 <- Sys.time()
  # STDIN pipe (not runfile-arg): headless 64-bit TNT on Hamilton launches the
  # curses UI when handed a runfile arg and yields no parseable stdout.
  out <- suppressWarnings(system2(TNT_EXE, input = script, stdout = TRUE, stderr = TRUE))
  .wall <- as.double(difftime(Sys.time(), .t0, units = "secs"))
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  rearr <- GrepNumLast(out, ".*Total rearrangements examined:\\s*([0-9,]+).*")
  trees <- tryCatch(ReadTntTree(file.path(wd, "out.tre")), error = function(e) NULL)
  finalR <- if (is.null(trees)) NA_real_ else {
    if (inherits(trees, "multiPhylo")) min(vapply(trees, TreeLength, double(1), d$phy))
    else TreeLength(trees, d$phy)
  }
  data.frame(engine = "TNT", seed = seed, final_len = finalR,
             rearrangements = rearr, wall = .wall, stringsAsFactors = FALSE)
}

# TS ratchet from the same start: ts_ratchet_search with production-like params
# (perturbProb 0.25, perturbMaxMoves 5, maxHits 1).
TsRatchet <- function(d, startTree, seed, nIter) {
  edge <- PhyloToKernelEdge(startTree, d)
  set.seed(seed)
  .t0 <- Sys.time()
  res <- TreeSearch:::ts_ratchet_search(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    nCycles = nIter, perturbProb = 0.25, maxHits = 1L,
    perturbMode = 0L, perturbMaxMoves = 5L)
  .wall <- as.double(difftime(Sys.time(), .t0, units = "secs"))
  resTree <- structure(list(edge = res$edge, Nnode = d$nTip - 1L,
                            tip.label = names(d$phy)), class = "phylo")
  data.frame(engine = "TS", seed = seed,
             final_len = TreeLength(resTree, d$phy),
             rearrangements = res$total_tbr_moves, wall = .wall,
             stringsAsFactors = FALSE)
}

cat(sprintf("RATCHET RACE | lib=%s | iter=%d | datasets {%s} | seeds {%s}\n",
            Sys.getenv("TS_LIB"), nIter, paste(dsN, collapse = ","),
            paste(seeds, collapse = ",")))

allRows <- list()
for (nm in dsN) {
  d <- prepareDataset(nm)
  set.seed(wagSeed)
  wag <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight, d$levels)
  wagTree <- Preorder(RenumberTips(structure(list(edge = wag$edge, Nnode = d$nTip - 1L,
               tip.label = names(d$phy)), class = "phylo"), names(d$phy)))
  startLen <- TreeLength(wagTree, d$phy)
  cat(sprintf("\n==== %s (%dt) Wagner(seed %d) start_len=%.0f ====\n",
              nm, d$nTip, wagSeed, startLen))
  rows <- list()
  for (s in seeds) {
    rows[[length(rows) + 1]] <- cbind(dataset = nm, tips = d$nTip,
        start_len = startLen, TntRatchet(d, wagTree, s, nIter))
    rows[[length(rows) + 1]] <- cbind(dataset = nm, tips = d$nTip,
        start_len = startLen, TsRatchet(d, wagTree, s, nIter))
  }
  dd <- do.call(rbind, rows)
  allRows[[nm]] <- dd
  print(dd[, c("engine", "seed", "final_len", "rearrangements", "wall")], row.names = FALSE)
}
res <- do.call(rbind, allRows)

cat(sprintf("\n=== PER-DATASET MEDIAN (ratchet from shared Wagner start, %d iters) ===\n", nIter))
agg <- do.call(rbind, lapply(split(res, res$dataset), function(z) {
  tnt <- z[z$engine == "TNT", ]; ts <- z[z$engine == "TS", ]
  data.frame(dataset = z$dataset[1], tips = z$tips[1], start_len = z$start_len[1],
    ts_final = median(ts$final_len), tnt_final = median(tnt$final_len),
    ts_rearr = round(median(ts$rearrangements)), tnt_rearr = round(median(tnt$rearrangements)),
    rearr_ratio = round(median(ts$rearrangements) / median(tnt$rearrangements), 2),
    ts_wall = round(median(ts$wall), 2), tnt_wall = round(median(tnt$wall), 2),
    wall_ratio = round(median(ts$wall) / median(tnt$wall), 2),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tips), ]
print(agg, row.names = FALSE)
cat("\n--- READ ---\n")
cat("PRIMARY metrics = score@fixed-iters (ts_final vs tnt_final) + wall (64-bit authoritative).\n")
cat("UNIT CAVEAT: TS rearrangements = total_tbr_moves (APPLIED moves) is NOT commensurable with TNT\n")
cat("'Total rearrangements examined' (EXAMINED candidates) -> rearr_ratio is NOT an efficiency ratio.\n")
cat("A clean ratchet-efficiency race needs TS to expose examined-candidates (RatchetResult lacks it).\n")
cat("So this is a COARSE score@budget + wall probe (advisor: order-of-magnitude only).\n")
outCsv <- Sys.getenv("RACE_OUT", "ratchet_race.csv")
write.csv(res, outCsv, row.names = FALSE)
cat(sprintf("\nrows -> %s\n", outCsv))
