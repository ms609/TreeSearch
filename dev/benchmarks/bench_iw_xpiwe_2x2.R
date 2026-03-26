#!/usr/bin/env Rscript
# 2x2 IW/XPIWE benchmark: TreeSearch vs TNT 1.6
# Compares standard IW (extended_iw=FALSE) and XPIWE (extended_iw=TRUE)
# for both engines, using k=3 to match test reference values.

cat("=== IW vs XPIWE benchmark (k=3) ===\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

AGENTLIB <- normalizePath("../../../.builds/TreeSearch-F")
.libPaths(c(AGENTLIB, .libPaths()))
library(TreeSearch.F)
.Internal(registerNamespace("TreeSearch", asNamespace("TreeSearch.F")))
library(TreeTools)
data("inapplicable.phyData", package = "TreeSearch")

TNT_EXE <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
stopifnot("TNT not found" = file.exists(TNT_EXE))
STAGING_DIR <- normalizePath("../../../.tnt-bench", mustWork = FALSE)
dir.create(STAGING_DIR, showWarnings = FALSE)

K <- 3

# --- TNT helpers ---
parse_tnt_score <- function(out_text) {
  m <- regmatches(out_text, regexpr("Best score: ([0-9]+[.][0-9]+)", out_text))
  if (length(m) == 1) return(as.numeric(sub("Best score: ", "", m)))
  NA_real_
}

parse_tnt_trees <- function(out_text) {
  m <- regmatches(out_text, regexpr("([0-9]+) trees? retained", out_text))
  if (length(m) == 1) return(as.integer(sub(" trees? retained", "", m)))
  NA_integer_
}

run_tnt <- function(data_file, timeout_s, seed, hits, reps, k, use_xpiwe) {
  hh <- timeout_s %/% 3600
  mm <- (timeout_s %% 3600) %/% 60
  ss <- timeout_s %% 60
  commands <- c(
    "mxram 1024;",
    sprintf("piwe=%d;", k),
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed)
  )
  if (use_xpiwe) commands <- c(commands, "xpiwe(;")
  commands <- c(commands,
    sprintf("timeout %d:%02d:%02d;", hh, mm, ss),
    sprintf("xmult=hits %d replic %d;", hits, reps),
    "best;", "quit;"
  )

  script_path <- file.path(STAGING_DIR, "tntbench.run")
  writeLines(commands, script_path)

  old_wd <- setwd(STAGING_DIR)
  on.exit(setwd(old_wd), add = TRUE)

  t0 <- proc.time()
  output <- withCallingHandlers(
    system2(TNT_EXE, args = "tntbench.run;",
            stdout = TRUE, stderr = TRUE, timeout = timeout_s + 60),
    warning = function(w) invokeRestart("muffleWarning")
  )
  wall_s <- as.double((proc.time() - t0)[3])

  output <- iconv(output, from = "", to = "UTF-8", sub = "")
  out_text <- paste(output, collapse = "\n")

  list(score = parse_tnt_score(out_text),
       n_trees = parse_tnt_trees(out_text),
       wall_s = wall_s)
}

# --- Parameters ---
DATASETS <- c("Vinther2008", "Sansom2010", "Sano2011", "Longrich2010",
              "Aria2015", "DeAssis2011", "Griswold1999", "Agnarsson2004",
              "Wortley2006", "Schulze2007")
TIMEOUT <- 30
SEED <- 1
HITS <- 5L
REPS <- 20L

# --- Prepare TNT data ---
for (nm in DATASETS) {
  ds <- inapplicable.phyData[[nm]]
  WriteTntCharacters(ds, filepath = file.path(STAGING_DIR, paste0(nm, ".tnt")))
}

# --- Run benchmark ---
rows <- list()
idx <- 0L

for (nm in DATASETS) {
  ds <- inapplicable.phyData[[nm]]
  n <- length(ds); nc <- sum(attr(ds, "weight"))
  data_file <- paste0(nm, ".tnt")

  cat(sprintf("%-20s (%dt, %dc)\n", nm, n, nc))

  # TNT IW (standard)
  tnt_iw <- run_tnt(data_file, TIMEOUT, SEED, HITS, REPS, K, use_xpiwe = FALSE)
  cat(sprintf("  TNT IW:    %10.5f (%.1fs)\n", tnt_iw$score, tnt_iw$wall_s))

  # TNT XPIWE
  tnt_xp <- run_tnt(data_file, TIMEOUT, SEED, HITS, REPS, K, use_xpiwe = TRUE)
  cat(sprintf("  TNT XPIWE: %10.5f (%.1fs)\n", tnt_xp$score, tnt_xp$wall_s))

  # TreeSearch IW (extended_iw = FALSE)
  set.seed(SEED)
  t0 <- proc.time()
  trees_iw <- MaximizeParsimony(ds, concavity = K, verbosity = 0,
                                 extended_iw = FALSE,
                                 targetHits = HITS, maxReplicates = REPS,
                                 maxSeconds = TIMEOUT)
  ts_iw_time <- as.double((proc.time() - t0)[3])
  ts_iw_score <- TreeLength(trees_iw[[1]], ds, concavity = K, extended_iw = FALSE)
  cat(sprintf("  TS IW:     %10.5f (%.1fs)\n", ts_iw_score, ts_iw_time))

  # TreeSearch XPIWE (extended_iw = TRUE)
  set.seed(SEED)
  t0 <- proc.time()
  trees_xp <- MaximizeParsimony(ds, concavity = K, verbosity = 0,
                                 extended_iw = TRUE,
                                 targetHits = HITS, maxReplicates = REPS,
                                 maxSeconds = TIMEOUT)
  ts_xp_time <- as.double((proc.time() - t0)[3])
  ts_xp_score <- TreeLength(trees_xp[[1]], ds, concavity = K, extended_iw = TRUE)
  cat(sprintf("  TS XPIWE:  %10.5f (%.1fs)\n\n", ts_xp_score, ts_xp_time))

  idx <- idx + 1L
  rows[[idx]] <- data.frame(
    dataset = nm, n_taxa = n, n_chars = nc,
    tnt_iw = tnt_iw$score, tnt_xpiwe = tnt_xp$score,
    ts_iw = ts_iw_score, ts_xpiwe = ts_xp_score,
    tnt_iw_time = round(tnt_iw$wall_s, 2),
    tnt_xp_time = round(tnt_xp$wall_s, 2),
    ts_iw_time = round(ts_iw_time, 2),
    ts_xp_time = round(ts_xp_time, 2),
    stringsAsFactors = FALSE
  )
}

results <- do.call(rbind, rows)
OUTFILE <- sprintf("iw_xpiwe_2x2_%s.csv",
                   format(Sys.time(), "%Y%m%d_%H%M"))
write.csv(results, OUTFILE, row.names = FALSE)

cat("\n=== Summary ===\n")
cat(sprintf("%-20s %10s %10s %10s %10s\n",
            "Dataset", "TNT_IW", "TNT_XPIWE", "TS_IW", "TS_XPIWE"))
cat(strrep("-", 65), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-20s %10.5f %10.5f %10.5f %10.5f\n",
              r$dataset, r$tnt_iw, r$tnt_xpiwe, r$ts_iw, r$ts_xpiwe))
}

ok_iw <- !is.na(results$tnt_iw)
ok_xp <- !is.na(results$tnt_xpiwe)
cat("\nIW:    TS better:", sum(results$ts_iw[ok_iw] < results$tnt_iw[ok_iw] - 1e-4),
    "| Tied:", sum(abs(results$ts_iw[ok_iw] - results$tnt_iw[ok_iw]) < 1e-4),
    "| TNT better:", sum(results$tnt_iw[ok_iw] < results$ts_iw[ok_iw] - 1e-4),
    "| TNT NA:", sum(!ok_iw), "\n")
cat("XPIWE: TS better:", sum(results$ts_xpiwe[ok_xp] < results$tnt_xpiwe[ok_xp] - 1e-4),
    "| Tied:", sum(abs(results$ts_xpiwe[ok_xp] - results$tnt_xpiwe[ok_xp]) < 1e-4),
    "| TNT better:", sum(results$tnt_xpiwe[ok_xp] < results$ts_xpiwe[ok_xp] - 1e-4),
    "| TNT NA:", sum(!ok_xp), "\n")

cat("\nResults saved to:", OUTFILE, "\n")
cat("Finished:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
