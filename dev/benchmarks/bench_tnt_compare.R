# TreeSearch vs TNT benchmark comparison
#
# Runs the same morphological datasets through both TreeSearch (C++ driven
# search) and TNT 1.6, comparing score quality and wall-clock time.
#
# Usage:
#   source("inst/benchmarks/bench_tnt_compare.R")
#   results <- tnt_compare_smoke()        # quick: 2 datasets, EW, 5s
#   results <- tnt_compare_full()         # all 14 datasets, EW+IW, 10s+30s
#   save_comparison(results)
#
# Prerequisites:
#   - TNT 1.6 at C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe
#   - TreeSearch installed (default library or specify lib_path)

library(TreeSearch)
library(TreeTools)

source("dev/benchmarks/bench_datasets.R")

# ---- Configuration ----

TNT_EXE <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
STAGING_DIR <- ".tnt-bench"

# ---- TNT helpers ----

#' Check that TNT is available
check_tnt <- function() {
  if (!file.exists(TNT_EXE)) {
    stop("TNT not found at ", TNT_EXE,
         "\nInstall TNT 1.6 or update TNT_EXE path.")
  }
  invisible(TRUE)
}

#' Export all benchmark datasets to TNT format
#' @param dataset_names Character vector of dataset names
#' @param dir Directory for staging files
#' @return Named character vector of file paths
export_datasets_tnt <- function(dataset_names = BENCHMARK_NAMES,
                                dir = STAGING_DIR) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- character(length(dataset_names))
  names(paths) <- dataset_names

  for (nm in dataset_names) {
    ds <- TreeSearch::inapplicable.phyData[[nm]]
    if (is.null(ds)) {
      warning("Dataset ", nm, " not found; skipping.")
      next
    }
    path <- file.path(dir, paste0(nm, ".tnt"))
    WriteTntCharacters(ds, filepath = path)
    paths[nm] <- normalizePath(path, winslash = "/")
  }
  paths
}

#' Write a TNT search script (.run file)
#'
#' @param data_file Basename of the .tnt data file (not full path)
#' @param weighting "EW" or "IW"
#' @param concavity Concavity constant for IW (ignored if EW)
#' @param timeout_s Timeout in seconds
#' @param seed Random seed
#' @param hits Convergence criterion (hits to best score)
#' @param reps Maximum replicates
#' @param dir Directory to write the script
#' @return Path to the .run script file
write_tnt_script <- function(data_file, weighting = "EW", concavity = 3,
                             timeout_s = 30, seed = 1, hits = 5, reps = 20,
                             xpiwe = FALSE,
                             dir = STAGING_DIR) {
  hh <- timeout_s %/% 3600
  mm <- (timeout_s %% 3600) %/% 60
  ss <- timeout_s %% 60
  timeout_str <- sprintf("%d:%02d:%02d", hh, mm, ss)

  # piwe must be set BEFORE proc so TNT reads data in IW mode.
  # xpiwe= (no arg) activates extended IW AFTER proc.
  # xpiwe=K with a number causes a "No command!" parse error in TNT 1.6.
  commands <- "mxram 1024;"
  if (weighting == "IW" || xpiwe) {
    commands <- c(commands, sprintf("piwe=%d;", concavity))
  }
  commands <- c(
    commands,
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed)
  )
  if (xpiwe) {
    # xpiwe( enables per-character concavity (Goloboff 2014).
    # Do NOT add ) — it's a separate sub-option that reverts to single-K.
    commands <- c(commands, "xpiwe(;")
  }

  commands <- c(
    commands,
    sprintf("timeout %s;", timeout_str),
    sprintf("xmult=hits %d replic %d;", hits, reps),
    "best;",
    "quit;"
  )

  # TNT parses the script filename as a command line — digits cause it to

  # split the name and try the alphabetic prefix as a command. Use a fixed
  # purely-alphabetic name and overwrite each time.
  script_name <- "tntbench.run"
  script_path <- file.path(dir, script_name)
  writeLines(commands, script_path)
  normalizePath(script_path, winslash = "/")
}

#' Run a TNT search and parse results
#'
#' @param script_path Path to .run script
#' @param dir Working directory (should contain the .tnt data file)
#' @param timeout_s Hard timeout for the system call (slightly longer than
#'   TNT's internal timeout to allow cleanup)
#' @return List with score, n_trees, wall_s, raw_output
run_tnt <- function(script_path, dir = STAGING_DIR, timeout_s = 120) {
  check_tnt()

  script_basename <- basename(script_path)
  abs_dir <- normalizePath(dir, winslash = "/")

  # setwd into the staging directory so TNT finds the data files,
  # then invoke TNT directly with the script name + trailing semicolon.
  old_wd <- setwd(dir)
  on.exit(setwd(old_wd), add = TRUE)

  t0 <- proc.time()
  output <- tryCatch(
    withCallingHandlers(
      system2(TNT_EXE, args = paste0(script_basename, ";"),
              stdout = TRUE, stderr = TRUE, timeout = timeout_s),
      warning = function(w) {
        # system2 warns on timeout; muffle it but keep the return value
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) paste("ERROR:", conditionMessage(e))
  )
  wall_s <- as.double((proc.time() - t0)[3])

  # Parse output — TNT emits non-UTF8 progress bar characters that break
  # regex on Windows. Strip non-ASCII bytes before parsing.
  output <- iconv(output, from = "", to = "UTF-8", sub = "")
  out_text <- paste(output, collapse = "\n")

  score <- NA_real_
  score_match <- regmatches(out_text,
                            regexpr("Best score:\\s+([0-9]+[.][0-9]+)", out_text))
  if (length(score_match) == 1) {
    score <- as.numeric(sub("Best score:\\s+", "", score_match))
  }

  n_trees <- NA_integer_
  trees_match <- regmatches(out_text,
                            regexpr("([0-9]+) trees? retained", out_text))
  if (length(trees_match) == 1) {
    n_trees <- as.integer(sub(" trees? retained", "", trees_match))
  }

  hits <- NA_integer_
  hits_match <- regmatches(out_text,
                           regexpr("Best score hit ([0-9]+) times?", out_text))
  if (length(hits_match) == 1) {
    hits <- as.integer(gsub("[^0-9]", "", hits_match))
  }

  rearrangements <- NA_real_
  rearr_match <- regmatches(
    out_text,
    regexpr("Total rearrangements examined:\\s+([0-9,]+)", out_text)
  )
  if (length(rearr_match) == 1) {
    rearrangements <- as.numeric(gsub("[^0-9]", "",
                                      sub("Total rearrangements examined:\\s+",
                                          "", rearr_match)))
  }

  list(
    score = score,
    n_trees = n_trees,
    hits = hits,
    rearrangements = rearrangements,
    wall_s = wall_s,
    raw_output = output
  )
}

# ---- TreeSearch helpers ----

#' Run a TreeSearch driven search
#'
#' @param ds Prepared dataset (from prepare_ts_data)
#' @param timeout_s Timeout in seconds
#' @param seed RNG seed
#' @param hits Convergence hits
#' @param reps Maximum replicates
#' @return List with score, n_trees, wall_s, replicates, hits
run_treesearch <- function(ds, timeout_s = 30, seed = 1, hits = 5,
                           reps = 20) {
  set.seed(seed)
  t0 <- proc.time()
  result <- tryCatch(
    TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = as.integer(reps),
      targetHits = as.integer(hits),
      maxSeconds = as.double(timeout_s),
      verbosity = 0L,
      nThreads = 1L
    ),
    error = function(e) {
      list(best_score = NA_real_, pool_size = NA_integer_,
           replicates = NA_integer_, hits_to_best = NA_integer_)
    }
  )
  wall_s <- as.double((proc.time() - t0)[3])

  list(
    score = result$best_score,
    n_trees = result$pool_size,
    wall_s = wall_s,
    replicates = result$replicates,
    hits = result$hits_to_best
  )
}

# ---- Comparison runner ----

#' Run a full comparison grid
#'
#' @param dataset_names Datasets to benchmark
#' @param weightings Character vector, subset of c("EW", "IW")
#' @param concavities Concavity values for IW (ignored for EW)
#' @param timeout_s Numeric vector of timeouts to test
#' @param seeds Numeric vector of seeds
#' @param hits Convergence criterion
#' @param reps Maximum replicates per search
#' @param use_fitch If TRUE, apply fitch_mode() to datasets before scoring.
#'   This makes TreeSearch treat inapplicable tokens as missing, matching
#'   TNT's default behavior.  Without this, TreeSearch uses the Brazeau
#'   et al. (2019) algorithm for inapplicable characters, which gives
#'   higher step counts than TNT's Fitch scorer.  See
#'   dev/benchmarks/iw_search_quality_analysis.md for details.
#' @return data.frame with one row per dataset × weighting × timeout × seed
run_comparison <- function(dataset_names = BENCHMARK_NAMES,
                           weightings = "EW",
                           concavities = 3,
                           timeout_s = 30,
                           seeds = 1:3,
                           hits = 5L,
                           reps = 20L,
                           use_fitch = FALSE) {
  check_tnt()

  if (use_fitch) cat("** Fitch mode: inapplicable treated as missing **\n")
  cat("Exporting datasets to TNT format...\n")
  data_files <- export_datasets_tnt(dataset_names)

  # Load datasets; optionally convert to Fitch-mode for fair TNT comparison
  datasets <- list()
  for (nm in dataset_names) {
    ds <- TreeSearch::inapplicable.phyData[[nm]]
    if (is.null(ds)) next
    if (use_fitch) ds <- fitch_mode(ds)
    datasets[[nm]] <- prepare_ts_data(ds)
  }

  # Keep only datasets that exported successfully
  ok <- nchar(data_files) > 0
  dataset_names <- dataset_names[ok]

  rows <- list()
  n_combos <- length(dataset_names) * length(weightings) * length(timeout_s) *
    length(seeds)
  idx <- 0L

  for (nm in dataset_names) {
    ds <- datasets[[nm]]
    data_basename <- paste0(nm, ".tnt")

    for (wt in weightings) {
      concs <- if (wt == "IW") concavities else NA_real_
      for (k in concs) {
        for (tmo in timeout_s) {
          for (seed in seeds) {
            idx <- idx + 1L
            cat(sprintf("[%d/%d] %s %s k=%s t=%ds seed=%d ... ",
                        idx, n_combos, nm, wt,
                        if (is.na(k)) "-" else k,
                        tmo, seed))

            # --- TNT ---
            script <- write_tnt_script(
              data_basename, weighting = wt, concavity = k,
              timeout_s = tmo, seed = seed, hits = hits, reps = reps
            )
            tnt_res <- run_tnt(script, timeout_s = tmo + 30)

            # --- TreeSearch ---
            conc_arg <- if (wt == "IW") k else Inf
            ts_res <- run_treesearch(
              ds, timeout_s = tmo, seed = seed, hits = hits, reps = reps
            )

            cat(sprintf("TNT=%.1f (%.1fs) TS=%.1f (%.1fs)\n",
                        tnt_res$score, tnt_res$wall_s,
                        ts_res$score, ts_res$wall_s))

            rows[[idx]] <- data.frame(
              dataset = nm,
              n_taxa = ds$n_taxa,
              n_chars = sum(ds$weight),
              weighting = wt,
              concavity = if (is.na(k)) NA_real_ else k,
              timeout_s = tmo,
              seed = seed,
              tnt_score = tnt_res$score,
              tnt_n_trees = tnt_res$n_trees,
              tnt_hits = tnt_res$hits,
              tnt_rearrangements = tnt_res$rearrangements,
              tnt_wall_s = round(tnt_res$wall_s, 3),
              ts_score = ts_res$score,
              ts_n_trees = ts_res$n_trees,
              ts_wall_s = round(ts_res$wall_s, 3),
              ts_replicates = ts_res$replicates,
              ts_hits = ts_res$hits,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
  }

  do.call(rbind, rows)
}

# ---- Convenience wrappers ----

#' Quick smoke test: 2 small datasets, EW only, 5s timeout
tnt_compare_smoke <- function() {
  run_comparison(
    dataset_names = c("Vinther2008", "Sansom2010"),
    weightings = "EW",
    timeout_s = 5,
    seeds = 1:2,
    hits = 3L,
    reps = 10L
  )
}

#' Medium comparison: 5 datasets spanning size range, EW, 10s
tnt_compare_medium <- function() {
  run_comparison(
    dataset_names = c("Vinther2008", "Aria2015", "Griswold1999",
                      "Agnarsson2004", "Zhu2013"),
    weightings = "EW",
    timeout_s = 10,
    seeds = 1:3,
    hits = 5L,
    reps = 20L
  )
}

#' Full comparison: all 14 datasets, EW + IW(k=3), 10s + 30s
#' @param use_fitch See [run_comparison()].
tnt_compare_full <- function(use_fitch = FALSE) {
  run_comparison(
    dataset_names = BENCHMARK_NAMES,
    weightings = c("EW", "IW"),
    concavities = 3,
    timeout_s = c(10, 30),
    seeds = 1:3,
    hits = 5L,
    reps = 20L,
    use_fitch = use_fitch
  )
}

# ---- Persistence ----

#' Save comparison results to CSV
save_comparison <- function(results,
                            file = sprintf("inst/benchmarks/tnt_compare_%s.csv",
                                           format(Sys.time(), "%Y%m%d_%H%M"))) {
  write.csv(results, file, row.names = FALSE)
  cat("Results saved to", file, "\n")
  invisible(file)
}
