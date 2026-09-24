# bench_tnt_settings.R
#
# TNT 1.6 Settings Survey: time-to-best-score across search configurations.
#
# METRIC: time-to-target (TTT) — wall-clock seconds for TNT to first reach
# the best known score B for each dataset. Censored (NA) when B not reached
# within TIMEOUT_S. This is TNT-vs-TNT, so relative wall-clock is valid.
#
# MACHINE METADATA (embedded at write time):
#   Hostname : DW-CZC429715G
#   CPU      : 12th Gen Intel(R) Core(TM) i7-12700
#   RAM      : 15.7 GB
#   TNT      : C:/Programs/Phylogeny/tnt/tnt.exe (v1.6, 32-bit)
#   Date     : 2026-06-17
#
# USAGE:
#   source("dev/benchmarks/bench_tnt_settings.R")
#   tnt_settings_validate()      # Zhu2013 x sect+fuse x seed=1 (smoke test)
#   results <- tnt_settings_full()
#   write.csv(results, "dev/benchmarks/tnt_settings_survey.csv", row.names=FALSE)
#
# ENV OVERRIDES:
#   TNT_EXE        path to tnt.exe  (default: C:/Programs/Phylogeny/tnt/tnt.exe)
#   TNT_TIMEOUT    per-run timeout in seconds          (default: 120)
#   TNT_SEEDS      seeds per (config,dataset)          (default: 5)
#   TNT_B_TIMEOUT  Phase-1 per-seed timeout in seconds (default: 300)
#   TNT_B_SEEDS    Phase-1 seeds to find B             (default: 10)
#   TNT_DATASETS   comma-separated dataset names       (default: 6 gap sets)
#
# REQUIRES: TreeSearch, TreeTools (with PhyDatToMatrix, MatrixToPhyDat,
#           WriteTntCharacters, inapplicable.phyData)

library(TreeSearch)
library(TreeTools)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

TNT_EXE    <- Sys.getenv("TNT_EXE",     "C:/Programs/Phylogeny/tnt/tnt.exe")
TIMEOUT_S  <- as.integer(Sys.getenv("TNT_TIMEOUT",   "120"))
N_SEEDS    <- as.integer(Sys.getenv("TNT_SEEDS",     "5"))
B_TIMEOUT  <- as.integer(Sys.getenv("TNT_B_TIMEOUT", "300"))
B_SEEDS    <- as.integer(Sys.getenv("TNT_B_SEEDS",   "10"))
STAGING    <- ".tnt-survey"

MACHINE <- list(
  hostname = "DW-CZC429715G",
  cpu      = "12th Gen Intel(R) Core(TM) i7-12700",
  ram_gb   = 15.7
)

GAP_DATASETS <- c("Wortley2006", "Eklund2004", "Zanol2014",
                  "Zhu2013", "Giles2015", "Dikow2009")

ALL_DATASETS <- c(
  "Longrich2010", "Vinther2008", "Sansom2010", "DeAssis2011",
  "Aria2015",     "Wortley2006", "Griswold1999", "Schulze2007",
  "Eklund2004",   "Agnarsson2004", "Zanol2014", "Zhu2013",
  "Giles2015",    "Dikow2009"
)

# xmult options per config.
# type="single"  -> xmult = <opts> giveupscore B replic 9999;
# type="level"   -> xmult = level N giveupscore B replic 9999;
# type="default" -> xmult = giveupscore B replic 9999;
# TNT 1.6 quirk: `fuse` inside `xmult =` prompts interactively for a count.
# Workaround: omit `fuse`/`nofuse` for configs that WANT fusing (TNT default
# has fuse=1), and use `nofuse` only where fusing must be disabled.
# The default also has drift=5; use `nodrift` to disable it.
CONFIGS <- list(
  "sect-only"    = list(type = "single", opts = "rss css xss nofuse noratchet nodrift"),
  "sect+fuse"    = list(type = "single", opts = "rss css xss noratchet nodrift"),
  "sect+ratchet" = list(type = "single", opts = "rss css xss ratchet 10 nodrift"),
  "sect+drift"   = list(type = "single", opts = "rss css xss drift 10 noratchet"),
  "all"          = list(type = "single", opts = "rss css xss ratchet 10 drift 10"),
  "ratchet-only" = list(type = "single", opts = "norss nocss noxss ratchet 10 nofuse nodrift"),
  "level0"       = list(type = "level",  level = 0L),
  "level1"       = list(type = "level",  level = 1L),
  "level2"       = list(type = "level",  level = 2L),
  "level3"       = list(type = "level",  level = 3L),
  "level4"       = list(type = "level",  level = 4L),
  "level5"       = list(type = "level",  level = 5L),
  "level10"      = list(type = "level",  level = 10L),
  "default"      = list(type = "default")
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

check_tnt <- function() {
  if (!file.exists(TNT_EXE))
    stop("TNT not found at ", TNT_EXE,
         ". Set env TNT_EXE to the correct path.")
}

#' Fitch-mode phyDat: replace inapplicable "-" tokens with "?" (missing)
to_fitch <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE)
  m[m == "-"] <- "?"
  MatrixToPhyDat(m)
}

#' Export one dataset to <dir>/<name>.tnt (Fitch mode).
#' Dispatches to export_nexus_dataset() for MorphoBank scaling sets.
#' Returns list(ntip, nchar)
export_dataset <- function(name, dir = STAGING) {
  if (name %in% names(SCALE_DATASETS)) return(export_nexus_dataset(name, dir))
  phy <- inapplicable.phyData[[name]]
  if (is.null(phy)) stop("Dataset not found: ", name)
  d <- to_fitch(phy)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  WriteTntCharacters(d, file.path(dir, paste0(name, ".tnt")))
  list(ntip = length(d), nchar = sum(attr(d, "weight")))
}

#' Build xmult run-line for a given config and target B.
#' Uses giveupscore B with hits 5 replic 100 (verified: hits<=5 works with
#' giveupscore; hits>=9999 silently disables it in TNT 1.6).
#' TTT is extracted from per-replicate output (BestScore column), which gives
#' exact timing regardless of whether giveupscore or timeout triggered the stop.
xmult_run_line <- function(cfg, B) {
  suffix <- sprintf("giveupscore %g hits 5 replic 100", B)
  switch(cfg$type,
    single  = sprintf("xmult = %s %s;", cfg$opts, suffix),
    level   = sprintf("xmult = level %d %s;", cfg$level, suffix),
    default = sprintf("xmult = %s;", suffix),
    stop("Unknown config type: ", cfg$type)
  )
}

#' Parse per-replicate output to find time of first BestScore <= B.
#' Returns list(ttb_s = seconds or NA, first_rep = integer or NA)
parse_ttt <- function(txt, B) {
  lines <- strsplit(txt, "\n", fixed = TRUE)[[1L]]
  # Per-replicate lines look like:
  #   "57    TBR        2           624         624          0:00:08    358,825,811"
  # Columns: rep  algor  tree  score  best_score  time  rearrangs
  # Interrupted lines may have "------" for score/best fields — skip those.
  pat <- "^\\s*(\\d+)\\s+\\S+\\s+\\d+\\s+(\\d+(?:\\.\\d+)?)\\s+(\\d+(?:\\.\\d+)?)\\s+(\\d+:\\d{2}:\\d{2})"
  for (ln in lines) {
    m <- regmatches(ln, regexec(pat, ln, perl = TRUE))[[1L]]
    if (length(m) < 5L) next
    best_sc <- suppressWarnings(as.numeric(m[4L]))  # m[1]=full, m[2]=rep, m[3]=score, m[4]=best, m[5]=time
    if (is.na(best_sc) || best_sc > B + 1e-6) next
    tp <- as.integer(strsplit(m[5L], ":")[[1L]])
    ttb_s <- tp[1L] * 3600L + tp[2L] * 60L + tp[3L]
    return(list(ttb_s = as.double(ttb_s), first_rep = as.integer(m[2L])))
  }
  list(ttb_s = NA_real_, first_rep = NA_integer_)
}

#' Write a survey script and return its path (filename: tntsurvey.run)
write_survey_script <- function(data_file, cfg, B, seed,
                                timeout_s = TIMEOUT_S,
                                dir = STAGING) {
  if (!is.finite(B)) stop("B must be finite; Phase 1 did not find a valid score.")
  hh <- timeout_s %/% 3600
  mm <- (timeout_s %% 3600) %/% 60
  ss <- timeout_s %% 60
  lines <- c(
    "mxram 1500;",
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed),
    "hold 1000;",
    sprintf("timeout %d:%02d:%02d;", hh, mm, ss),
    xmult_run_line(cfg, B),
    "best;",
    "quit;"
  )
  path <- file.path(dir, "tntsurvey.run")
  writeLines(lines, path)
  normalizePath(path, winslash = "/")
}

#' Write a Phase-1 script (find best score, no giveupscore target)
#' Uses TNT defaults (fuse=1, drift=5, rss+css) — no `fuse` keyword to avoid
#' the interactive prompt; hits 5 replic 10 gives 5 convergence confirmations.
write_phase1_script <- function(data_file, seed,
                                timeout_s = B_TIMEOUT,
                                dir = STAGING) {
  hh <- timeout_s %/% 3600
  mm <- (timeout_s %% 3600) %/% 60
  ss <- timeout_s %% 60
  lines <- c(
    "mxram 1500;",
    sprintf("proc %s;", data_file),
    sprintf("rseed %d;", seed),
    "hold 1000;",
    sprintf("timeout %d:%02d:%02d;", hh, mm, ss),
    "xmult = hits 5 replic 10;",
    "best;",
    "quit;"
  )
  path <- file.path(dir, "tntphaseone.run")
  writeLines(lines, path)
  normalizePath(path, winslash = "/")
}

#' Run a TNT script from dir; return list(score, rearr, wall_s, raw)
run_tnt <- function(script_path, dir = STAGING, hard_timeout_s = NULL) {
  check_tnt()
  if (is.null(hard_timeout_s)) hard_timeout_s <- TIMEOUT_S + 60L

  old_wd <- setwd(normalizePath(dir))
  on.exit(setwd(old_wd), add = TRUE)

  t0  <- Sys.time()
  raw <- tryCatch(
    withCallingHandlers(
      system2(TNT_EXE,
              args    = paste0(basename(script_path), ";"),
              stdout  = TRUE, stderr = TRUE,
              timeout = hard_timeout_s),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) character(0)
  )
  wall_s <- as.double(difftime(Sys.time(), t0, units = "secs"))

  txt <- paste(iconv(raw, from = "", to = "UTF-8", sub = ""), collapse = "\n")

  score <- NA_real_
  m_sc <- regmatches(txt, regexpr(
    "Best score(?:\\s+\\(TBR\\))?:\\s+[0-9]+\\.?[0-9]*", txt, perl = TRUE))
  if (length(m_sc) == 1L)
    score <- as.numeric(sub(".*:\\s+", "", m_sc))

  rearr <- NA_real_
  m_rr <- regmatches(txt, regexpr(
    "Total rearrangements examined:\\s+[0-9,]+", txt, perl = TRUE))
  if (length(m_rr) == 1L)
    rearr <- as.numeric(gsub("[^0-9]", "", sub(".*:\\s+", "", m_rr)))

  list(score = score, rearr = rearr, wall_s = wall_s, raw = txt)
}

# ---------------------------------------------------------------------------
# Phase 1: establish B per dataset
# ---------------------------------------------------------------------------

#' Find best achievable score B for each dataset using the thorough config.
#' Runs b_seeds seeds, each up to b_timeout_s seconds.
#' Returns named numeric vector: dataset -> B.
establish_B <- function(dataset_names = GAP_DATASETS,
                        b_seeds    = seq_len(B_SEEDS),
                        b_timeout  = B_TIMEOUT) {
  dir.create(STAGING, showWarnings = FALSE, recursive = TRUE)
  B_map <- setNames(rep(Inf, length(dataset_names)), dataset_names)

  for (nm in dataset_names) {
    info <- export_dataset(nm)
    cat(sprintf("Phase1 %s (%dt %dc):", nm, info$ntip, info$nchar))
    for (s in b_seeds) {
      script <- write_phase1_script(paste0(nm, ".tnt"), s,
                                    timeout_s = b_timeout)
      res <- run_tnt(script, hard_timeout_s = b_timeout + 60L)
      if (!is.na(res$score)) B_map[[nm]] <- min(B_map[[nm]], res$score)
      cat(sprintf(" %g(%.0fs)", res$score, res$wall_s))
    }
    cat(sprintf("  => B=%g\n", B_map[[nm]]))
  }
  B_map
}

# ---------------------------------------------------------------------------
# Phase 2: TTT per (config, dataset, seed)
# ---------------------------------------------------------------------------

#' Run the full settings survey and return a data frame.
#' If outfile is set, each row is appended to CSV immediately (crash recovery).
run_survey <- function(dataset_names = GAP_DATASETS,
                       B_map,
                       configs   = CONFIGS,
                       seeds     = seq_len(N_SEEDS),
                       timeout_s = TIMEOUT_S,
                       outfile   = NULL) {
  dir.create(STAGING, showWarnings = FALSE, recursive = TRUE)

  total  <- length(configs) * length(dataset_names) * length(seeds)
  idx    <- 0L
  rows   <- vector("list", total)
  wrote_header <- FALSE

  for (cfg_nm in names(configs)) {
    cfg <- configs[[cfg_nm]]
    for (nm in dataset_names) {
      B    <- B_map[[nm]]
      info <- export_dataset(nm)
      data_file <- paste0(nm, ".tnt")

      for (s in seeds) {
        idx <- idx + 1L
        cat(sprintf("[%d/%d] %-14s %-12s seed=%d B=%-6g ",
                    idx, total, cfg_nm, nm, s, B))

        script <- write_survey_script(data_file, cfg, B, s,
                                      timeout_s = timeout_s)
        res <- run_tnt(script, hard_timeout_s = timeout_s + 60L)

        # TTT determination:
        #   1. Primary: parse per-replicate "BestScore" column (exact, 1s resolution)
        #   2. Fallback: process wall_s when giveupscore fired mid-rep ("------")
        #      and per-replicate TTT is NA or 0 but final_score <= B
        reached_B_final <- isTRUE(!is.na(res$score) && res$score <= B + 1e-6)
        ttt     <- parse_ttt(res$raw, B)
        if (!is.na(ttt$ttb_s) && ttt$ttb_s > 0) {
          reached <- TRUE
          ttb     <- ttt$ttb_s
        } else if (reached_B_final) {
          # giveupscore triggered mid-rep; use process wall_s (accurate for <1s runs)
          reached <- TRUE
          ttb     <- round(res$wall_s, 3)
        } else {
          reached <- FALSE
          ttb     <- NA_real_
        }

        cat(sprintf("score=%-6s reached=%-5s ttt=%.1fs\n",
                    if (is.na(res$score)) "NA" else as.character(res$score),
                    reached, if (reached) ttb else res$wall_s))

        rows[[idx]] <- data.frame(
          machine     = MACHINE$hostname,
          cpu         = MACHINE$cpu,
          ram_gb      = MACHINE$ram_gb,
          config      = cfg_nm,
          dataset     = nm,
          ntip        = info$ntip,
          seed        = s,
          B           = B,
          reached_B   = reached,
          wall_s      = ttb,          # NA = censored (did not reach B)
          final_score = res$score,
          rearr       = res$rearr,
          stringsAsFactors = FALSE
        )

        if (!is.null(outfile)) {
          write.table(rows[[idx]], outfile,
                      append = wrote_header, sep = ",",
                      row.names = FALSE, col.names = !wrote_header,
                      quote = TRUE)
          wrote_header <- TRUE
        }
      }
    }
  }
  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Convenience entry points
# ---------------------------------------------------------------------------

#' Smoke test: Zhu2013 x sect+fuse x seed=1
#' Prints script, raw TNT output, and parsed result.
tnt_settings_validate <- function() {
  check_tnt()
  dir.create(STAGING, showWarnings = FALSE, recursive = TRUE)

  nm  <- "Zhu2013"
  cfg <- CONFIGS[["sect+fuse"]]
  cat("=== Validation run: Zhu2013 / sect+fuse / seed=1 ===\n\n")

  cat("--- Phase 1: finding B (3 seeds x 60s) ---\n")
  info <- export_dataset(nm)
  cat(sprintf("Dataset: %d tips, %d chars\n", info$ntip, info$nchar))

  best <- Inf
  for (s in 1:3) {
    sc <- write_phase1_script(paste0(nm, ".tnt"), s, timeout_s = 60L)
    r  <- run_tnt(sc, hard_timeout_s = 90L)
    cat(sprintf("  seed=%d  score=%g  wall=%.1fs\n", s, r$score, r$wall_s))
    if (!is.na(r$score)) best <- min(best, r$score)
  }
  cat(sprintf("  => B = %g\n\n", best))

  cat("--- Phase 2: giveupscore test ---\n")
  sc2 <- write_survey_script(paste0(nm, ".tnt"), cfg, best, seed = 1L,
                             timeout_s = 60L)
  cat("Script contents:\n")
  cat(readLines(file.path(STAGING, "tntsurvey.run")), sep = "\n")
  cat("\n\n")

  r2 <- run_tnt(sc2, hard_timeout_s = 90L)
  ttt <- parse_ttt(r2$raw, best)
  reached_final <- isTRUE(!is.na(r2$score) && r2$score <= best + 1e-6)
  if (!is.na(ttt$ttb_s) && ttt$ttb_s > 0) {
    reached <- TRUE; ttb_show <- ttt$ttb_s
  } else if (reached_final) {
    reached <- TRUE; ttb_show <- r2$wall_s
  } else {
    reached <- FALSE; ttb_show <- NA_real_
  }
  cat(sprintf("Result: final_score=%s  reached_B=%s  proc_wall=%.1fs  TTT=%s\n",
              if (is.na(r2$score)) "NA" else r2$score,
              reached, r2$wall_s,
              if (reached) sprintf("%.1fs", ttb_show) else "CENSORED"))
  cat("\nRaw TNT output:\n")
  cat(r2$raw)
  invisible(list(res = r2, ttt = ttt))
}

#' Full survey: 6 gap datasets, all configs, 5 seeds.
#' Set TNT_DATASETS env var (comma-separated) to override dataset list.
tnt_settings_full <- function(datasets   = NULL,
                               b_timeout  = B_TIMEOUT,
                               b_seeds    = seq_len(B_SEEDS),
                               run_timeout = TIMEOUT_S,
                               run_seeds  = seq_len(N_SEEDS)) {
  check_tnt()

  if (is.null(datasets)) {
    env_ds <- Sys.getenv("TNT_DATASETS", "")
    datasets <- if (nchar(env_ds) > 0)
      trimws(strsplit(env_ds, ",")[[1]])
    else
      GAP_DATASETS
  }

  cat("=== TNT 1.6 Settings Survey ===\n")
  cat(sprintf("Machine : %s\n", MACHINE$hostname))
  cat(sprintf("CPU     : %s\n", MACHINE$cpu))
  cat(sprintf("RAM     : %.1f GB\n", MACHINE$ram_gb))
  cat(sprintf("TNT     : %s\n", TNT_EXE))
  cat(sprintf("Date    : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Datasets: %s\n", paste(datasets, collapse = ", ")))
  cat(sprintf("Configs : %d\n", length(CONFIGS)))
  cat(sprintf("Seeds   : %d\n", length(run_seeds)))
  cat(sprintf("Timeout : %ds per run\n\n", run_timeout))

  B_map <- establish_B(datasets, b_seeds, b_timeout)
  cat("\nB values:\n")
  print(B_map)
  cat("\n")

  outfile <- file.path("dev/benchmarks", "tnt_settings_survey.csv")
  results <- run_survey(datasets, B_map, CONFIGS, run_seeds, run_timeout,
                        outfile = outfile)
  cat(sprintf("\nResults written incrementally to %s\n", outfile))
  results
}

# ---------------------------------------------------------------------------
# Scaling survey: larger MorphoBank datasets (above n=90 sector inflection)
# ---------------------------------------------------------------------------

NEOTRANS_DIR <- normalizePath("../neotrans/inst/projects", winslash = "/",
                               mustWork = FALSE)

SCALE_DATASETS <- list(
  "project691"  = list(path = file.path(NEOTRANS_DIR, "project691.nex")),
  "project4230" = list(path = file.path(NEOTRANS_DIR, "project4230.nex")),
  "project4103" = list(path = file.path(NEOTRANS_DIR, "project4103.nex")),
  "project3763" = list(path = file.path(NEOTRANS_DIR, "project3763.nex"))
)

#' Export a MorphoBank NEXUS dataset to <dir>/<name>.tnt (Fitch mode).
#' Parenthesised polymorphisms (e.g. "(0,1)") are recoded as "?" before reading.
export_nexus_dataset <- function(name, dir = STAGING) {
  info <- SCALE_DATASETS[[name]]
  if (is.null(info)) stop("Unknown scaling dataset: ", name)
  if (!requireNamespace("phangorn", quietly = TRUE))
    stop("phangorn required for NEXUS reading: install.packages('phangorn')")
  lines <- readLines(info$path)
  lines <- gsub("\\([0-9,]+\\)", "?", lines, perl = TRUE)
  tmp <- tempfile(fileext = ".nex")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines, tmp)
  phy <- phangorn::read.phyDat(tmp, format = "nexus", type = "STANDARD")
  d <- to_fitch(phy)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  WriteTntCharacters(d, file.path(dir, paste0(name, ".tnt")))
  list(ntip = length(d), nchar = sum(attr(d, "weight")))
}

#' Scaling survey: 4 MorphoBank datasets (103-205 taxa), all 14 configs, 3 seeds.
#' Uses longer timeouts than the gap-dataset survey (300s/run, 600s Phase-1).
tnt_scaling_full <- function(datasets    = names(SCALE_DATASETS),
                              b_timeout  = 600L,
                              b_seeds    = seq_len(5L),
                              run_timeout = 300L,
                              run_seeds  = seq_len(3L)) {
  check_tnt()
  if (!requireNamespace("phangorn", quietly = TRUE))
    stop("phangorn required: install.packages('phangorn')")

  cat("=== TNT 1.6 Scaling Survey ===\n")
  cat(sprintf("Machine : %s\n", MACHINE$hostname))
  cat(sprintf("CPU     : %s\n", MACHINE$cpu))
  cat(sprintf("RAM     : %.1f GB\n", MACHINE$ram_gb))
  cat(sprintf("TNT     : %s\n", TNT_EXE))
  cat(sprintf("Date    : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Datasets: %s\n", paste(datasets, collapse = ", ")))
  cat(sprintf("Configs : %d\n", length(CONFIGS)))
  cat(sprintf("Seeds   : %d\n", length(run_seeds)))
  cat(sprintf("Timeout : %ds per run\n\n", run_timeout))

  B_map <- establish_B(datasets, b_seeds, b_timeout)
  cat("\nB values:\n")
  print(B_map)
  cat("\n")

  outfile <- file.path("dev/benchmarks", "tnt_scaling_survey.csv")
  results <- run_survey(datasets, B_map, CONFIGS, run_seeds, run_timeout,
                        outfile = outfile)
  cat(sprintf("\nResults written incrementally to %s\n", outfile))
  results
}
