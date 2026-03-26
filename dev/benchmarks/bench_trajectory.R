#!/usr/bin/env Rscript
# T-251: TNT vs TreeSearch trajectory comparison
#
# Captures per-replicate search trajectories from both engines on the
# datasets where TNT has the largest score advantage. Focuses on:
#   - Score vs wall-clock time
#   - Rearrangements per improvement (TNT) vs phase cost per improvement (TS)
#   - Escape effectiveness (delta from ratchet/drift/sectorial)
#
# Usage:
#   source("dev/benchmarks/bench_trajectory.R")
#   results <- trajectory_compare()          # all 3 gap datasets, 30s
#   results <- trajectory_compare_quick()    # Wortley2006 only, 10s

library(TreeSearch)
library(TreeTools)
library(dplyr)

TNT_EXE <- "C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe"
STAGING_DIR <- ".tnt-bench"
dir.create(STAGING_DIR, showWarnings = FALSE, recursive = TRUE)

# Datasets with largest persistent gaps (from T-249):
# Geisler2001 +21, Zhu2013 +8, Wortley2006 +7, Conrad2008 +5, Zanol2014 +4
GAP_DATASETS <- c("Geisler2001", "Zhu2013", "Wortley2006")

# ---- Data preparation ----

prepare_dataset <- function(name) {
  ds <- inapplicable.phyData[[name]]
  # Convert inapplicable to missing to match TNT's default Fitch scoring
  mat <- PhyDatToMatrix(ds)
  mat[mat == "-"] <- "?"
  ds_clean <- MatrixToPhyDat(mat)

  # Export for TNT
  tnt_path <- file.path(STAGING_DIR, paste0(name, ".tnt"))
  WriteTntCharacters(ds_clean, filepath = tnt_path)

  # Prepare for TreeSearch C++ bridge
  at <- attributes(ds_clean)
  list(
    name = name,
    phyDat = ds_clean,
    contrast = at$contrast,
    tip_data = matrix(unlist(ds_clean, use.names = FALSE),
                      nrow = length(ds_clean), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(ds_clean),
    n_chars = sum(at$weight),
    tnt_file = paste0(name, ".tnt")
  )
}

# ---- TNT trajectory capture ----

run_tnt_trajectory <- function(data_file, timeout_s = 30, seed = 1,
                               hits = 10L, reps = 100L) {
  commands <- c(
    "mxram 1024;",
    sprintf("proc %s;", data_file),
    "hold 10000;",
    sprintf("rseed %d;", seed),
    sprintf("timeout %d:%02d:%02d;",
            timeout_s %/% 3600, (timeout_s %% 3600) %/% 60, timeout_s %% 60),
    sprintf("xmult=hits %d replic %d;", hits, reps),
    "best;",
    "quit;"
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
  parse_tnt_trajectory(output, wall_s)
}

parse_tnt_trajectory <- function(output, wall_s) {
  out_text <- paste(output, collapse = "\n")

  # TNT uses \r for progress bars — split on \r to get individual lines
  raw_text <- paste(output, collapse = "\n")
  all_lines <- unlist(strsplit(raw_text, "[\r\n]+"))
  all_lines <- trimws(all_lines)

  # Parse per-replicate lines:
  # "1     SECT       6           1301        1301         0:00:01    22,678,443"
  # "5     FUSE      20           ------      ------       0:00:04    100,410,686"
  # Score and Best Score fields can be "------"
  rep_pattern <- "(\\d+)\\s+(SECT|FUSE|RATCH|DRIFT|CSS|RAT|RAS|SPR|TBR|FUS)\\s+(\\d+)\\s+(-{2,}|\\d+)\\s+(-{2,}|\\d+)\\s+(\\d+:\\d+:\\d+)\\s+([0-9,]+)"

  reps <- list()
  for (line in all_lines) {
    m <- regmatches(line, gregexpr(rep_pattern, line, perl = TRUE))[[1]]
    for (match in m) {
      parts <- regmatches(match, regexec(rep_pattern, match, perl = TRUE))[[1]]
      if (length(parts) >= 8) {
        time_parts <- as.integer(strsplit(parts[7], ":")[[1]])
        secs <- time_parts[1] * 3600 + time_parts[2] * 60 + time_parts[3]
        reps[[length(reps) + 1]] <- data.frame(
          replicate = as.integer(parts[2]),
          algorithm = parts[3],
          trees = as.integer(parts[4]),
          score = if (grepl("-", parts[5])) NA_integer_ else as.integer(parts[5]),
          best_score = if (grepl("-", parts[6])) NA_integer_ else as.integer(parts[6]),
          time_s = secs,
          rearrangements = as.numeric(gsub(",", "", parts[8])),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # Parse totals (use raw_text which includes all \r-separated content)
  total_rearr <- NA_real_
  m <- regmatches(raw_text, regexpr("Total rearrangements examined:\\s+([0-9,]+)", raw_text))
  if (length(m) == 1) {
    total_rearr <- as.numeric(gsub("[^0-9]", "", sub("Total rearrangements examined:\\s+", "", m)))
  }

  best_score <- NA_real_
  m <- regmatches(raw_text, regexpr("Best score:\\s+([0-9.]+)", raw_text))
  if (length(m) == 1) best_score <- as.numeric(sub("Best score:\\s+", "", m))

  list(
    trajectory = if (length(reps) > 0) do.call(rbind, reps) else NULL,
    total_rearrangements = total_rearr,
    best_score = best_score,
    wall_s = wall_s,
    raw_output = output
  )
}

# ---- TreeSearch trajectory capture ----

run_ts_trajectory <- function(ds, timeout_s = 30, seed = 1,
                              hits = 10L, reps = 100L) {
  # Capture verbosity=2 output by redirecting Rprintf
  set.seed(seed)

  # Use a text connection to capture the C++ Rprintf output
  log_file <- tempfile(fileext = ".txt")
  t0 <- proc.time()

  # Capture C++ Rprintf output via output diversion
  log_con <- file(log_file, open = "wt")
  sink(log_con, type = "output")

  result <- tryCatch(
    TreeSearch:::ts_driven_search(
      ds$contrast, ds$tip_data, ds$weight, ds$levels,
      maxReplicates = as.integer(reps),
      targetHits = as.integer(hits),
      maxSeconds = as.double(timeout_s),
      verbosity = 2L,
      nThreads = 1L,
      # Match current default strategy
      ratchetCycles = 12L,
      ratchetPerturbProb = 0.25,
      driftCycles = 2L,
      nniFirst = TRUE,
      outerCycles = 1L,
      maxOuterResets = 2L,
      adaptiveLevel = TRUE
    ),
    finally = {
      sink(type = "output")
      close(log_con)
    }
  )
  wall_s <- as.double((proc.time() - t0)[3])

  log_lines <- readLines(log_file, warn = FALSE)
  unlink(log_file)

  parse_ts_trajectory(log_lines, result, wall_s)
}

parse_ts_trajectory <- function(log_lines, result, wall_s) {
  # Parse per-replicate, per-phase data from verbosity=2 output
  # Format: "  Phase score: NNNN [NNN ms total]"
  # Replicate headers: "Replicate N/M" or "Replicate N/M (best: N, pool: N, hits: N)"

  phases <- list()
  current_rep <- 0L
  cumulative_ms <- 0

  for (line in log_lines) {
    # Replicate header
    rep_match <- regmatches(line, regexec("Replicate (\\d+)/(\\d+)", line))[[1]]
    if (length(rep_match) >= 2) {
      current_rep <- as.integer(rep_match[2])
      next
    }

    # Phase line: "  Phase score: NNNN [NNN ms]" or "  Phase score: NNNN [NNN ms total]"
    phase_match <- regmatches(
      line,
      regexec("^\\s+(\\S+)\\s+.*score:\\s+(\\d+)\\s+\\[(\\d+\\.?\\d*)\\s+ms", line)
    )[[1]]
    if (length(phase_match) >= 4) {
      phase_name <- sub("_.*", "", phase_match[2])
      score <- as.integer(phase_match[3])
      ms <- as.numeric(phase_match[4])

      phases[[length(phases) + 1]] <- data.frame(
        replicate = current_rep,
        phase = phase_name,
        score = score,
        phase_ms = ms,
        stringsAsFactors = FALSE
      )
      next
    }

    # Wagner line: "  wag_rand+NNI tree score: NNNN [NNN ms]"
    wag_match <- regmatches(
      line,
      regexec("^\\s+wag.*score:\\s+(\\d+)\\s+\\[(\\d+\\.?\\d*)\\s+ms", line)
    )[[1]]
    if (length(wag_match) >= 3) {
      phases[[length(phases) + 1]] <- data.frame(
        replicate = current_rep,
        phase = "Wagner",
        score = as.integer(wag_match[2]),
        phase_ms = as.numeric(wag_match[3]),
        stringsAsFactors = FALSE
      )
      next
    }

    # Outer cycle reset line
    reset_match <- regmatches(
      line,
      regexec("Outer cycle improved.*\\((\\d+) -> (\\d+)\\)", line)
    )[[1]]
    if (length(reset_match) >= 3) {
      phases[[length(phases) + 1]] <- data.frame(
        replicate = current_rep,
        phase = "Reset",
        score = as.integer(reset_match[3]),
        phase_ms = 0,
        stringsAsFactors = FALSE
      )
    }
  }

  trajectory <- if (length(phases) > 0) do.call(rbind, phases) else NULL

  list(
    trajectory = trajectory,
    best_score = result$best_score,
    replicates = result$replicates,
    hits = result$hits_to_best,
    wall_s = wall_s,
    timings = result$timings,
    log_lines = log_lines
  )
}

# ---- Main comparison ----

trajectory_compare <- function(datasets = GAP_DATASETS,
                               timeout_s = 30, seeds = 1:3) {
  results <- list()

  for (nm in datasets) {
    cat(sprintf("\n=== %s ===\n", nm))
    ds <- prepare_dataset(nm)
    cat(sprintf("  %d taxa, %d chars\n", ds$n_taxa, ds$n_chars))

    for (seed in seeds) {
      cat(sprintf("  Seed %d: ", seed))
      key <- paste0(nm, "_s", seed)

      # TNT
      cat("TNT... ")
      tnt <- run_tnt_trajectory(ds$tnt_file, timeout_s = timeout_s,
                                seed = seed, hits = 10L, reps = 100L)
      cat(sprintf("%.0f (%.1fs, %.0fM rearr) | ", tnt$best_score,
                  tnt$wall_s, tnt$total_rearrangements / 1e6))

      # TreeSearch
      cat("TS... ")
      ts <- run_ts_trajectory(ds, timeout_s = timeout_s,
                              seed = seed, hits = 10L, reps = 100L)
      cat(sprintf("%.0f (%.1fs, %d reps)\n", ts$best_score,
                  ts$wall_s, ts$replicates))

      results[[key]] <- list(
        dataset = nm, seed = seed, n_taxa = ds$n_taxa, n_chars = ds$n_chars,
        tnt = tnt, ts = ts
      )
    }
  }

  results
}

trajectory_compare_quick <- function() {
  trajectory_compare(datasets = "Wortley2006", timeout_s = 10, seeds = 1:2)
}

# ---- Analysis helpers ----

summarize_trajectories <- function(results) {
  rows <- list()
  for (key in names(results)) {
    r <- results[[key]]
    tnt <- r$tnt
    ts <- r$ts

    # TNT trajectory summary
    tnt_traj <- tnt$trajectory
    tnt_n_reps <- if (!is.null(tnt_traj)) max(tnt_traj$replicate) else NA
    tnt_rearr_per_s <- if (!is.na(tnt$total_rearrangements) && tnt$wall_s > 0) {
      round(tnt$total_rearrangements / tnt$wall_s / 1e6, 1)
    } else NA

    # TreeSearch trajectory summary
    ts_traj <- ts$trajectory
    ts_n_phases <- if (!is.null(ts_traj)) nrow(ts_traj) else NA

    # Phase cost breakdown (ms)
    tm <- unlist(ts$timings)
    total_ms <- sum(tm)
    ratchet_pct <- round(100 * tm["ratchet_ms"] / total_ms, 1)
    tbr_pct <- round(100 * tm["tbr_ms"] / total_ms, 1)
    drift_pct <- round(100 * tm["drift_ms"] / total_ms, 1)
    xss_pct <- round(100 * tm["xss_ms"] / total_ms, 1)
    css_pct <- round(100 * tm["css_ms"] / total_ms, 1)

    rows[[key]] <- data.frame(
      dataset = r$dataset, seed = r$seed,
      n_taxa = r$n_taxa, n_chars = r$n_chars,
      tnt_score = tnt$best_score,
      tnt_wall_s = round(tnt$wall_s, 2),
      tnt_reps = tnt_n_reps,
      tnt_rearr_M = round(tnt$total_rearrangements / 1e6, 1),
      tnt_rearr_per_s_M = tnt_rearr_per_s,
      ts_score = ts$best_score,
      ts_wall_s = round(ts$wall_s, 2),
      ts_reps = ts$replicates,
      gap = ts$best_score - tnt$best_score,
      ratchet_pct = ratchet_pct, tbr_pct = tbr_pct,
      drift_pct = drift_pct, xss_pct = xss_pct, css_pct = css_pct,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# Extract per-replicate best score trajectory from TreeSearch log
ts_replicate_trajectory <- function(ts_result) {
  traj <- ts_result$trajectory
  if (is.null(traj)) return(NULL)

  # Get final score per replicate (last phase entry per replicate)
  library(dplyr)
  traj |>
    group_by(replicate) |>
    summarise(
      rep_score = last(score),
      total_phase_ms = sum(phase_ms),
      n_phases = n(),
      n_resets = sum(phase == "Reset"),
      .groups = "drop"
    ) |>
    mutate(
      best_so_far = cummin(rep_score),
      improved = rep_score < lag(best_so_far, default = Inf)
    )
}

# Compare escape effectiveness: how often does each perturbation phase
# actually improve the score?
ts_phase_effectiveness <- function(ts_result) {
  traj <- ts_result$trajectory
  if (is.null(traj)) return(NULL)

  # For each replicate, track score before and after each phase
  traj |>
    group_by(replicate) |>
    mutate(
      prev_score = lag(score, default = first(score)),
      delta = prev_score - score,  # positive = improvement
      improved = delta > 0
    ) |>
    ungroup() |>
    filter(phase != "Wagner", phase != "Reset") |>
    group_by(phase) |>
    summarise(
      n = n(),
      n_improved = sum(improved),
      hit_rate = round(mean(improved), 3),
      mean_delta = round(mean(delta[improved]), 1),
      total_ms = sum(phase_ms),
      ms_per_improvement = if (sum(improved) > 0) {
        round(sum(phase_ms) / sum(improved))
      } else NA_real_,
      .groups = "drop"
    ) |>
    arrange(desc(hit_rate))
}
