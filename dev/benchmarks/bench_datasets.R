# Benchmark dataset loading and scoring utilities
#
# Usage:
#   source("dev/benchmarks/bench_datasets.R")
#   datasets <- load_benchmark_datasets()
#   run_benchmark_suite(maxSeconds = 30, replicates = 5)

library(TreeSearch)
library(TreeTools)

# The 14 standard benchmark datasets (<=88 tips), ordered by tip count
BENCHMARK_NAMES <- c(
 "Longrich2010",   # 20 tips, 3 states, 45% missing
  "Vinther2008",    # 23 tips, 4 states, 21% missing
  "Sansom2010",     # 23 tips, 4 states, 40% missing
  "DeAssis2011",    # 33 tips, 3 states, 21% inapp
  "Aria2015",       # 35 tips, 6 states, 13% missing
  "Wortley2006",    # 37 tips, 8 states, 31% missing
  "Griswold1999",   # 43 tips, 6 states, 6% missing
  "Schulze2007",    # 52 tips, 3 states, 17% inapp
  "Eklund2004",     # 54 tips, 6 states, 30% missing
  "Agnarsson2004",  # 62 tips, 7 states, 6% missing
  "Zanol2014",      # 74 tips, 9 states, 17% inapp
  "Zhu2013",        # 75 tips, 4 states, 43% missing
  "Giles2015",      # 78 tips, 4 states, 42% missing
  "Dikow2009"       # 88 tips, 9 states, 0.4% missing
)

# Large-tree benchmark datasets (>= 100 tips).
# Loaded from dev/benchmarks/ rather than inapplicable.phyData.
LARGE_BENCHMARK_NAMES <- c(
  "mbank_X30754"    # 180 tips, 425 chars, 40% missing, 20% inapplicable
)

#' Convert a phyDat to Fitch-mode (inapplicable treated as missing)
#'
#' Modifies the contrast matrix so inapplicable tokens ("-") behave as
#' fully ambiguous ("?").  This makes TreeSearch's Fitch scorer match
#' TNT's default behavior, enabling fair score comparisons.
#'
#' Background: TreeSearch's CharacterLength/TreeLength use the Brazeau
#' et al. (2019) algorithm for inapplicable characters, which gives
#' higher (more correct) step counts than Fitch.  TNT always uses Fitch
#' (treating "-" as missing).  This function bridges the gap for
#' benchmarking purposes only.
#'
#' @param dataset A phyDat object
#' @return A phyDat with modified contrast (inapplicable rows set to all-1)
fitch_mode <- function(dataset) {
  contrast <- attr(dataset, "contrast")
  levels <- attr(dataset, "levels")
  inapp_col <- match("-", levels)
  if (is.na(inapp_col)) return(dataset)
  for (i in seq_len(nrow(contrast))) {
    if (contrast[i, inapp_col] == 1 && sum(contrast[i, ]) == 1) {
      contrast[i, ] <- 1
    }
  }
  attr(dataset, "contrast") <- contrast
  dataset
}

#' Prepare raw data for C++ bridge from a phyDat object
#' @param dataset A phyDat object
#' @return List with contrast, tip_data, weight, levels
prepare_ts_data <- function(dataset) {
  at <- attributes(dataset)
  list(
    contrast = at$contrast,
    tip_data = matrix(unlist(dataset, use.names = FALSE),
                      nrow = length(dataset), byrow = TRUE),
    weight = at$weight,
    levels = at$levels,
    n_taxa = length(dataset)
  )
}

#' Load all 14 standard benchmark datasets
#' @return Named list of prepared datasets (ready for C++ bridge)
load_benchmark_datasets <- function() {
  datasets <- list()
  for (nm in BENCHMARK_NAMES) {
    ds <- TreeSearch::inapplicable.phyData[[nm]]
    if (is.null(ds)) {
      warning("Dataset ", nm, " not found in inapplicable.phyData")
      next
    }
    datasets[[nm]] <- prepare_ts_data(ds)
  }
  datasets
}

#' Load large-tree benchmark datasets from dev/benchmarks/
#' @return Named list of prepared datasets (ready for C++ bridge)
load_large_benchmark_datasets <- function() {
  bench_dir <- "dev/benchmarks"
  datasets <- list()
  for (nm in LARGE_BENCHMARK_NAMES) {
    nex_path <- file.path(bench_dir, paste0(nm, ".nex"))
    if (!file.exists(nex_path)) {
      warning("Large dataset file not found: ", nex_path)
      next
    }
    phyDat <- TreeTools::ReadAsPhyDat(nex_path)
    datasets[[nm]] <- prepare_ts_data(phyDat)
  }
  datasets
}

#' Load all benchmark datasets (standard + large)
#' @return Named list of prepared datasets
load_all_benchmark_datasets <- function() {
  c(load_benchmark_datasets(), load_large_benchmark_datasets())
}

#' Characterize a benchmark dataset
#' @param ds Prepared dataset (from prepare_ts_data)
#' @return Data frame with one row of characteristics
characterize_dataset <- function(ds) {
  n_taxa <- ds$n_taxa
  n_patterns <- length(ds$weight)
  n_chars <- sum(ds$weight)
  lvls <- ds$levels
  contrast <- ds$contrast
  n_states <- ncol(contrast)
  inapp_idx <- which(lvls == "-")
  n_app_states <- n_states - length(inapp_idx)

  td <- ds$tip_data
  total_cells <- n_taxa * n_patterns

  n_inapp <- 0L
  n_missing <- 0L
  has_inapp <- length(inapp_idx) > 0
  for (i in seq_len(nrow(contrast))) {
    is_inapp <- has_inapp && contrast[i, inapp_idx] > 0.5
    cols_check <- setdiff(seq_len(n_states), inapp_idx)
    is_all <- all(contrast[i, cols_check] > 0.5)
    count <- sum(td == i)
    if (is_inapp && !is_all) n_inapp <- n_inapp + count
    if (is_all) n_missing <- n_missing + count
  }

  data.frame(
    n_taxa = n_taxa,
    n_chars = n_chars,
    n_patterns = n_patterns,
    pct_inapp = round(100 * n_inapp / total_cells, 1),
    n_app_states = n_app_states,
    pct_missing = round(100 * n_missing / total_cells, 1)
  )
}

#' Run a single benchmark: driven search on one dataset
#' @param name Dataset name (from BENCHMARK_NAMES)
#' @param maxSeconds Timeout in seconds
#' @param maxReplicates Maximum replicates
#' @param seed RNG seed
#' @param datasets Pre-loaded datasets (optional)
#' @return List with score, replicates, time, etc.
score_dataset <- function(name, maxSeconds = 10, maxReplicates = 20L,
                          seed = 42L, datasets = NULL) {
  if (is.null(datasets)) {
    ds <- prepare_ts_data(TreeSearch::inapplicable.phyData[[name]])
  } else {
    ds <- datasets[[name]]
  }
  if (is.null(ds)) stop("Dataset '", name, "' not found")

  set.seed(seed)
  t0 <- proc.time()
  result <- TreeSearch:::ts_driven_search(
    ds$contrast, ds$tip_data, ds$weight, ds$levels,
    maxReplicates = maxReplicates,
    targetHits = 5L,
    ratchetCycles = 5L,
    xssRounds = 1L,
    xssPartitions = 3L,
    fuseInterval = 5L,
    maxSeconds = maxSeconds,
    verbosity = 0L
  )
  elapsed <- (proc.time() - t0)[3]

  list(
    dataset = name,
    n_taxa = ds$n_taxa,
    best_score = result$best_score,
    replicates = result$replicates,
    pool_size = result$pool_size,
    hits_to_best = result$hits_to_best,
    timed_out = result$timed_out,
    elapsed = elapsed
  )
}

#' Run the full benchmark suite
#' @param maxSeconds Timeout per dataset
#' @param replicates Number of independent runs per dataset
#' @param seed Base seed (incremented per replicate)
#' @return Data frame with results
run_benchmark_suite <- function(maxSeconds = 30, replicates = 3L,
                                 seed = 42L) {
  datasets <- load_benchmark_datasets()
  results <- list()

  for (nm in names(datasets)) {
    for (rep in seq_len(replicates)) {
      cat(sprintf("[%s] rep %d/%d (timeout=%ds)...",
                  nm, rep, replicates, maxSeconds))
      res <- score_dataset(nm, maxSeconds = maxSeconds,
                            seed = seed + rep - 1L,
                            datasets = datasets)
      cat(sprintf(" score=%.0f reps=%d time=%.1fs\n",
                  res$best_score, res$replicates, res$elapsed))
      res$replicate <- rep
      results <- c(results, list(as.data.frame(res)))
    }
  }

  do.call(rbind, results)
}

#' Summarize benchmark results
#' @param results Data frame from run_benchmark_suite
#' @return Summary data frame (best score, median time, etc.)
summarize_benchmark <- function(results) {
  datasets <- unique(results$dataset)
  summaries <- list()

  for (nm in datasets) {
    sub <- results[results$dataset == nm, ]
    summaries <- c(summaries, list(data.frame(
      dataset = nm,
      n_taxa = sub$n_taxa[1],
      best_score = min(sub$best_score),
      median_score = median(sub$best_score),
      median_time = round(median(sub$elapsed), 2),
      median_reps = median(sub$replicates),
      stringsAsFactors = FALSE
    )))
  }

  do.call(rbind, summaries)
}

# ===========================================================================
# MorphoBank external benchmark datasets (neotrans corpus)
# ===========================================================================

# Hard-coded path to the neotrans matrices directory.
# The neotrans repo is a sibling of the TreeSearch source tree under GitHub/.
# This is a git submodule, so the path is stable.
NEOTRANS_MATRICES_DIR <- local({
  # Try from TreeSearch source root (getwd() == TreeSearch-a/)
  candidates <- c(
    file.path(getwd(), "..", "neotrans", "inst", "matrices"),
    # From dev/benchmarks/ (when sourcing directly)
    file.path(getwd(), "..", "..", "neotrans", "inst", "matrices")
  )
  for (d in candidates) {
    d_norm <- normalizePath(d, mustWork = FALSE)
    if (dir.exists(d_norm)) return(d_norm)
  }
  # Return the most likely path even if it doesn't exist yet
  normalizePath(candidates[1], mustWork = FALSE)
})

# Minimum taxon count for benchmarking. Matrices below this size are
# trivially solved in milliseconds and contribute no useful signal.
MBANK_MIN_NTAX <- 20L

# Fixed 25-matrix training sample, selected for diversity across size tiers.
# Chosen via max-min distance on standardized (ntax, nchar, pct_missing,
# pct_inapp) within each tier: 7 small, 7 medium, 7 large, 4 xlarge.
# Do not modify: results are only comparable when the same sample is used.
# Used for Fitch-track benchmarking (all scoring modes, including matrices
# with zero inapplicable coding).
MBANK_FIXED_SAMPLE <- c(
  # Small (20-30 taxa)
  "project532", "project2346", "project2451", "project4501",
  "project944", "project971_(1)", "project2762",
  # Medium (31-60 taxa)
  "project826", "project561", "project571", "project4146_(3)",
  "project3688", "project4049", "project423",
  # Large (61-120 taxa)
  "project4286", "project4359", "project4397", "project2084_(1)",
  "project2771", "project2184", "project3938",
  # XLarge (121+ taxa)
  "syab07201", "project4133", "project804", "project4284"
)

# Fixed 20-matrix Brazeau-track sample: training matrices with pct_inapp >= 4%
# (where the Brazeau et al. 2019 three-pass algorithm materially differs from
# Fitch scoring). Selected via max-min distance on standardized (ntax, nchar,
# pct_inapp) within each tier: 5 small, 6 medium, 6 large, 3 xlarge.
# Matrices with < 20 characters excluded (pathological for search benchmarking).
# Do not modify: results are only comparable when the same sample is used.
MBANK_BRAZEAU_SAMPLE <- c(
  # Small (20-30 taxa)
  "project4182", "project2346", "project906", "project4112", "project537",
  # Medium (31-60 taxa)
  "project709", "project561", "project2359", "project4761",
  "project4146_(3)", "project4867",
  # Large (61-120 taxa)
  "project4286", "project3512_(2)", "project2084_(1)", "project2086",
  "project2771", "project3938",
  # XLarge (121+ taxa)
  "project4103", "project804", "syab07204"
)

# Minimum inapplicable percentage for Brazeau-track benchmarking.
# Below this threshold, Brazeau and Fitch scoring produce near-identical
# results, so Brazeau-specific benchmarking adds no information.
MBANK_BRAZEAU_INAPP_THRESHOLD <- 4.0

#' Load the MorphoBank matrix catalogue
#'
#' Reads the pre-built catalogue CSV from dev/benchmarks/mbank_catalogue.csv.
#' Filters to usable matrices (parse_ok, ntax >= MBANK_MIN_NTAX) and
#' optionally excludes redundant multi-matrix duplicates.
#'
#' @param include_redundant If FALSE (default), exclude rows flagged
#'   as redundant in the catalogue.
#' @return Data frame with one row per matrix.
load_mbank_catalogue <- function(include_redundant = FALSE) {
  # Find the catalogue CSV
  cat_candidates <- c(
    file.path(getwd(), "dev", "benchmarks", "mbank_catalogue.csv"),
    file.path(getwd(), "mbank_catalogue.csv")
  )
  cat_path <- NULL
  for (p in cat_candidates) {
    if (file.exists(p)) { cat_path <- p; break }
  }
  if (is.null(cat_path)) {
    stop("mbank_catalogue.csv not found. Run build_mbank_catalogue.R first.")
  }

  cat <- read.csv(cat_path, stringsAsFactors = FALSE)

  # Filter to usable matrices
  cat <- cat[cat$parse_ok & !is.na(cat$ntax) & cat$ntax >= MBANK_MIN_NTAX, ]

  # Exclude redundant multi-matrix duplicates (if column exists)
  if (!include_redundant && "dedup_drop" %in% names(cat)) {
    cat <- cat[!cat$dedup_drop, ]
  }

  # Add tier classification
  cat$tier <- cut(cat$ntax,
                  breaks = c(0, 30, 60, 120, Inf),
                  labels = c("small", "medium", "large", "xlarge"))

  rownames(cat) <- cat$key
  cat
}

#' Load MorphoBank datasets by key
#'
#' Reads .nex files from the neotrans matrices directory and prepares them
#' for the C++ bridge.
#'
#' @param catalogue Data frame from load_mbank_catalogue().
#' @param keys Character vector of matrix keys to load.
#' @param verbose If TRUE, print progress.
#' @return Named list of prepared datasets.
load_mbank_datasets <- function(catalogue, keys, verbose = TRUE) {
  if (!dir.exists(NEOTRANS_MATRICES_DIR)) {
    stop("Neotrans matrices directory not found: ", NEOTRANS_MATRICES_DIR,
         "\nIs the neotrans repo checked out?")
  }

  datasets <- list()
  for (k in keys) {
    if (!k %in% catalogue$key) {
      warning("Key '", k, "' not in catalogue; skipping.")
      next
    }
    row <- catalogue[catalogue$key == k, ]
    nex_path <- file.path(NEOTRANS_MATRICES_DIR, row$filename)
    if (!file.exists(nex_path)) {
      warning("File not found: ", nex_path, "; skipping.")
      next
    }
    if (verbose) {
      cat(sprintf("  Loading %s (%d taxa, %d chars)...\n",
                  k, row$ntax, row$nchar))
    }
    tryCatch({
      pd <- suppressWarnings(TreeTools::ReadAsPhyDat(nex_path))
      datasets[[k]] <- prepare_ts_data(pd)
    }, error = function(e) {
      warning("Failed to load ", k, ": ", conditionMessage(e))
    })
  }
  datasets
}

#' Load a stratified sample of MorphoBank datasets
#'
#' Draws a reproducible stratified sample from the training or validation
#' split, with equal representation from each size tier.
#'
#' @param catalogue Data frame from load_mbank_catalogue().
#' @param n Total number of matrices to sample (approximately).
#' @param seed RNG seed for reproducibility.
#' @param split "training" (default) or "validation".
#' @param tier Optional: restrict to a specific tier ("small", "medium",
#'   "large", "xlarge").
#' @param verbose If TRUE, print summary of what was loaded.
#' @return Named list of prepared datasets.
load_mbank_sample <- function(catalogue, n = 25L, seed = 7193L,
                              split = "training", tier = NULL,
                              verbose = TRUE) {
  pool <- catalogue[catalogue$split == split, ]
  if (!is.null(tier)) {
    pool <- pool[pool$tier == tier, ]
  }
  if (nrow(pool) == 0) {
    stop("No matrices in the ", split, " split",
         if (!is.null(tier)) paste0(" (tier: ", tier, ")") else "")
  }

  # Stratified sampling: allocate n proportionally across tiers
  tier_counts <- table(pool$tier)
  tier_counts <- tier_counts[tier_counts > 0]
  n_per_tier <- round(n * tier_counts / sum(tier_counts))
  # Ensure at least 1 per tier if tier has matrices
  n_per_tier <- pmax(n_per_tier, 1L)

  set.seed(seed)
  selected <- character(0)
  for (t in names(n_per_tier)) {
    tier_pool <- pool[pool$tier == t, ]
    k <- min(n_per_tier[t], nrow(tier_pool))
    selected <- c(selected, sample(tier_pool$key, k))
  }

  if (verbose) {
    cat(sprintf("MorphoBank %s sample: %d matrices from %d tiers\n",
                split, length(selected), length(n_per_tier)))
    for (t in names(n_per_tier)) {
      cat(sprintf("  %s: %d selected (of %d available)\n",
                  t, sum(pool$tier[pool$key %in% selected] == t),
                  sum(pool$tier == t)))
    }
  }

  load_mbank_datasets(catalogue, selected, verbose = verbose)
}

#' Load all MorphoBank datasets for a given split
#'
#' @param catalogue Data frame from load_mbank_catalogue().
#' @param split "training" or "validation".
#' @param verbose If TRUE, print progress.
#' @return Named list of prepared datasets.
load_mbank_split <- function(catalogue, split = "training", verbose = TRUE) {
  pool <- catalogue[catalogue$split == split, ]
  if (verbose) {
    cat(sprintf("Loading all %d %s matrices...\n", nrow(pool), split))
  }
  load_mbank_datasets(catalogue, pool$key, verbose = verbose)
}

# ===========================================================================
# Brazeau-track benchmark utilities
# ===========================================================================

#' Filter catalogue to matrices with meaningful inapplicable coding
#'
#' Returns rows where Brazeau et al. (2019) scoring produces materially
#' different results from Fitch.  Use for Brazeau-track benchmarking.
#'
#' @param catalogue Data frame from load_mbank_catalogue().
#' @param threshold Minimum pct_inapp to include (default:
#'   MBANK_BRAZEAU_INAPP_THRESHOLD = 4%).
#' @param min_chars Minimum character count (default: 20). Matrices below
#'   this are pathological for search benchmarking.
#' @return Filtered data frame.
has_meaningful_inapp <- function(catalogue,
                                 threshold = MBANK_BRAZEAU_INAPP_THRESHOLD,
                                 min_chars = 20L) {
  catalogue[catalogue$pct_inapp >= threshold &
              catalogue$nchar >= min_chars, ]
}

#' Load the fixed 20-matrix Brazeau-track sample
#'
#' These are training matrices with pct_inapp >= 4%, selected for diversity
#' via max-min distance on (ntax, nchar, pct_inapp).  Use for benchmarking
#' search strategies under Brazeau scoring.
#'
#' @param catalogue Data frame from load_mbank_catalogue().
#' @param verbose If TRUE, print summary.
#' @return Named list of prepared datasets.
load_mbank_brazeau_sample <- function(catalogue, verbose = TRUE) {
  if (verbose) {
    cat(sprintf("Brazeau-track sample: %d fixed matrices (pct_inapp >= %.0f%%)\n",
                length(MBANK_BRAZEAU_SAMPLE), MBANK_BRAZEAU_INAPP_THRESHOLD))
  }
  load_mbank_datasets(catalogue, MBANK_BRAZEAU_SAMPLE, verbose = verbose)
}
