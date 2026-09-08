# B2 collapse-aggressive SPEED test, per-cell SLURM driver (2026-06-23).
#
# Question (the one the correctness gates can't answer): does TS_COLLAPSE_AGGRESSIVE
# reach the known optimum SOONER (wall-clock) than plain search, on the corpus's
# CHAR-POOR regime — where collapse density at/near the optimum is high (zero-length
# branches persist) and search is expensive?  The ≤88-tip morphological roster is
# uniformly character-RICH (density ~0 at the optimum) so it cannot exhibit the
# effect; the mbank training corpus has a substantial char-poor tail (14% of
# training datasets have <1.5 chars/tip) that can.
#
# Design (mirrors fuse_efficiency.R, advisor-reviewed):
# * WALL-CLOCK is the only fair currency; record the anytime (elapsed, best_score)
#   trajectory via progressCallback, paired by seed, ON vs OFF.
# * Both arms share an IDENTICAL recipe; only TS_COLLAPSE_AGGRESSIVE differs.
# * nThreads=1 (wall-clock + collapse dynamics nondeterministic in parallel).
# * gaps `-` -> `?` (corpus protocol; makes matrices NA-free so the aggressive
#   min-length-0 criterion is active rather than falling back to conservative).
# * n>=32 seeds per (key, arm) — the n=1/n=3 timing-noise trap is real (a 12x
#   fuse "speedup" and a "more-evals" collapse signal both vanished at n>=20).
#
# Cell = (arm, key, seed) over expand.grid(B2_ARMS x B2_KEYS x B2_SEEDS), 0-based
# $SLURM_ARRAY_TASK_ID.  One partial RDS per cell -> pull -> b2_speed_analyze.
#
# Env: TS_LIB, B2_KEYS (space-sep catalogue keys), B2_SEEDS (count),
#      B2_BUDGET (per-cell maxSeconds), B2_STRATEGY, MBANK_DIR, MBANK_CAT,
#      PARTIAL_DIR.

suppressMessages({
  ts_lib <- Sys.getenv("TS_LIB", ".agent-fuse")
  if (nzchar(ts_lib) && dir.exists(ts_lib))
    .libPaths(c(normalizePath(ts_lib, winslash = "/"), .libPaths()))
  library(TreeSearch)
  library(TreeTools)
})
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a

mbank_dir <- Sys.getenv("MBANK_DIR", "/nobackup/pjjg18/neotrans/inst/matrices")
mbank_cat <- Sys.getenv("MBANK_CAT", "/nobackup/pjjg18/floor/mbank_catalogue.csv")
budget_s  <- as.double(Sys.getenv("B2_BUDGET", "90"))
strategy  <- Sys.getenv("B2_STRATEGY", "thorough")
keys      <- strsplit(trimws(Sys.getenv("B2_KEYS", "Vinther2008")), "\\s+")[[1]]
nSeeds    <- as.integer(Sys.getenv("B2_SEEDS", "32"))
arms      <- c("off", "on")
pdir      <- Sys.getenv("PARTIAL_DIR", "b2_speed_partials")

# --- loader (mirrors floor_cell.R: catalogue key -> nex -> phyDat -> fitch) ---
load_matrix <- function(key) {
  e <- new.env(); utils::data("inapplicable.phyData", package = "TreeSearch", envir = e)
  dl <- e[["inapplicable.phyData"]]
  if (key %in% names(dl)) return(dl[[key]])
  cat0 <- utils::read.csv(mbank_cat, stringsAsFactors = FALSE)
  r <- cat0[cat0$key == key, , drop = FALSE]
  if (nrow(r) != 1L) stop("key not unique in catalogue: ", key)
  fp <- file.path(mbank_dir, r$filename[1])
  if (!file.exists(fp)) stop("matrix not found: ", fp)
  TreeTools::ReadAsPhyDat(fp)
}
fitch_convert <- function(phy) {
  m <- PhyDatToMatrix(phy, ambigNA = FALSE); m[m == "-"] <- "?"; MatrixToPhyDat(m)
}

# One timed run; records the (elapsed, best) trajectory + final score/wall.
runCell <- function(phy, seed, aggr) {
  if (aggr) Sys.setenv(TS_COLLAPSE_AGGRESSIVE = "1") else Sys.unsetenv("TS_COLLAPSE_AGGRESSIVE")
  on.exit(Sys.unsetenv("TS_COLLAPSE_AGGRESSIVE"), add = TRUE)
  rows <- list()
  cb <- function(info) {
    rows[[length(rows) + 1L]] <<- c(elapsed = as.double(info$elapsed %||% NA),
                                    best = as.double(info$best_score %||% NA))
  }
  set.seed(seed)
  t0 <- Sys.time()
  invisible(suppressWarnings(MaximizeParsimony(
    phy, strategy = strategy, maxReplicates = 5000L, maxSeconds = budget_s,
    nThreads = 1L, verbosity = 0L, progressCallback = cb)))
  wall <- as.double(difftime(Sys.time(), t0, units = "secs"))
  tr <- if (length(rows)) do.call(rbind, rows) else matrix(NA_real_, 1, 2,
            dimnames = list(NULL, c("elapsed", "best")))
  tr <- tr[!is.na(tr[, "elapsed"]) & !is.na(tr[, "best"]), , drop = FALSE]
  ord <- order(tr[, "elapsed"])
  list(arm = if (aggr) "on" else "off", seed = seed, wall = wall,
       elapsed = tr[ord, "elapsed"], best = cummin(tr[ord, "best"]),
       finalBest = if (nrow(tr)) min(tr[, "best"]) else NA_real_)
}

# --- cell dispatch ---
args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(if (length(args)) args[[1]] else Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
grid <- expand.grid(arm = arms, key = keys, seed = seq_len(nSeeds),
                    stringsAsFactors = FALSE)
cell <- grid[idx + 1L, ]
phy <- fitch_convert(load_matrix(cell$key))
res <- runCell(phy, cell$seed, aggr = (cell$arm == "on"))
res$key <- cell$key; res$ntax <- length(phy); res$nchar <- sum(attr(phy, "weight"))
dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(pdir, sprintf("b2_%05d.rds", idx)))
cat(sprintf("[cell %d] %s key=%s seed=%d final=%g wall=%.1f\n",
            idx, res$arm, cell$key, cell$seed, res$finalBest, res$wall))
