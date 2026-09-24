# lever_b_curve.R -- Mission A / Lever B Stage-1: post-first-hit redundancy curve.
#
# QUESTION: after the optimum is FIRST hit, is the remaining wall CONVERTIBLE to
#   new distinct MPTs (curve climbs / bursty => Lever B viable) or FLAT (=> the
#   only breadth-safe fix is stopping earlier, which cuts breadth => Lever B
#   REFUTED for the breadth constraint)?  AND: how much of the post-hit breadth
#   comes from the MAIN LOOP (ratchet/sector/drift/fuse reps) vs the terminal
#   MPT-enumerator (the cheap TBR plateau walk from each pool seed)?
#
# METHOD (NO engine change):
#   * uncapped pool (poolMaxSize huge) + pool_suboptimal 0 + nThreads 1
#       => NO eviction => pool_size == count_at_best EXACTLY, and
#          hits_to_best - pool_size is an EXACT cumulative dedup-reject count.
#     (The naive capped proxy lies here: once the pool caps, eviction pins
#      count_at_best while hits climb -> a fake dedup signal in exactly this
#      regime.  Uncapping removes that confound. It perturbs split_freq/fuse
#      donor feedback, so this is a CHARACTERISATION run, not the production one.)
#   * targetHits huge => never converge-stop => we observe the FULL over-search
#     regime (all maxReplicates), which is the wall Lever B wants to reclaim.
#   * progressCallback captures per-rep "replicate" rows (PRE-enumerator) and the
#     terminal "done" row (POST-enumerator).  The two gaps decompose breadth:
#       (pre_enum  - count@first_hit) = MAIN-LOOP plateau conversion
#       (post_enum - pre_enum)        = ENUMERATOR conversion
#
# Usage:
#   TS_LIB=<lib> Rscript lever_b_curve.R <matrix.nex|named> <ew|iw> <strategy> \
#                                        <maxRep> <seed> <outdir> [poolCap] [maxSec]
suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-disc"),
                                              winslash = "/"))
  library(TreeTools)
})

a        <- commandArgs(trailingOnly = TRUE)
dset     <- a[[1]]
mode     <- if (length(a) >= 2) a[[2]] else "ew"
strategy <- if (length(a) >= 3) a[[3]] else "thorough"
maxRep   <- if (length(a) >= 4) as.integer(a[[4]]) else 200L
seed     <- if (length(a) >= 5) as.integer(a[[5]]) else 1L
outdir   <- if (length(a) >= 6) a[[6]] else "."
poolCap  <- if (length(a) >= 7) as.integer(a[[7]]) else 20000L
maxSec   <- if (length(a) >= 8) as.double(a[[8]]) else 0     # 0 = pure replicate-bound
tHits    <- if (length(a) >= 9) as.integer(a[[9]]) else NA_integer_  # NA = never-stop (observe full)

if (file.exists(dset)) {
  phy0 <- ReadAsPhyDat(dset); nm <- basename(dset)
} else {
  data("inapplicable.phyData", package = "TreeSearch"); phy0 <- inapplicable.phyData[[dset]]
  nm <- dset
}
# EW-Fitch convention: inapplicable -> missing (matches the kernel/reach harnesses).
m <- PhyDatToMatrix(phy0, ambigNA = FALSE); m[m == "-"] <- "?"; phy <- MatrixToPhyDat(m)
conc <- if (mode == "iw") 10 else Inf
nTip <- length(phy); nChar <- sum(attr(phy, "weight"))
cat(sprintf("=== LEVER-B CURVE %s mode=%s strat=%s nTip=%d nChar=%d maxRep=%d seed=%d poolCap=%d maxSec=%g targetHits=%s ===\n",
            nm, mode, strategy, nTip, nChar, maxRep, seed, poolCap, maxSec,
            if (is.na(tHits)) "inf" else as.character(tHits)))

# --- capture the progress stream ---
buf <- new.env(parent = emptyenv()); buf$rows <- list()
cb <- function(info) {
  buf$rows[[length(buf$rows) + 1L]] <- list(
    phase   = as.character(info$phase),
    rep     = as.integer(info$replicate),
    best    = as.double(info$best_score),
    pool    = as.integer(info$pool_size),
    hits    = as.integer(info$hits_to_best),
    elapsed = as.double(info$elapsed))
}

set.seed(seed)
t0 <- Sys.time()
r <- suppressWarnings(MaximizeParsimony(
  dataset       = phy,
  concavity     = conc,
  strategy      = strategy,
  maxReplicates = maxRep,
  targetHits    = if (is.na(tHits)) as.integer(min(.Machine$integer.max, maxRep * 1000L)) else tHits,
  maxSeconds    = maxSec,
  nThreads      = 1L,
  verbosity     = 0L,
  poolMaxSize   = poolCap,          # DOTS override -> survives strategy preset; ~uncapped
  progressCallback = cb))
wall <- as.double(difftime(Sys.time(), t0, units = "secs"))

# --- assemble per-rep frame ---
df <- do.call(rbind, lapply(buf$rows, function(x)
  data.frame(x, stringsAsFactors = FALSE)))
if (is.null(df) || !nrow(df)) { cat("!! no progress rows captured\n"); quit(status = 1) }

finalBest <- min(df$best[is.finite(df$best)])
repRows   <- df[df$phase == "replicate", , drop = FALSE]
doneRows  <- df[df$phase == "done", , drop = FALSE]

# --- derived metrics (uncapped => pool == count_at_best, hits-pool == dedup) ---
hitMask   <- repRows$best <= finalBest + 1e-9
firstHit  <- if (any(hitMask)) min(repRows$rep[hitMask]) else NA_integer_
idxHit    <- if (!is.na(firstHit)) match(firstHit, repRows$rep) else NA_integer_
atHit     <- if (!is.na(idxHit)) repRows$pool[idxHit] else NA_integer_
preEnum   <- if (nrow(repRows)) repRows$pool[nrow(repRows)] else NA_integer_
postEnum  <- if (nrow(doneRows)) doneRows$pool[nrow(doneRows)] else length(r)
totReps   <- if (nrow(repRows)) max(repRows$rep) else NA_integer_
postReps  <- totReps - firstHit

dedupHit  <- if (!is.na(idxHit)) repRows$hits[idxHit] - atHit else NA_real_
dedupPre  <- if (nrow(repRows)) repRows$hits[nrow(repRows)] - preEnum else NA_real_

mainGain  <- preEnum  - atHit      # MPTs the redundant main-loop reps added at plateau
enumGain  <- postEnum - preEnum    # MPTs the terminal enumerator added
perRepGain <- if (!is.na(postReps) && postReps > 0) mainGain / postReps else NA_real_

# per-phase timers (attributes on result, cumulative totals) if surfaced
tk <- function(k) { v <- attr(r, k); if (is.null(v)) NA_real_ else as.double(v) }
phases <- c("wagner_ms","tbr_ms","xss_ms","rss_ms","css_ms","ratchet_ms",
            "drift_ms","final_tbr_ms","fuse_ms")
phMs <- setNames(vapply(phases, tk, numeric(1)), phases)

# --- write outputs ---
tag <- sprintf("%s_%s_%s_r%d_p%d_t%s_s%d", gsub("[^A-Za-z0-9]+", "", nm), mode, strategy,
               maxRep, poolCap, if (is.na(tHits)) "inf" else as.character(tHits), seed)
outcsv <- file.path(outdir, sprintf("leverb_%s.csv", tag))
utils::write.csv(df, outcsv, row.names = FALSE)

# score dynamics: optimum reached, and whether it improved AFTER the first apparent best
imp_reps <- repRows$rep[c(TRUE, diff(repRows$best) < -1e-9) & is.finite(repRows$best)]
stopped_early <- if (!is.na(totReps)) totReps < maxRep else NA
cat(sprintf("OPTIMUM=%.4g  reached@rep=%s (%.0fs)  score-improved@reps={%s}  stopped_early=%s (reps %s/%d)\n",
            finalBest, firstHit, if (!is.na(idxHit)) repRows$elapsed[idxHit] else NA_real_,
            paste(imp_reps, collapse=","), stopped_early, totReps, maxRep))
cat(sprintf("first_hit_rep=%s  count@first_hit=%s  pre_enum=%s  post_enum=%s  final_ntree=%d  total_reps=%s  wall=%.1fs\n",
            firstHit, atHit, preEnum, postEnum, length(r), totReps, wall))
cat(sprintf("BREADTH SPLIT:  main-loop(plateau)=%s over %s post-hit reps  |  enumerator=%s\n",
            mainGain, postReps, enumGain))
cat(sprintf("DEDUP (cumulative, exact): @first_hit=%s  @pre_enum=%s  (plateau delta=%s)\n",
            dedupHit, dedupPre, dedupPre - dedupHit))
cat(sprintf("VERDICT-HINT: main/enum = %s/%s ; plateau new-MPT/rep = %.3f  (~0 => flat/refute; >0 => convertible)\n",
            mainGain, enumGain, perRepGain))
if (any(is.finite(phMs))) {
  tot <- sum(phMs, na.rm = TRUE)
  cat("phase-ms (cumulative totals):\n")
  for (p in phases) if (is.finite(phMs[[p]]))
    cat(sprintf("   %-13s %10.0f ms  (%.1f%%)\n", p, phMs[[p]],
                if (tot > 0) 100 * phMs[[p]] / tot else 0))
}
cat(sprintf("wrote %s\n", outcsv))
