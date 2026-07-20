# disc_nstates.R -- vary-n_states DISCRIMINATOR (Mission B, gate #2).
#
# Question (advisor-reframed): does TS's per-candidate TBR scoring cost track
# the (bit-sliced, in n_states words/block) inner loop + gather, or is it
# dominated by n_states-INDEPENDENT fixed per-candidate overhead?
#   flat at corpus n_states (~9) -> fixed-overhead-bound -> packing dead, screen is the lever.
#   proportional to n_states     -> per-word-bound (do NOT auto-build packing:
#                                   TS==TNT bit density; re-check corpus heterogeneity).
#
# Method: hold nChar & nTip fixed, vary n_states. For each config: generate
# random USER data, CONVERGE to a kernel local optimum (loose-cutoff transient),
# then time a NEAR-OPTIMAL round of TBR at that optimum (tight cutoff = the
# regime the ~47 ns fingerprint was measured in). Per-candidate SPR & REROOT ns
# come from the committed TS_IW_TIMING chrono; min-of-runs is the load-robust
# true-cost estimator (box is loaded; load only ever adds time).
#
# Usage: TS_LIB=.agent-disc Rscript dev/profiling/reeval/disc_nstates.R \
#          [nTip=100] [nChar=141] [reps=12] [states=2,4,8,9,16,32]

suppressMessages({
  library(TreeSearch, lib.loc = normalizePath(Sys.getenv("TS_LIB", ".agent-disc"),
                                              winslash = "/"))
  library(TreeTools)
})

args   <- commandArgs(trailingOnly = TRUE)
nTip   <- if (length(args) >= 1) as.integer(args[[1]]) else 100L
nChar  <- if (length(args) >= 2) as.integer(args[[2]]) else 141L
reps   <- if (length(args) >= 3) as.integer(args[[3]]) else 12L
states <- if (length(args) >= 4) as.integer(strsplit(args[[4]], ",")[[1]]) else
          c(2L, 4L, 8L, 9L, 16L, 32L)

makeData <- function(nState, seed) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  m <- matrix(sample(0:(nState - 1L), nTip * nChar, replace = TRUE),
              nrow = nTip, dimnames = list(tips, NULL))
  phy <- phangorn::phyDat(m, type = "USER", levels = as.character(0:(nState - 1L)))
  at  <- attributes(phy)
  list(phy = phy, contrast = at$contrast,
       tip_data = matrix(unlist(phy, use.names = FALSE),
                         nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels,
       nTip = length(phy), labels = names(phy),
       nPat = length(at$weight))
}

runDiag <- function(edge, d) {
  TreeSearch:::ts_tbr_diagnostics(
    edge, d$contrast, d$tip_data, d$weight, d$levels,
    maxHits = 1L, acceptEqual = FALSE, maxChanges = 0L)
}

parseIWT <- function(msg) {
  line <- grep("^IWT ", msg, value = TRUE)
  if (!length(line)) return(NULL)
  line <- line[[length(line)]]
  g <- function(pat) {
    m <- regmatches(line, regexec(pat, line))[[1]]
    if (length(m) >= 2) as.numeric(m[[2]]) else NA_real_
  }
  list(clips   = g("clips=([0-9.]+)"),
       n_spr   = g("SPR n=([0-9.]+)"),
       spr_ns  = g("SPR n=[0-9.]+ ([0-9.]+)ns"),
       n_rer   = g("REROOT n=([0-9.]+)"),
       rer_ns  = g("REROOT n=[0-9.]+ ([0-9.]+)ns"),
       flat    = g("flat=([0-9]+)"))
}

cat(sprintf("=== vary-n_states discriminator: nTip=%d nChar=%d reps=%d ===\n",
            nTip, nChar, nChar))
cat(sprintf("%-8s %-6s %-7s %-9s %-9s %-9s %-9s %-7s\n",
            "nStates", "nPat", "clips", "SPRns", "REROOTns", "totns",
            "SPR/9*", "flat?"))

rows <- list()
for (k in states) {
  d <- makeData(k, seed = 20260714L + k)
  set.seed(100L + k)
  start <- RandomTree(d$phy, root = TRUE)
  edge0 <- Preorder(RenumberTips(start, d$labels))[["edge"]]
  Sys.unsetenv("TS_IW_TIMING")
  conv  <- runDiag(edge0, d)                       # converge (untimed)
  edgeC <- conv$edge

  Sys.setenv(TS_IW_TIMING = "1")
  spr <- rer <- clp <- nsp <- nre <- fl <- rep(NA_real_, reps)
  for (r in seq_len(reps)) {
    msg <- capture.output(.tmp <- runDiag(edgeC, d), type = "message")
    p <- parseIWT(msg)
    if (!is.null(p)) {
      spr[r] <- p$spr_ns; rer[r] <- p$rer_ns; clp[r] <- p$clips
      nsp[r] <- p$n_spr;  nre[r] <- p$n_rer;  fl[r]  <- p$flat
    }
  }
  Sys.unsetenv("TS_IW_TIMING")
  sprM <- suppressWarnings(min(spr, na.rm = TRUE))
  rerM <- suppressWarnings(min(rer, na.rm = TRUE))
  # per-candidate blended cost, weighted by evaluated counts (SPR + REROOT)
  nspM <- stats::median(nsp, na.rm = TRUE); nreM <- stats::median(nre, na.rm = TRUE)
  totM <- (sprM * nspM + rerM * nreM) / (nspM + nreM)
  rows[[as.character(k)]] <- data.frame(
    nStates = k, nPat = d$nPat, clips = stats::median(clp, na.rm = TRUE),
    SPRns = sprM, REROOTns = rerM, totns = totM,
    flat = stats::median(fl, na.rm = TRUE))
  cat(sprintf("%-8d %-6d %-7.0f %-9.2f %-9.2f %-9.2f %-9.2f %-7.0f\n",
              k, d$nPat, stats::median(clp, na.rm = TRUE),
              sprM, rerM, totM, totM / (k / 9), stats::median(fl, na.rm = TRUE)))
}

res <- do.call(rbind, rows)
# Shape test: regress total per-candidate ns on n_states. slope~0 => flat (screen);
# strong positive slope, intercept small => proportional (per-word-bound).
fit <- lm(totns ~ nStates, data = res)
cat(sprintf("\nSHAPE: totns = %.2f + %.3f * nStates   (R^2=%.3f)\n",
            coef(fit)[[1]], coef(fit)[[2]], summary(fit)$r.squared))
cat(sprintf("       ratio totns(nStates=%d)/totns(nStates=%d) = %.2fx  (proportional would be %.2fx)\n",
            max(states), min(states),
            res$totns[which.max(res$nStates)] / res$totns[which.min(res$nStates)],
            max(states) / min(states)))
saveRDS(res, file.path("dev/profiling/reeval",
                       sprintf("disc_nstates_%dt.rds", nTip)))
cat("* SPR/9 col = totns normalised to n_states=9 (flat if ~equal across rows).\n")
