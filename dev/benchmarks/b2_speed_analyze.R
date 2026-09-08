# Analyse the B2 collapse-aggressive SPEED A/B (run locally after pulling
# b2_speed_partials).  Paired (by key,seed) ON-vs-OFF wall-clock-to-target on the
# corpus char-poor regime + char-rich controls.  Mirrors fuse_efficiency.R.
#
# Verdict logic ([[mission-wallclock-to-optimum]]): the mission metric is
# WALL-CLOCK-to-optimum.  Per (key,seed) the target is the best score reached by
# EITHER arm (a shared, censoring-free target); time-to-target per arm; paired.
# A real speed win = ON reaches the target sooner on a meaningful fraction at
# matched wall-clock, WITHOUT a final-score regression.  Expect controls (char-
# rich) to be null; the question is whether char-poor datasets show a win.
#
# Usage: Rscript dev/benchmarks/b2_speed_analyze.R [partialDir]

args <- commandArgs(trailingOnly = TRUE)
pdir <- if (length(args)) args[[1]] else "dev/benchmarks/b2_speed_pull"
fs <- list.files(pdir, pattern = "^b2_.*\\.rds$", full.names = TRUE)
if (!length(fs)) stop("no partials in ", pdir)
parts <- lapply(fs, readRDS)

key  <- vapply(parts, function(p) p$key %||% NA_character_, character(1))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L) b else a
key  <- vapply(parts, function(p) if (is.null(p$key)) NA_character_ else p$key, character(1))
arm  <- vapply(parts, function(p) p$arm, character(1))
seed <- vapply(parts, function(p) p$seed, numeric(1))
ntax <- vapply(parts, function(p) p$ntax %||% NA, numeric(1))
nchr <- vapply(parts, function(p) p$nchar %||% NA, numeric(1))
fin  <- vapply(parts, function(p) p$finalBest %||% NA, numeric(1))

bestAt <- function(p, t) { ix <- which(p$elapsed <= t); if (!length(ix)) Inf else p$best[ix[length(ix)]] }
timeTo <- function(p, tgt) { hit <- which(p$best <= tgt + 1e-6); if (!length(hit)) NA_real_ else p$elapsed[hit[1]] }

cap <- max(vapply(parts, function(p) if (length(p$elapsed)) max(p$elapsed) else 0, numeric(1)))
ukeys <- unique(key[!is.na(key)])
cat(sprintf("%-42s %4s %5s %6s %7s %7s %7s %6s\n",
            "key","tax","chr","target","medT_OFF","medT_ON","faster%","dScore"))
summ <- list()
for (k in ukeys) {
  ik <- which(key == k)
  seeds <- intersect(seed[ik][arm[ik]=="off"], seed[ik][arm[ik]=="on"])
  tOff <- tOn <- dFin <- numeric(0)
  for (s in seeds) {
    po <- parts[[ik[arm[ik]=="off" & seed[ik]==s][1]]]
    pn <- parts[[ik[arm[ik]=="on"  & seed[ik]==s][1]]]
    if (is.null(po)||is.null(pn)||!length(po$elapsed)||!length(pn$elapsed)) next
    tgt <- max(min(po$best), min(pn$best))         # shared reachable target (worse of the two minima)
    to <- timeTo(po, tgt); tn <- timeTo(pn, tgt)
    cO <- ifelse(is.na(to), cap, to); cN <- ifelse(is.na(tn), cap, tn)
    tOff <- c(tOff, cO); tOn <- c(tOn, cN); dFin <- c(dFin, min(pn$best) - min(po$best))
  }
  if (!length(tOff)) next
  faster <- mean(tOn < tOff)
  summ[[k]] <- data.frame(key=k, tax=ntax[ik][1], chr=nchr[ik][1],
                          medOFF=median(tOff), medON=median(tOn), faster=faster,
                          dScore=mean(dFin), n=length(tOff))
  cat(sprintf("%-42s %4d %5d %6s %7.1f %7.1f %6.0f%% %+6.2f\n",
              substr(k,1,42), as.integer(ntax[ik][1]), as.integer(nchr[ik][1]),
              "pair", median(tOff), median(tOn), 100*faster, mean(dFin)))
}
S <- do.call(rbind, summ)
S$ratio <- S$chr / S$tax
cat("\n=== summary (char-poor ratio<1.5 = test; char-rich = control) ===\n")
poor <- S[S$ratio < 1.5, ]; rich <- S[S$ratio >= 1.5, ]
cat(sprintf("CHAR-POOR (n=%d keys): median speedup(OFF-ON)=%.2fs, mean faster-frac=%.0f%%, mean dScore=%+.2f (>0 ON worse)\n",
            nrow(poor), median(poor$medOFF - poor$medON), 100*mean(poor$faster), mean(poor$dScore)))
cat(sprintf("CHAR-RICH (n=%d keys): median speedup(OFF-ON)=%.2fs, mean faster-frac=%.0f%%, mean dScore=%+.2f\n",
            nrow(rich), median(rich$medOFF - rich$medON), 100*mean(rich$faster), mean(rich$dScore)))
cat("\nREAD: char-poor faster-frac >>50% & speedup>0 & dScore~0 => collapse is a real speed lever on that regime.\n")
cat("     faster-frac ~50% / speedup~0 => no win even where density is high (lever doesn't convert to wall-clock).\n")
saveRDS(S, file.path(pdir, "b2_speed_summary.rds"))
