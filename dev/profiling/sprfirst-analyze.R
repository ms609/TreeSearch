#!/usr/bin/env Rscript
# Aggregate the sprFirst handoff gate (Job 1).  Reads sprfirst_*.csv from a
# directory (arg 1, default ./sprfirst_out) and prints, per key sorted by tip
# count, the mean end-of-first-TBR score & wall for each arm, plus the two
# decision deltas: spr_ex - nni (does exact SPR warmup beat NNI warmup?) and
# spr_ex - spr_un (did the exactness fix change the warmup at all?).
#
# GO test (advisor): on the LARGE end, spr_ex must reach a strictly better
# handoff score than nni (at comparable or lower wall) for there to be anything
# worth escalating to the full-pipeline frontier.  Washed on the large end -> the
# outer cycles wash it too -> NO-GO.
dir <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "sprfirst_out"
fs <- Sys.glob(file.path(dir, "sprfirst_*.csv"))
if (!length(fs)) stop("no sprfirst_*.csv in ", dir)
d <- do.call(rbind, lapply(fs, read.csv))
arms <- c("nni", "spr_ex", "spr_un", "none")

keys <- unique(d[order(d$nTip), "key"])
cat(sprintf("%-13s %5s %5s | %-28s | %-28s | %8s %8s\n",
            "key", "nTip", "nChr", "mean handoff score", "mean wall (s)",
            "dScr", "dWall"))
cat(strrep("-", 104), "\n")
go <- data.frame()
for (k in keys) {
  dk <- d[d$key == k, ]
  sc <- sapply(arms, function(a) mean(dk$score[dk$arm == a]))
  wl <- sapply(arms, function(a) mean(dk$wall_s[dk$arm == a]))
  dScore <- sc["spr_ex"] - sc["nni"]; dWall <- wl["spr_ex"] - wl["nni"]
  cat(sprintf("%-13s %5d %5d | %6.0f %6.0f %6.0f %6.0f | %6.1f %6.1f %6.1f %6.1f | %+8.1f %+7.1f%%\n",
              k, dk$nTip[1], dk$nChar[1],
              sc["nni"], sc["spr_ex"], sc["spr_un"], sc["none"],
              wl["nni"], wl["spr_ex"], wl["spr_un"], wl["none"],
              dScore, 100 * dWall / wl["nni"]))
  go <- rbind(go, data.frame(key = k, nTip = dk$nTip[1], dScore = dScore,
                             dWall = dWall, exactness = sc["spr_ex"] - sc["spr_un"]))
}
cat("\ncols: score = nni | spr_ex | spr_un | none ; wall likewise ; dScr = spr_ex-nni (neg = SPR better)\n")
cat("\n== VERDICT INPUTS ==\n")
large <- go[go$nTip >= 150, ]
cat(sprintf("LARGE keys (>=150t): spr_ex beats nni on %d/%d (dScore<0); mean dScore=%+.2f\n",
            sum(large$dScore < -1e-9), nrow(large), mean(large$dScore)))
cat(sprintf("exactness effect (spr_ex-spr_un) mean=%+.2f  -> if ~0, the fix reopened nothing\n",
            mean(go$exactness)))
cat(sprintf("ALL keys: spr_ex beats nni on %d/%d\n", sum(go$dScore < -1e-9), nrow(go)))
