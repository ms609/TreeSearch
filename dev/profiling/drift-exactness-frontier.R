#!/usr/bin/env Rscript
# Authoritative matched-wall frontier analysis for the Hamilton gate (30 seeds).
# Score achievable by wall <= W: per seed+scorer, min score among runs with
# wall_s <= W (score is monotone non-increasing in nCycles => in wall); mean and
# best over the 30 seeds. Reports union vs exact and the crossover per dataset.
dir <- "dev/profiling"
files <- list.files(dir, pattern = "^drift-exactness-gate-hamilton-.*\\.csv$",
                    full.names = TRUE)
df <- do.call(rbind, lapply(files, read.csv))
cat(sprintf("rows=%d datasets=%s seeds=%d nCycles=%s\n\n",
            nrow(df), paste(unique(df$dataset), collapse=","),
            length(unique(df$seed)), paste(sort(unique(df$nCycles)), collapse=",")))

# Raw mean score + mean wall by dataset x nCycles x scorer (matched-cycles view)
cat("== mean score / mean wall(s) by dataset x nCycles x scorer ==\n")
ag <- aggregate(cbind(score, wall_s) ~ dataset + nCycles + scorer, df, mean)
ag <- ag[order(ag$dataset, ag$nCycles, ag$scorer), ]
ag$score <- round(ag$score, 2); ag$wall_s <- round(ag$wall_s, 3)
print(ag, row.names = FALSE)

# Matched-WALL frontier
scoreAtWall <- function(d, W) {
  seeds <- unique(d$seed)
  v <- sapply(seeds, function(s){
    x <- d[d$seed==s & d$wall_s<=W, ]; if (nrow(x)==0) NA else min(x$score)
  })
  c(mean=mean(v, na.rm=TRUE), best=min(v, na.rm=TRUE), cov=mean(!is.na(v)))
}
Ws <- c(0.1,0.2,0.4,0.8,1.6,3.2,6.4,12.8)
for (nm in unique(df$dataset)) {
  cat(sprintf("\n== %s : matched-WALL frontier (mean score over 30 seeds; * = coverage<1) ==\n", nm))
  cat(sprintf("  %6s %10s %10s   %s\n","wall","union","exact","winner(margin)"))
  for (W in Ws) {
    u <- scoreAtWall(df[df$dataset==nm & df$scorer=="union",], W)
    e <- scoreAtWall(df[df$dataset==nm & df$scorer=="exact",], W)
    mk <- function(z) if (z["cov"]<1) "*" else " "
    win <- if (any(is.na(c(u["mean"],e["mean"])))) "-" else
           if (abs(u["mean"]-e["mean"])<1e-9) "tie" else
           if (e["mean"]<u["mean"]) sprintf("EXACT (%.2f)", u["mean"]-e["mean"]) else
                                    sprintf("union (%.2f)", e["mean"]-u["mean"])
    cat(sprintf("  %6.2f %9.2f%s %9.2f%s   %s\n", W, u["mean"], mk(u), e["mean"], mk(e), win))
  }
  # saturated reach (best score at max budget, over seeds)
  du <- df[df$dataset==nm & df$scorer=="union",]; de <- df[df$dataset==nm & df$scorer=="exact",]
  cat(sprintf("  saturated best-of-30 (any budget): union=%.0f exact=%.0f | mean@maxNc: union=%.2f exact=%.2f\n",
      min(du$score), min(de$score),
      mean(du$score[du$nCycles==max(du$nCycles)]), mean(de$score[de$nCycles==max(de$nCycles)])))
}
