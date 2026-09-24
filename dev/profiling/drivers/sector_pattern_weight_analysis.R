d <- readRDS("C:/Users/pjjg18/GitHub/TreeSearch/dev/profiling/sector_pattern_killtest.rds")

cat("==== Work-weighted block_reduction (TBR cost proxy) ====\n")
for (w in c("count","lin","quad")) {
  d$cost <- switch(w, count = 1, lin = d$sector_tips, quad = d$sector_tips^2)
  agg <- tapply(seq_len(nrow(d)), d$dataset, function(i)
    weighted.mean(d$block_reduction[i], d$cost[i]))
  cat(sprintf("  [%-5s] ", w))
  for (nm in names(agg)) cat(sprintf("%s=%.2f  ", nm, agg[nm]))
  cat("\n")
}

cat("\n==== Also: fractional block reduction (compact/orig) work-weighted ====\n")
d$frac_remaining <- d$compact_blocks / d$orig_blocks
for (w in c("lin","quad")) {
  d$cost <- if (w=="lin") d$sector_tips else d$sector_tips^2
  agg <- tapply(seq_len(nrow(d)), d$dataset, function(i)
    weighted.mean(d$frac_remaining[i], d$cost[i]))
  cat(sprintf("  [%-5s] frac of blocks REMAINING (1.0 = no win): ", w))
  for (nm in names(agg)) cat(sprintf("%s=%.2f  ", nm, agg[nm]))
  cat("\n")
}

cat("\n==== Top-cost decile of sectors (by sector_tips) per dataset ====\n")
for (nm in unique(d$dataset)) {
  dd <- d[d$dataset == nm, ]
  thr <- quantile(dd$sector_tips, 0.9)
  top <- dd[dd$sector_tips >= thr, ]
  cat(sprintf("  %-7s top-decile (tips>=%.0f, n=%d): med tips=%.0f, med n_inf=%.0f, orig=%d, med compact=%.1f, med reduction=%.1f\n",
              nm, thr, nrow(top), median(top$sector_tips), median(top$n_inf),
              top$orig_blocks[1], median(top$compact_blocks), median(top$block_reduction)))
}

cat("\n==== Largest single sector per dataset ====\n")
for (nm in unique(d$dataset)) {
  dd <- d[d$dataset == nm, ]
  big <- dd[which.max(dd$sector_tips), ]
  cat(sprintf("  %-7s biggest: tips=%d, n_inf=%d/%d, blocks %d->%d (reduction %d)\n",
              nm, big$sector_tips, big$n_inf, big$npat, big$orig_blocks, big$compact_blocks, big$block_reduction))
}
