# Aggregate the mission-gate cells into the wall-vs-escape verdict.
#
# The MISSION GATE: sectorial escape is a wall-clock win ONLY if some config
# reaches a better score than the ratchet it would displace at <= ratchet's wall.
# For each dataset we pool all T0-seeds x sectorial-RNG-seeds (n per config), and
# report median/best escape + median wall per config, then test whether ANY
# sectorial config Pareto-beats ratchet_ref (cheaper-or-equal wall AND better
# median escape). A clean "no" retires the sectorial-escape thread to the budget
# axis; a "yes" is a real cheap-escape finding worth a src walk-up follow-up.

indir <- Sys.getenv("OUT_DIR", "dev/benchmarks/missiongate/out")
files <- list.files(indir, pattern = "^cell_.*\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
D <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))

cfg <- function(d) {
  ifelse(d$arm == "ratchet_ref", "ratchet_ref",
    ifelse(d$arm == "baseline", sprintf("baseline r=%d", d$rounds),
      sprintf("selectem pk=%d ras=%d r=%d", d$rssPicks, d$rasStarts, d$rounds)))
}
D$cfg <- cfg(D)

agg <- do.call(rbind, lapply(split(D, list(D$dataset, D$cfg), drop = TRUE), function(g) {
  data.frame(dataset = g$dataset[1], arm = g$arm[1], cfg = g$cfg[1],
    n = nrow(g), start = g$start[1], tnt_sect = g$tnt_sect[1],
    esc_med = median(g$escape), esc_best = min(g$escape),
    wall_med = round(median(g$wall_s), 3), Mcand_med = round(median(g$Mcand), 2),
    stringsAsFactors = FALSE)
}))

cat("\n===== Per-dataset wall-vs-escape (median over all seeds; escape<0 = better than T0) =====\n")
for (ds in unique(agg$dataset)) {
  a <- agg[agg$dataset == ds, ]
  a <- a[order(a$wall_med), ]
  rr <- a[a$arm == "ratchet_ref", ]
  cat(sprintf("\n-- %s (start=%.0f, tnt_sect target=%+d) --\n",
              ds, a$start[1], as.integer(a$tnt_sect[1] - a$start[1])))
  cat(sprintf("   %-26s %4s %8s %8s %8s %8s\n", "config", "n", "esc_med", "esc_best", "wall_s", "Mcand"))
  for (i in seq_len(nrow(a))) {
    flag <- ""
    if (nrow(rr) == 1 && a$arm[i] != "ratchet_ref" &&
        a$wall_med[i] <= rr$wall_med * 1.05 && a$esc_med[i] < rr$esc_med) flag <- "  <<BEATS-RATCHET"
    cat(sprintf("   %-26s %4d %+8d %+8d %8.2f %8.1f%s\n",
                a$cfg[i], a$n[i], as.integer(a$esc_med[i]), as.integer(a$esc_best[i]),
                a$wall_med[i], a$Mcand_med[i], flag))
  }
}

# Mission-gate verdict: any sectorial config with wall <= ratchet AND better median escape?
cat("\n===== MISSION-GATE VERDICT =====\n")
win <- 0L
for (ds in unique(agg$dataset)) {
  a <- agg[agg$dataset == ds, ]
  rr <- a[a$arm == "ratchet_ref", ]
  if (nrow(rr) != 1) { cat(sprintf("%-12s: no ratchet_ref\n", ds)); next }
  cand <- a[a$arm != "ratchet_ref" & a$wall_med <= rr$wall_med * 1.05 & a$esc_med < rr$esc_med, ]
  if (nrow(cand)) {
    win <- win + 1L
    b <- cand[which.min(cand$esc_med), ]
    cat(sprintf("%-12s: CHEAP-ESCAPE — %s esc_med=%+d @ %.2fs beats ratchet esc_med=%+d @ %.2fs\n",
                ds, b$cfg, as.integer(b$esc_med), b$wall_med, as.integer(rr$esc_med), rr$wall_med))
  } else {
    # best sectorial regardless of wall, for context
    sct <- a[a$arm != "ratchet_ref", ]; bb <- sct[which.min(sct$esc_med), ]
    cat(sprintf("%-12s: ratchet-dominated — best sectorial %s esc_med=%+d @ %.2fs vs ratchet %+d @ %.2fs\n",
                ds, bb$cfg, as.integer(bb$esc_med), bb$wall_med,
                as.integer(rr$esc_med), rr$wall_med))
  }
}
cat(sprintf("\nDatasets with a cheap sectorial escape beating ratchet: %d / %d\n",
            win, length(unique(agg$dataset))))
cat(if (win == 0)
  "=> MISSION-DEAD: no sectorial config escapes cheaper than ratchet. Thread retires to budget axis.\n"
  else "=> LIVE LEAD: cheap sectorial escape exists; warrants the src walk-up follow-up.\n")

write.csv(agg, file.path(indir, "aggregate.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s\n", file.path(indir, "aggregate.csv")))
