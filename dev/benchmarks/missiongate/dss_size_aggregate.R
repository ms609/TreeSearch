# Aggregate the DSS size-axis probe into the build-or-negative verdict.
#
# METRIC NOTE (load-bearing; see advisor reconciliation 2026-07-16): steps/Mcand is a
# BIASED axis for the small-vs-large comparison. candidates_evaluated is a raw count of
# scored topologies; it does not know a maxS=12 candidate is scored on a ~13-tip REDUCED
# dataset. In this engine small-sector candidates are in fact ~2x MORE expensive per
# WALL than ratchet's (per-pick reduced-dataset build + 3x RAS overhead spread over
# fewer candidates), so steps/Mcand and steps/wall RANK THE ARMS DIFFERENTLY. rev5 could
# lean on steps/Mcand because (a) its selectem used maxS>=0.65n+collapse (candidate scale
# near ratchet's) and (b) the bias ran AGAINST its negative (conservative). Using the same
# instrument to claim sectorial is BETTER feeds the bias. So this aggregate leads with the
# WALL axis (steps/sec) and REACH-TO-TARGET (best escape), and reports steps/Mcand only as
# a flagged, non-adjudicating column. NEVER diff Mcand across configs (rev2: goes negative).
#
# DECISION on the WALL axis:
#   - small sectors do NOT beat ratchet on steps/sec AND do NOT reach the target =>
#     NEGATIVE (rev5's verdict holds + extends to the small-maxS regime it never measured);
#     diameter-limiting is strictly MORE local => shallower reach => cannot help. No build.
#   - small sectors beat ratchet on steps/sec on the HARD targets AND reach them => LEAD.

indir <- Sys.getenv("OUT_DIR", "dev/benchmarks/missiongate/dss_out")
files <- list.files(indir, pattern = "^cell_.*\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
D <- do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
D$cfg <- ifelse(D$arm == "ratchet_ref", "ratchet_ref",
         ifelse(D$arm == "selectem_large", sprintf("selectem_large r=%d", D$rounds),
                sprintf("smallmax maxS=%d r=%d", D$maxS, D$rounds)))

agg <- do.call(rbind, lapply(split(D, list(D$dataset, D$cfg), drop = TRUE), function(g) {
  em <- median(g$escape); wm <- median(g$wall_s); mc <- median(g$Mcand)
  data.frame(dataset = g$dataset[1], arm = g$arm[1], cfg = g$cfg[1], maxS = g$maxS[1],
    n = nrow(g), start = g$start[1], tnt_sect = g$tnt_sect[1],
    esc_med = em, esc_best = min(g$escape), wall_med = round(wm, 3),
    Mcand_med = round(mc, 3),
    st_per_sec  = if (wm > 0) round(em / wm, 1) else NA,
    st_per_Mc   = if (mc > 0) round(em / mc, 2) else NA,   # BIASED — see header
    stringsAsFactors = FALSE)
}))

verdict <- character(0)
for (ds in unique(agg$dataset)) {
  a <- agg[agg$dataset == ds, ]
  rr <- a[a$arm == "ratchet_ref", ]
  r_sec <- if (nrow(rr)) rr$st_per_sec else NA; r_esc <- if (nrow(rr)) rr$esc_med else NA
  r_reach <- if (nrow(rr)) rr$esc_best else NA; r_wall <- if (nrow(rr)) rr$wall_med else NA
  tgt <- a$tnt_sect[1] - a$start[1]
  a <- a[order(match(a$arm, c("ratchet_ref","smallmax","selectem_large")), a$maxS, a$cfg), ]
  cat(sprintf("\n== %s (start=%.0f, tnt target=%+d; ratchet esc=%+d reach=%+d @ %.2fs = %s st/s) — WALL AXIS ==\n",
      ds, a$start[1], as.integer(tgt), as.integer(r_esc), as.integer(r_reach), r_wall,
      ifelse(is.na(r_sec),"NA",sprintf("%.1f",r_sec))))
  cat(sprintf("   %-24s %4s %8s %8s %8s %11s %11s\n",
      "config","n","esc_med","reach","wall_s","st/s(WALL)","st/Mc[bias]"))
  for (i in seq_len(nrow(a))) {
    beats_wall <- a$arm[i] != "ratchet_ref" && !is.na(r_sec) &&
                  a$wall_med[i] <= r_wall * 1.05 && a$esc_med[i] < r_esc
    cat(sprintf("   %-24s %4d %+8d %+8d %8.2f %11.1f %11.2f%s\n",
        a$cfg[i], a$n[i], as.integer(a$esc_med[i]), as.integer(a$esc_best[i]),
        a$wall_med[i], a$st_per_sec[i], a$st_per_Mc[i],
        if (beats_wall) "  <<BEATS-RATCHET@WALL" else ""))
  }
  sm <- a[a$arm == "smallmax", ]
  sm_reach <- if (nrow(sm)) min(sm$esc_best) else NA       # deepest small-sector escape
  sm_wall_win <- any(sm$wall_med <= r_wall * 1.05 & sm$esc_med < r_esc, na.rm = TRUE)
  reaches <- !is.na(sm_reach) && !is.na(tgt) && sm_reach <= tgt
  verdict[ds] <- if (sm_wall_win && reaches) "WALL-WIN+REACHES" else
                 if (sm_wall_win) "WALL-WIN-BUT-SHALLOW" else "RATCHET-DOMINATES"
  cat(sprintf("   => %s: deepest smallmax reach=%+d (target %+d, %s); beats ratchet@wall=%s\n",
      ds, as.integer(sm_reach), as.integer(tgt),
      if (reaches) "REACHES" else "PLATEAUS", sm_wall_win))
}

cat("\n===== BUILD-OR-NEGATIVE VERDICT (WALL axis) =====\n")
for (ds in names(verdict)) cat(sprintf("  %-12s %s\n", ds, verdict[ds]))
nwin <- sum(verdict == "WALL-WIN+REACHES")
cat(sprintf("\nDatasets where small sectors beat ratchet on wall AND reach the target: %d / %d\n",
            nwin, length(verdict)))
cat(if (nwin == 0)
  paste0("=> NEGATIVE: on the WALL axis small sectors do not both beat ratchet and reach\n",
         "   the target. The steps/Mcand 'win' is a candidate-count artifact. Diameter-\n",
         "   limiting is strictly more local (shallower) => cannot help. dss MISSION-DEAD.\n")
  else "=> LEAD on the wall axis; the diameter build may be warranted — reconsult before building.\n")

write.csv(agg, file.path(indir, "dss_aggregate.csv"), row.names = FALSE)
cat(sprintf("\nWrote %s\n", file.path(indir, "dss_aggregate.csv")))
