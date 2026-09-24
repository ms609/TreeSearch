# tbr_throughput_race.R — DISCRIMINATOR #1 (advisor 2026-06-19): is the
# gold-standard 2–3× per-candidate gap (T-P5h2) actually IN TBR, or does it live
# in the ratchet/sectorial trajectory?
#
# The framing.R whole-search throughput is CONFOUNDED (TS runs to its 50-rep cap
# while TNT xmult stops at hits 10; phase mixes differ).  This driver removes
# that confound: SAME operator (strict-descent TBR), SAME starting tree, run to
# convergence on BOTH engines — only the search trajectory can differ.
#
#   rate    = rearrangements / engine_wall   (Mcand-or-rearr / s)
#   thrpt   = tnt_rate / ts_rate
#
# TNT wall includes startup/proc/tread (INFLATES it) ⇒ tnt rate is a LOWER bound,
# so a measured thrpt is conservative.  Aimed at the LOSERS (Zanol, Zhu) + Giles
# as a winner control.
#
#   thrpt ≈ framing's 2–2.5×  ⇒ gap is IN TBR  → bookkeeping candidate (L3b/cluster)
#   thrpt ≈ 1                  ⇒ gap is TRAJECTORY (ratchet/sectorial) → L3b wrong lever
#
# Env: TS_LIB (default .agent-l1 = HEAD 00d73d6a), TNT_EXE, RACE_DATASETS,
#      RACE_SEEDS, RACE_WAGSEED.
source("dev/benchmarks/tbr_shared_start_lib.R")

dsN   <- strsplit(trimws(Sys.getenv("RACE_DATASETS",
           "Zanol2014 Zhu2013 Giles2015")), "\\s+")[[1]]
seeds <- as.integer(strsplit(trimws(Sys.getenv("RACE_SEEDS", "1 2 3")), "\\s+")[[1]])
wagSeed <- as.integer(Sys.getenv("RACE_WAGSEED", "11"))

cat(sprintf("TBR THROUGHPUT RACE | lib=%s | datasets {%s} | seeds {%s}\n",
            Sys.getenv("TS_LIB", ".agent-l1"), paste(dsN, collapse = ","),
            paste(seeds, collapse = ",")))

allRows <- list()
for (nm in dsN) {
  d <- prepareDataset(nm)
  # Poor Wagner start ⇒ long climb ⇒ swap-dominated wall (clean rate).
  set.seed(wagSeed)
  wag <- TreeSearch:::ts_random_wagner_tree(d$contrast, d$tip_data, d$weight,
                                            d$levels)
  wagTree <- Preorder(RenumberTips(
    structure(list(edge = wag$edge, Nnode = d$nTip - 1L,
                   tip.label = names(d$phy)), class = "phylo"), names(d$phy)))
  startLen <- TreeLength(wagTree, d$phy)
  cat(sprintf("\n==== %s (%dt) | Wagner(seed %d) start_len = %.0f ====\n",
              nm, d$nTip, wagSeed, startLen))
  rows <- list()
  for (s in seeds) {
    rows[[length(rows) + 1]] <- cbind(dataset = nm, tips = d$nTip,
        TntTbr(d, wagTree, seed = s, mulpars = FALSE, hold = 1))
    rows[[length(rows) + 1]] <- cbind(dataset = nm, tips = d$nTip,
        TsTbr(d, wagTree, seed = s, acceptEqual = FALSE)$row)
  }
  dd <- do.call(rbind, rows)
  dd$rate <- dd$rearrangements / dd$wall / 1e6          # Mcand-or-rearr / s
  allRows[[nm]] <- dd
  print(dd[, c("engine", "seed", "final_len", "rearrangements", "wall", "rate")],
        row.names = FALSE)
}
res <- do.call(rbind, allRows)

cat("\n=== PER-DATASET MEDIAN (strict-descent TBR from shared Wagner start) ===\n")
agg <- do.call(rbind, lapply(split(res, res$dataset), function(z) {
  tnt <- z[z$engine == "TNT", ]; ts <- z[z$engine == "TS", ]
  data.frame(dataset = z$dataset[1], tips = z$tips[1],
    ts_final = median(ts$final_len), tnt_final = median(tnt$final_len),
    ts_cand = round(median(ts$rearrangements) / 1e6, 1),
    tnt_rearr = round(median(tnt$rearrangements) / 1e6, 1),
    ts_rate = round(median(ts$rate), 2), tnt_rate = round(median(tnt$rate), 2),
    throughput = round(median(tnt$rate) / median(ts$rate), 2),
    ts_wall = round(median(ts$wall), 2), tnt_wall = round(median(tnt$wall), 2),
    stringsAsFactors = FALSE)
}))
agg <- agg[order(agg$tips), ]
print(agg, row.names = FALSE)

cat("\n--- READ ---\n")
cat("score parity (ts_final ~ tnt_final) is the validity gate: large gap = apples-to-oranges.\n")
cat("throughput ~ framing 2-2.5x => per-candidate gap is IN TBR (bookkeeping lever live).\n")
cat("throughput ~ 1            => gap is ratchet/sectorial TRAJECTORY, not TBR per-candidate.\n")

outCsv <- Sys.getenv("RACE_OUT", "dev/profiling/tbr_throughput_race.csv")
write.csv(res, outCsv, row.names = FALSE)
cat(sprintf("\nrows -> %s\n", outCsv))
