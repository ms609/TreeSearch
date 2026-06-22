# IW element isolation (recipe-agnostic).
#
# Runs ts_bench_tbr_phases (a pure TBR neighbourhood scan, NO moves applied) in
# EW vs IW on the IDENTICAL start tree + candidate set, so the per-element deltas
# are the cost IW *adds over EW on the same machinery*, independent of how much
# wall any recipe happens to spend in a phase.
#
# Elements isolated (all per-UNIT, so recipe-proportion-free):
#   per-candidate  = time_indirect_us*1000 / n_candidates  (ns/candidate)
#       EW: SIMD any_hit_reduce + branchless popcount
#       IW: SAME any_hit_reduce + scalar ctz-gather of iw_delta[pattern_index]
#       => the IW delta here = the gather-vs-popcount cost (lever #3)
#   per-clip       = time_clip_incr_us / n_clips            (us/clip)
#       IW adds: extract_char_steps (O(nodes*blocks)) + compute_iw base
#                + precompute_iw_delta (2 div/pattern)      (levers #1,#2)
#   full-rescore   = time_full_rescore_us                   (us, 1 call)
#       IW adds: extract_char_steps + compute_iw (div/pattern) (lever #4)
#
# has_na is forced FALSE (recode "-"->"?") to isolate the PURE-IW directional
# path (EW kernel + IW gather) -- the cleanest EW-vs-IW contrast, matching the
# EW program's round-7 standard-Fitch methodology. NA+IW is a separate element.
#
# Usage: Rscript dev/benchmarks/bench_iw_element_isolation.R [lib] [concavity] [reps]
suppressMessages(suppressWarnings({
  args <- commandArgs(trailingOnly = TRUE)
  LIB  <- if (length(args) >= 1) args[[1]] else
    "C:/Users/pjjg18/GitHub/worktrees/TreeSearch/confident-gates-0f627e/.agent-tbr"
  CONC <- if (length(args) >= 2) as.numeric(args[[2]]) else 10  # IW default k=10
  REPS <- if (length(args) >= 3) as.integer(args[[3]]) else 7L
  library(TreeSearch, lib.loc = LIB)
  library(TreeTools)
}))
data("inapplicable.phyData", package = "TreeSearch")

# Recode "-" -> "?" (drop inapplicable level; map the "-" token to all-states)
# so make_dataset sees inapp_state = -1 => pure Fitch / pure-IW path.
recode_na_off <- function(phy) {
  at <- attributes(phy); ct <- at$contrast; lv <- at$levels; al <- at$allLevels
  inapp_col <- which(lv == "-"); dash_row <- which(al == "-")
  if (length(dash_row)) ct[dash_row, ] <- 1
  if (length(inapp_col)) { ct <- ct[, -inapp_col, drop = FALSE]; lv <- lv[-inapp_col] }
  attr(phy, "contrast") <- ct; attr(phy, "levels") <- lv
  phy
}

prep <- function(phy) {
  at <- attributes(phy)
  list(contrast = at$contrast,
       tip = matrix(unlist(phy, use.names = FALSE),
                    nrow = length(phy), byrow = TRUE),
       weight = at$weight, levels = at$levels, n = length(phy),
       ns = ncol(at$contrast))
}

bench_one <- function(p, edge, conc) {
  # median over REPS independent scans for timing stability
  fr <- ci <- ind <- numeric(REPS); nc <- ncand <- NA_integer_
  for (r in seq_len(REPS)) {
    set.seed(1000 + r)
    res <- TreeSearch:::ts_bench_tbr_phases(
      edge, p$contrast, p$tip, p$weight, p$levels, integer(0), conc)
    fr[r]  <- res$time_full_rescore_us
    ci[r]  <- res$time_clip_incr_us
    ind[r] <- res$time_indirect_us
    nc     <- res$n_clips
    ncand  <- res$n_candidates
  }
  list(full = median(fr), clip = median(ci), indirect = median(ind),
       n_clips = nc, n_cand = ncand)
}

DATASETS <- c("Vinther2008", "Wortley2006", "Zanol2014", "Giles2015", "Dikow2009")

cat(sprintf("IW element isolation | lib=%s | concavity=%g | reps=%d\n",
            basename(LIB), CONC, REPS))
cat(sprintf("%-13s %3s %2s | %5s %8s | %-22s | %-22s | %-18s\n",
            "dataset", "n", "ns", "clips", "cands",
            "per-cand ns (EW/IW/Δ)", "per-clip us (EW/IW/Δ)",
            "full us (EW/IW)"))

for (nm in DATASETS) {
  phy <- recode_na_off(inapplicable.phyData[[nm]])
  p <- prep(phy)
  set.seed(7294)
  edge <- RandomTree(names(inapplicable.phyData[[nm]]), root = TRUE)$edge
  ew <- bench_one(p, edge, -1.0)
  iw <- bench_one(p, edge, CONC)

  # per-unit
  pc_ew <- ew$indirect * 1000 / ew$n_cand
  pc_iw <- iw$indirect * 1000 / iw$n_cand
  cl_ew <- ew$clip / ew$n_clips
  cl_iw <- iw$clip / iw$n_clips

  cat(sprintf("%-13s %3d %2d | %5d %8d | %6.1f/%6.1f/%+5.1f (%.2fx) | %6.1f/%6.1f/%+5.1f (%.2fx) | %6.0f/%6.0f (%.2fx)\n",
              nm, p$n, p$ns, ew$n_clips, ew$n_cand,
              pc_ew, pc_iw, pc_iw - pc_ew, pc_iw / pc_ew,
              cl_ew, cl_iw, cl_iw - cl_ew, cl_iw / cl_ew,
              ew$full, iw$full, iw$full / ew$full))
}
