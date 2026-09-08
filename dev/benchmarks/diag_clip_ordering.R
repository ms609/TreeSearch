# diag_clip_ordering.R
#
# Diagnostic script for the size-weighted TBR clip ordering experiment.
#
# Purpose: Characterise baseline (random) TBR clip ordering behaviour to test
# whether the small-clip-first hypothesis holds empirically.
#
# For each dataset and seed, builds a random Wagner starting tree, runs
# ts_tbr_diagnostics() to convergence, and accumulates per-pass records.
# Produces three summary tables:
#
#   1. Accepted clip size breakdown by bucket (tips / small / large).
#      Key question: are tip clips over-represented in accepted moves
#      relative to their uniform expectation?
#
#   2. Clips tried before acceptance (productive passes).
#      Key question: is n_clips_tried typically large enough that a
#      small-first ordering could meaningfully reduce it?
#
#   3. Evaluation budget split: productive vs null passes.
#      Key question: what fraction of TBR work is "wasted" in null passes?
#
# Usage: Rscript dev/benchmarks/diag_clip_ordering.R [lib_path]
#   lib_path defaults to ".agent-wc"

args <- commandArgs(trailingOnly = TRUE)
lib_path <- if (length(args) >= 1) args[1] else ".agent-wc"

library(TreeSearch, lib.loc = lib_path)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

DATASETS <- c("Vinther2008", "Agnarsson2004", "Zhu2013", "Dikow2009")
SEEDS <- c(1847L, 2956L, 3712L, 4519L, 5823L, 6401L, 7238L, 8145L, 9032L, 9871L)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

prepare <- function(name) {
  ds <- TreeSearch::inapplicable.phyData[[name]]
  at <- attributes(ds)
  list(
    name     = name,
    contrast = at$contrast,
    tip_data = matrix(unlist(ds, use.names = FALSE),
                      nrow = length(ds), byrow = TRUE),
    weight   = at$weight,
    levels   = at$levels,
    n_taxa   = length(ds)
  )
}

# Bucket label for a clip of subtree size s given n_tip.
# Tip: s == 1
# Small: 2 <= s <= floor(sqrt(n_tip))
# Large: s > floor(sqrt(n_tip))
clip_bucket <- function(s, n_tip) {
  sq <- floor(sqrt(n_tip))
  ifelse(s == 1, "tip",
         ifelse(s <= sq, "small", "large"))
}

# Expected fraction of clips in each bucket for a binary rooted tree with
# n_tip leaves.  Total clips = 2*(n_tip-1).
#   Tip clips (s==1)   : exactly n_tip
#   Non-tip clips       : n_tip - 2
#   Among non-tip, sizes 2..n_tip-1.  Approximate uniform distribution:
#     small (2..floor(sqrt)) : floor(sqrt)-1 sizes
#     large (floor(sqrt)+1..n_tip-1): n_tip-1-floor(sqrt) sizes
# This is approximate (not all sizes appear equally often), but adequate
# for comparison against observed acceptance fractions.
expected_bucket_fracs <- function(n_tip) {
  n_clips   <- 2L * (n_tip - 1L)
  sq        <- floor(sqrt(n_tip))
  n_tip_c   <- n_tip          # tip clips
  n_nontip  <- n_tip - 2L     # non-tip clips
  n_small_c <- sq - 1L        # sizes 2..sq (approximate, may be 0)
  n_large_c <- n_nontip - n_small_c
  list(
    tip   = n_tip_c   / n_clips,
    small = max(0, n_small_c) / n_clips,
    large = max(0, n_large_c) / n_clips
  )
}

# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

cat("Collecting TBR pass diagnostics (", length(SEEDS), "seeds per dataset)...\n\n",
    sep = "")

all_records <- list()

for (dname in DATASETS) {
  d     <- prepare(dname)
  n_tip <- d$n_taxa
  sq    <- floor(sqrt(n_tip))
  exp   <- expected_bucket_fracs(n_tip)

  cat(sprintf("Dataset: %-15s  n_tip=%d  sqrt_n=%d  total_clips=%d\n",
              dname, n_tip, sq, 2L*(n_tip-1L)))

  ds_records <- vector("list", length(SEEDS))

  for (i in seq_along(SEEDS)) {
    set.seed(SEEDS[i])

    # Random Wagner starting tree
    wag <- TreeSearch:::ts_random_wagner_tree(
      d$contrast, d$tip_data, d$weight, d$levels
    )

    # TBR to convergence with per-pass diagnostics (default clip_order = RANDOM)
    res <- TreeSearch:::ts_tbr_diagnostics(
      wag$edge, d$contrast, d$tip_data, d$weight, d$levels
    )

    passes <- res$passes
    passes$dataset     <- dname
    passes$seed        <- SEEDS[i]
    passes$n_tip       <- n_tip
    passes$n_clips     <- 2L * (n_tip - 1L)
    passes$final_score <- res$score
    passes$bucket      <- clip_bucket(passes$accepted_clip_size, n_tip)
    # bucket is only meaningful for productive passes; set NA for null passes
    passes$bucket[!passes$productive] <- NA_character_

    ds_records[[i]] <- passes
  }

  all_records[[dname]] <- do.call(rbind, ds_records)
  recs <- all_records[[dname]]
  prod <- recs[recs$productive, ]
  null <- recs[!recs$productive, ]

  cat(sprintf("  Passes: %d  productive=%d (%.0f%%)  null=%d (%.0f%%)\n",
    nrow(recs), nrow(prod), 100*nrow(prod)/nrow(recs),
    nrow(null), 100*nrow(null)/nrow(recs)))

  if (nrow(prod) > 0) {
    tip_obs <- mean(prod$accepted_clip_size == 1)
    tip_exp <- exp$tip
    enrich  <- tip_obs / tip_exp
    cat(sprintf("  Tip-clip acceptance: observed=%.0f%%  expected=%.0f%%  enrichment=%.2fx\n",
      100*tip_obs, 100*tip_exp, enrich))
    cat(sprintf("  Clips tried before accept: median=%d  mean=%.1f  (out of %d clips)\n",
      median(prod$n_clips_tried), mean(prod$n_clips_tried), 2L*(n_tip-1L)))
    cat(sprintf("  Final score range: %.0f – %.0f\n",
      min(recs$final_score), max(recs$final_score)))
  }
  cat("\n")
}

combined <- do.call(rbind, all_records)
prod_all <- combined[combined$productive, ]

# ---------------------------------------------------------------------------
# Table 1: Accepted clip size bucket breakdown
# ---------------------------------------------------------------------------

cat("=== Table 1: Accepted clip size breakdown (productive passes only) ===\n\n")

fmt_pct <- function(x) sprintf("%.1f%%", 100 * x)

bucket_tbl <- do.call(rbind, lapply(DATASETS, function(dname) {
  p     <- prod_all[prod_all$dataset == dname, ]
  n_tip <- p$n_tip[1]
  exp   <- expected_bucket_fracs(n_tip)
  tot   <- nrow(p)

  tip_obs   <- mean(p$accepted_clip_size == 1)
  small_obs <- mean(p$accepted_clip_size > 1 &
                    p$accepted_clip_size <= floor(sqrt(n_tip)))
  large_obs <- mean(p$accepted_clip_size > floor(sqrt(n_tip)))

  data.frame(
    dataset         = dname,
    n_tip           = n_tip,
    n_prod_passes   = tot,
    tip_obs         = fmt_pct(tip_obs),
    tip_exp         = fmt_pct(exp$tip),
    tip_enrichment  = round(tip_obs / exp$tip, 2),
    small_obs       = fmt_pct(small_obs),
    small_exp       = fmt_pct(exp$small),
    large_obs       = fmt_pct(large_obs),
    large_exp       = fmt_pct(exp$large)
  )
}))

print(bucket_tbl, row.names = FALSE)

# ---------------------------------------------------------------------------
# Table 2: Clips tried before acceptance
# ---------------------------------------------------------------------------

cat("\n=== Table 2: Clips tried in productive passes ===\n")
cat("(n_clips_tried includes the accepted clip itself; 1 = first clip accepted)\n\n")

tried_tbl <- do.call(rbind, lapply(DATASETS, function(dname) {
  p        <- prod_all[prod_all$dataset == dname, ]
  n_clips  <- p$n_clips[1]
  tried    <- p$n_clips_tried

  data.frame(
    dataset          = dname,
    n_clips          = n_clips,
    n_prod_passes    = nrow(p),
    pct_first_clip   = fmt_pct(mean(tried == 1)),
    pct_within_5     = fmt_pct(mean(tried <= 5)),
    pct_within_10pct = fmt_pct(mean(tried <= 0.1 * n_clips)),
    median_tried     = median(tried),
    mean_tried       = round(mean(tried), 1),
    median_position  = round(median(tried) / n_clips, 2)
  )
}))

print(tried_tbl, row.names = FALSE)

# ---------------------------------------------------------------------------
# Table 3: Evaluation budget — productive vs null passes
# ---------------------------------------------------------------------------

cat("\n=== Table 3: Evaluation budget by pass type ===\n\n")

eval_tbl <- do.call(rbind, lapply(DATASETS, function(dname) {
  d    <- combined[combined$dataset == dname, ]
  prod <- d[d$productive, ]
  null <- d[!d$productive, ]
  tot  <- sum(d$n_candidates_evaluated)

  data.frame(
    dataset           = dname,
    n_prod_passes     = nrow(prod),
    n_null_passes     = nrow(null),
    pct_evals_prod    = fmt_pct(sum(prod$n_candidates_evaluated) / tot),
    pct_evals_null    = fmt_pct(sum(null$n_candidates_evaluated) / tot),
    med_evals_prod    = if (nrow(prod) > 0) median(prod$n_candidates_evaluated) else NA_real_,
    med_evals_null    = if (nrow(null) > 0) median(null$n_candidates_evaluated) else NA_real_
  )
}))

print(eval_tbl, row.names = FALSE)

# ---------------------------------------------------------------------------
# Hypothesis assessment
# ---------------------------------------------------------------------------

cat("\n=== Hypothesis assessment ===\n")
cat("H: small clips (s=1) are over-represented in accepted moves,\n")
cat("   AND n_clips_tried is large enough that ordering would help.\n\n")

for (dname in DATASETS) {
  p        <- prod_all[prod_all$dataset == dname, ]
  n_clips  <- p$n_clips[1]
  n_tip    <- p$n_tip[1]
  enrich   <- (mean(p$accepted_clip_size == 1)) /
              (expected_bucket_fracs(n_tip)$tip)
  med_pos  <- median(p$n_clips_tried) / n_clips  # fraction of clips needed

  # Potential saving if tips-first: E[position of accepted tip clip in random
  # order] - E[position in tips-first order]. Very roughly:
  # random E[pos] ≈ n_clips/2; tips-first E[pos] ≈ n_tip/2.
  # saving_fraction ≈ (n_clips/2 - n_tip/2) / n_clips = (1 - n_tip/n_clips)/2 ≈ 0.25
  # But only beneficial if tip clips ARE more commonly accepted (enrich > 1).

  verdict <- if (enrich >= 2.0 && med_pos >= 0.25) {
    "STRONGLY SUPPORTS ordering (high enrichment + late acceptance)"
  } else if (enrich >= 1.5 && med_pos >= 0.15) {
    "SUPPORTS ordering (moderate enrichment + moderate position)"
  } else if (enrich >= 1.5) {
    "PARTIAL (enrichment, but acceptance mostly in first few clips)"
  } else if (enrich < 0.8) {
    "CONTRADICTS hypothesis (large clips accepted more often)"
  } else {
    "NEUTRAL (no consistent tip-clip enrichment)"
  }

  cat(sprintf("  %-15s: tip enrichment=%.2fx  median_pos=%.2f  -> %s\n",
    dname, enrich, med_pos, verdict))
}

cat("\nDone.\n")
