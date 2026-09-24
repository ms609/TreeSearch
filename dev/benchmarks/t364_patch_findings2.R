#!/usr/bin/env Rscript
# Second surgical patch: add the measured MAGNITUDES to the T-384 row, now that the
# 25-matrix battery and the production-budget sweep have landed.  Same
# match-exactly-once discipline as patch_findings.R.
f <- Sys.getenv("FINDINGS", "C:/Users/pjjg18/GitHub/TreeSearch/dev/red-team/findings.md")
raw <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

old <- "**Corollary: this row's \"roughly 1-in-11\" understates T-364 badly** — that is the common-geometry figure; on the minority geometry `AdditionTree()` violated up to 81% of the time."

new <- paste0(
  "**Corollary: this row's \"roughly 1-in-11\" understates T-364 badly** — that is the common-geometry figure; on the minority geometry `AdditionTree()` violated up to **94%** of the time (21 matrices x 200 identical addition orders per arm; the arm-1 violating and arm-2 complement-rooted distributions match in median AND range, `in0` 81.2% [71.0-94.0], `no0` 6.0% [0.0-11.5], arm 3 zero throughout). ",
  "**MAGNITUDE, measured on the 25-matrix `MBANK_FIXED_SAMPLE` (training only), 47 matrix x shape pairs x 3 seeds, three arms interleaved on one node per cell, zero censored runs.** Complement-without-reroot vs pre-fix wall: at `maxReplicates = 24`, `in0` 1.010 / `no0` 0.985 / pooled 0.995 — no cost; at the production `maxReplicates = 96`, `in0` **2.531** (11/13 slower, 11 >10% slower, sign p = 0.022), `no0` 1.053, pooled **1.211** (19/8, 16 >10%, p = 0.052). ",
  "**The ratio GROWS with the budget, which is the discriminating prediction:** budget exhaustion predicts growth, uniform per-unit slowness predicts invariance. Replicates run on `in0` at budget 96: pre-fix **15**, merged **17**, complement-only **65**; budget-24 exhaustion rate on `in0` 50% / 47% / **95.5%**. `merged/complement-only` on `in0` at budget 96 = **0.382** (1/12, p = 0.0034) — the re-root removes the cost entirely. ",
  "**And it reconciles the old measurement rather than contradicting it:** at budget 24 the complement-only arm's wall is flat but its SCORE is worse on 6 of 22 `in0` matrices, 0 better (p = 0.031), with the merged arm beating it on the same 6 (p = 0.031); given more budget it recovers the score and pays 2.53x to do so. The 1.43x was measured as \"wall to match the old score\" (score held fixed); this battery holds budget fixed. Same phenomenon, two ways of holding things constant, and the historical 1.43x sits between the pooled 1.21x and the stratum 2.53x — an exact match is not expected, the T-214 battery being 3 matrices of 10-15 tips at an unrecorded budget. ",
  "**THE MERGED FIX ITSELF COSTS NOTHING:** merged vs pre-fix wall median **1.000** (22 slower / 22 faster / 3 tie, 8 of 47 >10% slower but balanced by 9 >10% faster, p = 1) at budget 24, and **0.967** (10/17, 4 of 27 >10% slower, p = 0.25) at budget 96; score +0 with 2 worse / 2 better / 43 tie (p = 1). Compliance of returned trees: merged 100% of 2634, pre-fix 99.96% — **one violating tree in 2751**, the rejection sampler exhausting its 100 attempts, so the search path was very nearly but not perfectly shielded. Read the 1.000 as two effects cancelling (pre-fix pays reshuffles ~1/(1-p) builds per start, complement-only pays blocking, merged pays neither), not as \"nothing happens\". Harness and full tables: `dev/benchmarks/t364_*`, findings write-up `dev/benchmarks/t364_FINDINGS.md`.")

n <- gregexpr(old, raw, fixed = TRUE)[[1]]
if (identical(n[1], -1L)) stop("PATTERN NOT FOUND -- did patch_findings.R run?")
if (length(n) != 1L) stop("matched ", length(n), " times, expected 1")
raw <- sub(old, new, raw, fixed = TRUE)

con <- file(f, open = "wb")            # binary: preserve LF, do not re-introduce CRLF
writeLines(raw, con, sep = "\n", useBytes = TRUE)
close(con)
cat("patch 2 applied\n")
