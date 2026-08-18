#!/usr/bin/env Rscript
# Fourth surgical patch: patch_findings.R wrote the signature correction from the
# 7-matrix LOCAL pilot, before the 25-matrix battery and the budget-96 sweep landed.
# The corpus result is stronger and supersedes it: the recorded signature's SIGN
# FLIPS with the replicate budget, which is a cleaner refutation than "cheaper per
# replicate", and the durable invariant is TBR work per replicate halving.
f <- Sys.getenv("FINDINGS", "C:/Users/pjjg18/GitHub/TreeSearch/dev/red-team/findings.md")
raw <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

old <- paste0(
  "Cumulative per-phase totals do rise (ratchet 1.20x, xss 1.50x, final_tbr 1.24x, ",
  "arm 2 vs pre-fix), but those totals are summed over however many replicates ran, ",
  "and a blocked arm runs MORE of them. Normalised PER REPLICATE, arm 2 is CHEAPER ",
  "in nearly every phase \u2014 tbr 0.61x, wagner 0.70x, ratchet 0.83x, fuse 0.45x ",
  "\u2014 which is what rejected moves must do: less work per replicate, not more. ",
  "Measured from the engine's own `timings` attribute (microsecond resolution, so ",
  "not a clock artefact), 13 cells x 3 arms. The replicate counts are the whole ",
  "effect: median `reps_done` 17 pre-fix / **24 = budget exhausted** / 16 merged.")

new <- paste0(
  "**The test's SIGN FLIPS WITH THE REPLICATE BUDGET, on the same code and the same ",
  "defect**, which is what disqualifies it. Cumulative per-phase complement-only vs ",
  "pre-fix: at `maxReplicates = 24` the blocked arm is *faster* in 6 of 7 live ",
  "phases (ratchet 0.79x, tbr 0.44x, wagner 0.39x, fuse 0.24x; only xss 1.26x up), ",
  "while at `maxReplicates = 96` it is *slower* in 6 of 7 (ratchet 1.60x, xss ",
  "1.95x, rss 1.67x, final_tbr 1.65x, fuse 2.10x). A diagnostic that reverses its ",
  "verdict when a budget setting changes cannot be evidence of a defect. The ",
  "arithmetic reason: phase totals are summed over however many replicates ran, and ",
  "the replicate COUNT is precisely what this defect changes -- budget 24, ",
  "`reps_done` 17 pre-fix / **24 = exhausted** / 16 merged; budget 96 on the ",
  "affected geometry, **15 / 65 / 17**. **The invariant that IS diagnostic, holding ",
  "at both budgets, is TBR work per replicate roughly HALVING** in the blocked arm ",
  "-- 0.44x cumulative at budget 24, 0.498x per-replicate at budget 96 -- which is ",
  "what rejected moves must do: less work per replicate, not more. Per-replicate ",
  "total wall drops too (0.032 s vs 0.058 pre-fix / 0.044 merged). Measured from the ",
  "engine's own `timings` attribute, microsecond resolution, so not a clock ",
  "artefact. **Never compare cumulative phase totals across runs whose replicate ",
  "counts differ.**")

n <- gregexpr(old, raw, fixed = TRUE)[[1]]
if (identical(n[1], -1L)) stop("PATTERN NOT FOUND -- did patch_findings.R run?")
if (length(n) != 1L) stop("matched ", length(n), " times, expected 1")
raw <- sub(old, new, raw, fixed = TRUE)

con <- file(f, open = "wb")            # binary: preserve LF
writeLines(raw, con, sep = "\n", useBytes = TRUE)
close(con)
cat("patch 4 applied\n")
