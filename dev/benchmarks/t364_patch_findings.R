#!/usr/bin/env Rscript
# Surgical, verified literal patch of two sentences in the T-384 row of
# dev/red-team/findings.md.  Literal (fixed = TRUE) and asserted to match exactly
# once each, so it cannot silently no-op or hit the wrong row -- the file is shared
# with a concurrent session actively editing other rows.
f <- Sys.getenv("FINDINGS", "C:/Users/pjjg18/GitHub/TreeSearch/dev/red-team/findings.md")
x <- readLines(f, warn = FALSE, encoding = "UTF-8")
raw <- paste(x, collapse = "\n")

reps <- list(
  list(
    old = "Diagnostic signature worth reusing: every phase (TBR, XSS, ratchet) slower on a score-IDENTICAL trajectory — 226→220 and 233→223 at both builds — which means moves being rejected, not a different search being run.",
    new = "Diagnostic signature — CORRECTED 2026-07-29, do not reuse the original form: the claim as filed was \"every phase (TBR, XSS, ratchet) slower on a score-IDENTICAL trajectory ⇒ moves being rejected\". Cumulative per-phase totals do rise (ratchet 1.20x, xss 1.50x, final_tbr 1.24x, arm 2 vs pre-fix), but those totals are summed over however many replicates ran, and a blocked arm runs MORE of them. Normalized PER REPLICATE, arm 2 is CHEAPER in nearly every phase — tbr 0.61x, wagner 0.70x, ratchet 0.83x, fuse 0.45x — which is what rejected moves must do: less work per replicate, not more. Measured from the engine's own `timings` attribute (microsecond resolution, so not a clock artefact), 13 cells x 3 arms. The replicate counts are the whole effect: median `reps_done` 17 pre-fix / **24 = budget exhausted** / 16 merged. **Correct signature: MORE replicates, each CHEAPER, identical score** — a blocked replicate never registers a hit on the best score, so the convergence rule never trips and the search burns its entire budget."
  ),
  list(
    old = "Not shown identical (different battery), but the mechanism matches: **re-measure before citing 1.43x.**",
    new = "**RE-MEASURED 2026-07-29 — mechanism CONFIRMED; the 1.43x is REATTRIBUTED, not retracted.** Three-arm battery (`bbcca1ba` pre-fix / `7685bf07` complement-only / `796a29d3` merged, all predating this fix so the defect is still live in arm 2), harness `dev/benchmarks/t364_*`. Decisive result: **arm 1's violating set and arm 2's complement-rooted set are the SAME SET** — on the T-364 test case the recorded 35/400 = 8.75% pre-fix violation rate reproduces to the digit, and arm 2's complement-only rate on the same case is also 35/400, on the same addition orders (arm 3: 0/400). An order that used to skip the constraint ran an effectively *unconstrained* search — cheap, converging early — and under complement-enforcement-alone the identical order runs a fully move-blocked one. That substitution is the wall gap, and it needs no accidental Wagner restarts: the old attribution was wrong twice over, since `AdditionTree()` never retried (`has_posthoc = FALSE`) and the Wagner build is ~1-2 ms, far too small to source a 1.43x. So the cost was real, and belongs to complement-WITHOUT-reroot — the configuration `355c4196` shipped — not to the merged fix and not to the constraint being honoured. Do not write \"retracted\": complement enforcement is not free for any producer that hands an unrerooted tree to `map_constraint_nodes`, which the fuse/sector/parallel paths listed above still do. **Prevalence, which decides how often this fires and which no user controls:** the canonical mask is always the tip-0-EXCLUDING side, so what matters is which side of the constraint holds `names(dataset)[1]` — i.e. matrix row order. Tip 0 in the large side (probability ~1 - k/n for a k-taxon group) gives **6-11.5% complement-rooted** (median 8.5%, reproducing this row's 9.5%); tip 0 in the small side gives **71-81%**. **Corollary: this row's \"roughly 1-in-11\" understates T-364 badly** — that is the common-geometry figure; on the minority geometry `AdditionTree()` violated up to 81% of the time."
  )
)

for (r in reps) {
  n <- length(gregexpr(r$old, raw, fixed = TRUE)[[1]])
  if (identical(gregexpr(r$old, raw, fixed = TRUE)[[1]][1], -1L)) {
    stop("PATTERN NOT FOUND (file changed under us?): ", substr(r$old, 1, 60))
  }
  if (n != 1L) stop("pattern matched ", n, " times, expected 1: ", substr(r$old, 1, 60))
  raw <- sub(r$old, r$new, raw, fixed = TRUE)
}

writeLines(raw, f, useBytes = FALSE)
cat("patched OK; both sentences replaced exactly once\n")
