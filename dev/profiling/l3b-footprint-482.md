# L3b edge-set changed-view footprint at LARGE N (project5432)

**Date:** 2026-07-13
**Lever:** L3b (#6) — maintain the exact directional insertion edge set
`E[D] = combine(prelim[D], up[D])` INCREMENTALLY across clips (patch the views
as the clip point moves) instead of the current per-clip from-scratch O(n)
recompute (`compute_insertion_edge_sets`, `src/ts_fitch.cpp:521`).

**Question (the explicit reopen condition):** L3b was measured DEAD at 37–75t
(footprint 41–68% of edges; GO needs `<~0.3`), but the footprint FELL with
tree size (0.68→0.46→0.41 at 37→74→75t) and the finding named large N as its
reopen condition. project5432 = 482 tips = 6× past anything measured. Does the
footprint cross `<0.3` at 482t?

## Verdict

**FEASIBLE at scale — the "dead-at-scale" finding is refuted.** The per-clip
changed-view footprint collapses with tree size. At 482t the **amortized
(mean) `fp_ref` = 0.177** — the wall-relevant statistic, since total patch work
is set by the mean, not the median — which clears the 0.3 GO line and sits well
under the 0.41–0.68 small-N dead regime. The median is 0.05, but that is the
least wall-relevant number (the median clip is a cherry; see the skew note
below), so the verdict rests on the mean. 84% of clips change fewer than 30% of
the tree's directional edge views. The June "DEAD / no cross-clip locality"
conclusion is an artefact of small N and does **not** hold at 482t; the reopen
condition is met.

**Feasibility bracket [0.18, 0.33].** Because each clip is an independent
perturbation of a fixed per-pass base tree, a real L3b must patch full→divided
and then get *back* to full for the next clip. If the restore is a cheap
undo-log memcpy and the changed frontier is discovered by the patched up-pass
itself, the combine work stays at ~mean(`fp_ref`) = **0.18**. If the round-trip
write cost counts, it is ~2× — corroborated directly: mean(`fp_prev`), the
consecutive-clip "slide" cost, = **0.329** at 482t (1.86× mean `fp_ref`). So the
honest amortized footprint is **[0.18, 0.33]**, straddling the 0.3 line. The
provenance of the exact 0.3 threshold (one-trip vs round-trip) is gone from git,
so bracketing is the honest call rather than a spurious re-derivation. Either
way this is a decisive step down from the 0.41–0.68 that condemned it at
small N.

## Method

- **Metric `fp_ref` (primary):** for each clip in a TBR pass, the fraction of
  surviving non-root in-tree nodes `D` whose `E[D]` in the clip's DIVIDED tree
  differs (word-for-word) from the intact base tree's `E[D]`. The base-tree
  reference is recomputed once per pass and each clip is diffed against it.
  This is exactly the per-clip patch cost the L3b incremental lever would pay
  instead of the full recompute — order-independent, self-contained per clip
  (no dependence on accepts), so no warm-up-to-convergence is needed.
  `fp_ref < ~0.3` ⇒ patch beats recompute (GO).
- **Metric `fp_prev` (secondary):** same, but diffed against the *previous
  recorded clip's* divided tree. Under the default RANDOM clip order this is
  the union of two independent perturbations, so it runs ~2× higher than
  `fp_ref` and is a pessimistic bound — reported only as a cross-check.
- **Instrumentation:** flag-gated `-DTS_EDGESET_FOOTPRINT`
  (`src/ts_edgeset_footprint.h`; reference-capture + record hooks in
  `src/ts_tbr.cpp`; probe export `ts_tbr_fp_probe` in `src/ts_rcpp.cpp`,
  registered in `src/TreeSearch-init.c`). Flag OFF ⇒ hooks compile out and the
  probe returns empty vectors (guardrail: the driver aborts on an empty result,
  so a silent no-compile can't masquerade as a plausible number).
- **Config (identical across all N — the cross-N trend is the verdict):**
  pure equal-weights Fitch (inapplicable `-` recoded to `?` so `has_na=false`
  and the exact directional edge-set path is exercised); Wagner-addition start
  (`AdditionTree`); `tbr_search` with `unrooted=false` (no reroot tip-sweep, so
  482t stays light) and `maxChanges=0` (inner convergence); default RANDOM clip
  order; production smaller-side clip filter active (`clip_size > n/2` skipped)
  — so the measured clips are exactly the set TS actually scores.
- **Subsampling:** project5432 (482t × 189 chars) subsampled to N random taxa
  (`set.seed(N)`). Driver: `dev/profiling/drivers/l3b_footprint.R`.
  Raw samples + summary: `dev/profiling/l3b-footprint-482.rds`.

## Results

`fp_ref` = fraction of surviving non-root nodes whose `E[D]` changes vs the
intact base tree, per clip. Denominator `edges_med` = median surviving
non-root nodes.

| N   | clips  | edges | fp_ref median | fp_ref mean | p75  | p90  | frac<0.3 | fp_prev median |
|-----|--------|-------|---------------|-------------|------|------|----------|----------------|
| 75  |   665  |  144  | **0.303**     | 0.406       | 0.71 | 0.97 | 0.49     | 0.634          |
| 120 |  1517  |  234  | **0.254**     | 0.450       | 0.97 | 1.00 | 0.53     | 0.949          |
| 240 |  3989  |  474  | **0.078**     | 0.243       | 0.32 | 0.79 | 0.74     | 0.309          |
| 482 | 12127  |  956  | **0.049**     | 0.177       | 0.20 | 0.68 | 0.84     | 0.186          |

**Trend:** `fp_ref` falls monotonically and steeply with N. The amortized mean
is at/above the 0.3 line at 75–120t (0.41/0.45 — the small-N "dead" regime) and
falls to **0.24 at 240t and 0.18 at 482t**; the median (cherry-dominated) falls
faster, from 0.30 to 0.05.

- On the wall-relevant **mean**, the GO line is cleared by 482t (0.177) and is
  borderline at 240t (0.24). On the median it is crossed by 240t. Even the
  pessimistic slide cost mean(`fp_prev`) = 0.33 at 482t is only just above 0.3.
- **Why it falls:** a random smaller-side clip removes a subtree that is a
  shrinking fraction of a growing tree, and the directional Fitch combine
  (intersect-else-union) *attenuates* the perturbation as it propagates — so at
  large N the changed views are tightly localized around the clip.
- **Right-skewed tail (honest):** ~16% of clips at 482t still change most views
  (p90 ≈ 0.68, p99 ≈ 1.0) — the large smaller-side subtrees / near-root clips.
  A real L3b implementation would patch cheaply for the 84% majority and could
  fall back to full recompute above a footprint threshold, keeping the
  amortized cost at or below the 0.18 mean.
- Run-to-run RNG noise in the search trajectory moves the 482t median in
  0.05–0.06; both replicates sit an order of magnitude below 0.3. Verdict is
  robust.

## Honest ceiling (do not over-read a GO)

This measurement decides **FEASIBILITY, not magnitude.** Even with fp_ref ≪ 0.3,
L3b is bounded to **~1.5× wall at best**: the abapprox A/B showed the edge-set
machinery is only ~17 ns of the ~51 ns/candidate, so eliminating the recompute
cannot exceed that share. The bigger per-move lever is elsewhere (a plain-EW
fast path). L3b is also **ORTHOGONAL to project5432 *reach*** (cold-start basin
capture) — a feasible L3b is a general large-N *wall* lever, not a 1943-reach
lever.

**Net:** L3b moves from "DEAD (small-N footprint)" to "FEASIBLE at scale, modest
(~≤1.5× wall) lever — build only if a large-N wall win is prioritized." The
per-clip from-scratch recompute is genuinely wasteful at 482t; incremental
patching would recover most of it. Whether the ~≤1.5× is worth the
implementation lift (incremental up-pass maintenance + a footprint-threshold
fallback) is a prioritization call, not a feasibility one.

## Reproduce

```
# from a worktree on cpp-search, with the footprint flag in src/Makevars.win:
Rscript -e 'Rcpp::compileAttributes(".")'
R CMD INSTALL --library=.agent-l3b --no-docs --no-multiarch .
Rscript dev/profiling/drivers/l3b_footprint.R           # writes the RDS
```
Guardrail: confirm the build log shows `-DTS_EDGESET_FOOTPRINT` and the driver
does not abort on an empty `fp_ref`.
