# EW/NA bail-cutoff tighten — drop the +1 (`TS_TBR_LOOSE_CUTOFF`) — 2026-07-16

**Status:** LANDED default-tight (byte-identical, guaranteed non-negative). Local per-candidate
ns is **flat** (no measurable win on project970/project2668); the tighter bound is the correct one
regardless (drift already ships it) and may pay in the raw-TBR @482t plateau regime — flagged for
the EPYC A/B, not blocking.

## The lever

The EW/NA integer bail cutoff was `best_candidate - divided_length + 1`; scorers bail on
`extra_steps >= cutoff`. The `+1` lets an **equal-length, non-improving** candidate accumulate its
full `extra` before rejection (since `candidate == best` is not `< best`), i.e. it scans **every
block** — the plateau deep-scan flagged in `kernel-speed-baseline.md`. Dropping the `+1`
(slack 0) makes it bail the moment its running `extra` first reaches the current minimum.

Chosen by the post-packing lever-hunt workflow (`kernel-reroot-simd-hunt`) as the top pick, with
**3/3 adversarial verifiers confirming** byte-identity and non-wash — the strongest consensus of any
candidate lever.

## Why byte-identical (verified, and empirically gated)

Acceptance is strict `<` at every site (`best_candidate` updates only on a true improvement), so:
- a bailed candidate returns `>= cutoff = best-divided`, giving `candidate >= best_candidate` →
  rejected under `<` either way (no false accept, no lost move);
- a true improver (`extra < best-divided`) never reaches the tighter cutoff → computed exactly;
- an equal candidate bails one block earlier but is rejected identically.

The `+1` was **not** protecting `acceptEqual` (that path reads the recorded `best_candidate`
post-scan, never an individual scan return) and the **drift path already ships slack 0**
(`ts_drift.cpp`) — a live production instance. The x4 shared-bail breaks only when all four lanes
bust, so no improver is under-reported.

Sites tightened (`ts_tbr.cpp`): `spr_scan_plain_ew` template, SPR general, reroot x4, reroot scalar.
**Left `+1`:** `spr_scan_plain_iw` (cutoff dead on the IW path — IW bounds on the float
`best_candidate`), the IW-float `iw_cutoff`, and the `INT_MAX` init. Missing a site is still exact —
it only leaves savings on the table. Gated `TS_TBR_LOOSE_CUTOFF=1` reverts to `+1` for the A/B.

## Evidence

- **Byte-identity (must):** tight (default) vs loose on project970 EW + IW(k=10) + Sansom2010 NA —
  identical `score` + `n_candidates_evaluated` (38.9M / 44.8M / 30.6K) + `n_clips_tried` + final edge.
  (`scratchpad/cutoff_validate.R`.)
- **Local ns A/B (blocks-scanned is hardware-independent, so this is representative):**
  project970 loose 60.86 → tight 61.05 ns/cand; project2668 34.90 → 34.80 — **flat (noise).**
  The exact-equal-full-scan population is small: dedup (`seen_vp_hashes`) already skips equivalent
  "put-it-back" rerootings, and worse candidates already bail early under char cost-ordering.

## Disposition

Landed default-tight because it is **byte-identical and guaranteed non-negative** (fewer-or-equal
blocks scanned, always — unlike a prefetch hint it cannot regress), it is the correct tighter bound
matching the shipping drift path, and it may pay in the raw-TBR @482t plateau regime the workflow
flagged (project5432) even though it is a no-op at 157-196 tips. To be measured alongside the
prefetch lever on the Hamilton EPYC raw-TBR wall A/B; not blocking.
