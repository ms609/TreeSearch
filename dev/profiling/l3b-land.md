# L3b (lever #6) — incremental cross-clip edge-set maintenance: LANDED (Scheme 1)

**Date:** 2026-07-14 · **Branch:** `claude/lever6-edgeset` (off #7 tip `5d7e7c2c`) ·
**Status:** BUILT + oracle-validated + wall-positive at scale. Gated OFF by default
(`TS_L3B_INCREMENTAL`); supervisor decides default-on. NOT pushed.

## TL;DR

An EXACT (bit-identical) incremental replacement for the per-clip from-scratch
`compute_insertion_edge_sets` recompute, gated to the plain directional EW/IW TBR
regime. On real project5432 data it delivers a **raw-`tbr_search` wall speedup of
1.36–1.38× at 482t** (both single-rooting and the production `unrooted` reroot
config), rising from a wash at ≤75t — exactly the predicted large-N ceiling. Score
and trajectory are bit-identical to production; the per-clip oracle is clean on real
multi-state/ambiguous data including the ratchet `upweight_mask`/`collapsed` context.

**But the end-to-end `MaximizeParsimony` win is a WASH** (measured, 1.01× at both 120t
and 240t), because the sectorial ~30% of wall is gated out AND ratchet's short passes
amortize the per-pass O(n_node) base recompute poorly. **So the honest verdict is: an
exact, correct, large-N *raw-TBR* wall lever that does NOT (as built) translate into a
net *mission-wall* win** in full `MaximizeParsimony` at the sizes measured. Per the
spec ("a clean STOP with recorded numbers is a valid result"), this is a recorded
near-negative on the mission metric: **land default-OFF as a validated building block;
do NOT flip default-ON.** The unlock is the base-incremental-update follow-on (removes
the per-pass overhead that eats the win — see caveats). This is the measurement-first
outcome the spec asked for.

## Why Scheme 1, not the Euler "Config C" the design doc locked

The spec (`dev/plans/2026-07-14-lever6-incremental-edgeset-land.md`) names "Config C"
= the Euler/DFS cross-clip incremental kernel from
`dev/plans/2026-06-19-lever3-incremental-edgeset.md`. I built **Scheme 1**
(patch-from-full-per-clip) instead, deliberately and with advisor sign-off. The
reasoning (recorded here so this does not read as quietly skipping Config C):

- **The evidence that REOPENED the lever supports Scheme 1, not Euler.** The 482t
  feasibility (`dev/profiling/l3b-footprint-482.md`) measures `fp_ref` — the per-clip
  *divided-vs-intact-base* footprint (0.18 mean at 482t). That is exactly Scheme 1's
  patch cost. The Euler route's actual cost-driver is the *per-descend-step delta*,
  which was measured ≈1.14–1.24× the footprint at 37–75t — that was the small-N
  *kill* — and was **never re-measured at 482t**. Building Config C would be building
  the harder kernel on the weaker evidence.
- **The TBR clip loop is already Scheme-1-native.** It restores the tree to a fixed
  per-pass base after *every* clip (`TS_REVERT_CHECK`, ts_tbr.cpp Phase 2), so each
  clip is an independent perturbation of one intact tree. Scheme 1 fits without
  ripping out restore-per-clip or the clip shuffle.
- **Scheme 1 keeps the shuffle ⇒ bit-identical trajectory vs production (config A).**
  So the ship test is a clean pure-wall C-vs-A with NO B-vs-A candidate cost — a
  stronger gate than the spec's C-vs-A-minus-B-vs-A. Euler would have paid the ~1.2×
  DFS-order candidate cost.
- **Lower bug surface.** One patch + undo per clip; no sweeping suppress-node state
  machine carried across a DFS tour.

Decision tree (advisor): Scheme 1 wins wall at scale → ship (done — it does). Had it
been a node-count win but wall-null, the next step would have been to re-measure the
Euler descend-delta at 482t before attempting Config C. Had it been node-marginal,
the whole lever would be dead at scale — STOP. Neither fallback was needed.

## Design (what was built)

`patch_insertion_edge_sets` (src/ts_fitch.cpp) patches the working `up[]`/`edge_set[]`
(== a pristine per-pass intact base at entry) to the current divided (clipped) tree,
touching only the changed frontier around the clip, and records every touched node id
in a `changed` list.

- **Shared combine kernel.** The from-scratch and incremental paths share one
  `ts_fitch_combine` (factored out of `compute_insertion_edge_sets`), so patched words
  are byte-identical to from-scratch — the oracle equality rests on this.
- **Per-pass base.** At each pass top (tree intact, restored there after every clip),
  `compute_insertion_edge_sets` computes the pristine `edge_set_base`/`up_base` once
  (O(n_node), negligible vs the O(n_node)×n_clips it replaces), then memcpy → working
  buffers. Recomputed every pass, so it survives accepted moves and inter-pass reroots.
- **Discovery = worklist incremental up-pass (correct-by-construction).** Seeds:
  (A) path nodes nz→root have changed prelim but *invariant* up, so only their
  `edge_set` is refreshed; (B) each path node's off-path sibling + nz's two divided
  children {ns, W} have changed up. Propagate down each seed, recompute `up[D]`,
  compare word-for-word to `up_base` (the combine-attenuation frontier), stop on match.
  Every up-change traces back to a prelim-change (path) or the suppress-node topology
  move via the up[parent] chain, so discovery is complete; a tree's single-parent
  structure + disjoint seed subtrees make it single-visit (no visited set).
- **Changed-list = undo log.** After the clip's scan, memcpy base→buffer over `changed`
  restores the base for the next clip. No separate save buffer.
- **Gate** (`l3b_active`, ts_tbr.cpp): `TS_L3B_INCREMENTAL` set AND `use_directional`
  (EW/IW, no NA) AND no sector_mask / constraint / tabu / pool / b2 probe AND
  n_tip ≥ 4. The `unrooted` reroot loop (`do_reroot`) IS supported — see below. Any
  other config falls back to the from-scratch recompute. The patched `edge_set_buf`
  is the SAME buffer #7's `spr_scan_plain_ew/iw` templates, the general loop, and the
  reroot `vroot_cache` all read.

### `do_reroot` is supported (and is not a footprint penalty)

Neither the spec nor the advisor asked to exclude the unrooted reroot path, and the
unconstrained production driven search runs `unrooted=true` (`TBRParams.unrooted`
defaults true; ts_driven.cpp leaves it), so excluding it would have made the lever a
no-op in production. It is compatible: the tree is still restored to the pass-start
base after every clip; reroots happen only *between* passes (base recomputed at each
pass top after the reroot's `full_rescore`); and `try_root_edge_moves` uses its own
buffers, never `edge_set_buf`. `do_reroot` only relaxes the smaller-side clip filter,
so it also processes larger clips — but a larger clip leaves a *smaller* divided tree,
so its footprint is equal-or-lower (measured: fp/mean_changed for `unrooted=1` ≤
`unrooted=0` at every size). A footprint-threshold fallback for a p99 tail was
therefore **not needed** (per-clip patch cost is ~constant in N: ~50–58 nodes at
40–160t while from-scratch grows linearly 71→310).

## Correctness (the mandatory gate — fully passed)

Per-clip **oracle** (ts_tbr.cpp): under `TS_L3B_ORACLE` or an asserts-on build, every
clip re-runs `compute_insertion_edge_sets` into a reference buffer and asserts
`edge_set_buf` equals it word-for-word for every divided in-tree non-root node;
`Rf_error` aborts on the first mismatch. Results — **zero mismatches** across:

- **Random EW** (validate_l3b.R): 4 sizes (12/25/50/80t) × 8 seeds × `unrooted`∈{F,T}
  = 64 configs. C-vs-A **bit-identical** (score + full edge matrix): 0/64 differ.
- **Real multi-state / ambiguous data** (the representative cases): panel
  Wortley2006/Zanol2014/Zhu2013 + project5432 subsampled to 120t & 240t, both
  `unrooted`. Oracle mismatch = 0 everywhere; A==C score everywhere.
- **Canonical completeness oracles**: `tbr_oracle.R` (150 trees, do_reroot path,
  EW) and `tbr_oracle_iw.R` (80 trees, IW) — 0 missing moves with l3b+oracle on.

Coverage includes near-root and `nz==root` clips (handled by the root degree-2 special
case), single-rooting and `do_reroot`, EW and IW.

**Harness gotcha found & fixed:** an unseeded `AdditionTree(dat)` start made every
process build a *different* start (R's RNG is time/PID-seeded pre-`set.seed`), so raw
A-vs-C score comparison looked like a bug (488 vs 491). Seeding before the start made
A and C bit-identical and stable. Always seed the start when comparing A-vs-C. The
per-clip oracle (not score-equality) is the true bit-identity gate on nondeterministic
real-data runs.

## Wall (C-vs-A, raw `tbr_search` to convergence, fixed start ⇒ identical trajectory)

`-O2`, this machine, min-of-reps ms. C = `TS_L3B_INCREMENTAL=1`, A = default. Scores
identical in every row.

| dataset            | n_tip | A (ms) | C (ms) | speedup | unrooted |
|--------------------|-------|--------|--------|---------|----------|
| Wortley2006        | 37    | 3.6    | 3.7    | 0.97×   | 0        |
| Zanol2014          | 74    | 26.9   | 25.4   | 1.06×   | 0        |
| Zhu2013            | 75    | 14.5   | 13.1   | 1.11×   | 0        |
| project5432 (sub)  | 120   | 66.3   | 60.3   | 1.10×   | 0        |
| project5432 (sub)  | 240   | 338.1  | 243.9  | 1.39×   | 0        |
| **project5432**    | **482** | **3005.7** | **2172.2** | **1.38×** | 0 |
| **project5432**    | **482** | **2899.2** | **2131.7** | **1.36×** | 1 (reroot) |

The win is a wash at ≤75t (small-N regime: footprint ≈ 1.0) and grows to ~1.4× at
mission scale — matching the honest ceiling (edge-set machinery ~17 of ~51 ns/candidate
⇒ ~1.5× is the raw-TBR bound) and the `fp_ref` trend (mine falls 0.55→0.34→0.19 over
40→80→160t, tracking `l3b-footprint-482.md`). It holds in the production `unrooted=1`
config (1.36×).

### End-to-end `MaximizeParsimony` (ratchet + sector + driven) — DILUTED

The raw-`tbr_search` number above is NOT the end-to-end mission wall. Two effects
dilute it in a full search:

1. **Coverage** (confirmed by `TS_L3B_STATS`): l3b fires in the driven-TBR and
   **ratchet** phases (ratchet's inner `tbr_search` is l3b-eligible — no
   sector_mask, `upweight_mask` doesn't gate it). But the **sectorial workhorse
   (~30% of wall) is gated OUT** (`sector_mask != null`), so ~30% of the search sees
   no speedup.
2. **Per-pass base-recompute amortization.** The pristine base is recomputed once
   per pass = one from-scratch-equivalent (O(n_node)). A pass ends at its first
   accept (ts_tbr.cpp:2790 `break`), so a pass of k clips nets a win iff
   **k > 1/(1−fp_ref)** (≈1.2 at 482t fp=0.18, but ≈5 at 75t fp=0.8). Ratchet's
   short high-fp perturbation-recovery passes (10–40 clips each at 75t, many of them)
   therefore amortize the O(n_node) base poorly — at small N this can erase the win.

**Oracle re-confirmed in this context** (Zhu2013 EW, ratchet+sector, `TS_L3B_ORACLE`):
0 mismatches across the whole full-search run — the `upweight_mask`/`collapsed`/
interleaved-accept contexts my `ts_tbr_diagnostics` validation never hit are clean.

Measured full `MaximizeParsimony` C-vs-A (project5432, EW, deterministic start/seed,
score identical in every row):

| N   | maxReplicates | A (s) | C (s) | end-to-end speedup |
|-----|---------------|-------|-------|--------------------|
| 120 | 3             | 19.03 | 18.81 | 1.01× (wash)       |
| 240 | 2             | 59.84 | 59.51 | 1.01× (wash)       |

**Both are a wash — even at 240t, where the raw-TBR lever is 1.39×.**

Why (discriminated, not guessed): two hypotheses could explain the wash — (a)
per-pass base-recompute amortization across ratchet's many short passes, or (b) the
per-clip win being structurally small in-search because ratchet perturbs to
higher-homoplasy trees with a larger footprint. The `TS_L3B_STATS` in-search
`fp_changed` at 240t MP settles it: **0.244** (changed/clip 77.2, edges/clip 316.4) —
essentially the *clean* 240t footprint (`l3b-footprint-482.md` reports 0.243), so
(b) is refuted: **the per-clip win is intact.** The culprit is (a) + the sector
exclusion: sector-off alone predicts ~1.24× (`1/(0.3 + 0.7/1.39)`); the drop from
1.24× to 1.01× is the per-pass O(n_node) base recompute, which MP's short
ratchet-recovery passes amortize far worse than my single Wagner→convergence raw
measurement (dominated by long early passes). So the recompute overhead, not the
patch, eats the win in the full search.

## Honest caveats / open items

- **Raw-`tbr_search` win is real; end-to-end is diluted (MEASURED, not estimated).**
  1.36–1.38× is a single TBR search to convergence. Full `MaximizeParsimony` is a
  wash at 120t and diluted at scale — see the End-to-end table above and its two
  causes (sector ~30% gated out; per-pass base-recompute amortization). The honest
  framing: **this is a large-N raw-TBR wall lever, not (yet) a net mission-wall
  lever.** To make it mission-relevant, in rough priority:
  1. **Incrementally update the base after an accept** instead of the per-pass full
     `compute_insertion_edge_sets`. The in-search fp=0.244 measurement confirms the
     per-clip patch is intact, so this recompute IS the dominant end-to-end overhead;
     the accept mutates only a local region, so an O(footprint) base update should
     recover most of the sector-diluted ~1.24× ceiling — the highest-value follow-on.
     **Caveat: this is itself a correctness-critical kernel** — incremental edge-set
     maintenance across a *TBR accept* (two dirty prelim paths, nz→root and nx→root,
     plus possible reroot) is harder than the clip patch and needs its own oracle.
  2. **Gate to large N only** (e.g. n_tip ≥ ~150) so small-N searches, where the win
     is a wash/loss, keep the from-scratch path. Cheap, safe, one line.
  3. **Extend to the sector path** (the ~30% currently excluded) — larger follow-on.
  - **Hamilton full-search 482t: NOT worth running yet.** The end-to-end wash at
    120t AND 240t means a full `MaximizeParsimony` 482t run would very likely also be
    a wash (the dilution is structural, not size-threshold), so the Hamilton deploy
    cost is not justified until the base-incremental-update (1) lands and re-opens a
    plausible end-to-end win. Re-run it then, at the shipped ratchet, multi-seed.
- **Memory:** +2 arrays (`edge_set_base`, `up_base`, each n_node×total_words) beyond
  the existing `edge_set_buf`/`edge_set_up`. Fine at 482t; note for very large N.
- **Orthogonal to project5432 *reach*** (cold-start basin capture) — a general large-N
  wall lever, as the feasibility doc stated.

## Reproduce

```
# build (headers touched -> avoid stale .o):
CCACHE_DISABLE=1 R CMD INSTALL --preclean --no-multiarch --no-docs --library=.agent-l3b .
# correctness (oracle aborts on any mismatch):
TS_L3B_INCREMENTAL=1 TS_L3B_ORACLE=1 Rscript dev/benchmarks/tbr_oracle.R 150 12
TS_L3B_INCREMENTAL=1 TS_L3B_ORACLE=1 Rscript dev/benchmarks/tbr_oracle_iw.R 80 12 3
# wall + footprint: see scratchpad real_wall.R / fp_one.R (TS_L3B_STATS=1 for fp_changed).
```
Flags: `TS_L3B_INCREMENTAL` (on), `TS_L3B_ORACLE` (per-clip assert; auto in asserts-on
builds), `TS_L3B_STATS` (mean patched-node fraction, oracle-populated).

## Files (worktree branch `claude/lever6-edgeset`, NOT pushed)

- `src/ts_fitch.cpp` — `ts_fitch_combine` (factored), `patch_insertion_edge_sets`.
- `src/ts_fitch.h` — declaration.
- `src/ts_tbr.cpp` — `l3b_active` gate, per-pass base, per-clip patch + oracle + restore, stats.
