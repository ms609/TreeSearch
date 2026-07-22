# Lever #6 — incremental cross-clip edge-set maintenance (L3b "Config C"), LAND

**Status:** SPEC, ready to build — but this is the BIGGER, correctness-critical lift. **Land #7 first.**
**Written 2026-07-14** from the per-move investigation (memory `tnt-per-move-kernel-gap`).

## 0. Mission (and the honest ceiling — read before committing effort)
Land an EXACT speedup by maintaining the directional edge set `E[D] = combine(prelim[D], up[D])`
**incrementally as the clip point sweeps a TBR pass** (Euler tour), instead of the current per-clip
from-scratch O(n) recompute (`compute_insertion_edge_sets`, src/ts_fitch.cpp).

- **Measured ceiling: ~1.2–1.5× wall** — the edge-set machinery is only ~17 ns of the ~51 ns/candidate
  (the abapprox A/B). This is a *large-N EW* lever and **orthogonal to 5432 reach**.
- **This is a substantial, correctness-critical kernel** (the suppress-node incremental transition). It
  is NOT a prototype-polish. **Gate the decision:** confirm a ~1.5× large-N EW wall win is worth the
  build before committing; if the C-vs-B saving turns out ≤ the B-vs-A candidate cost (see §4), it's a
  net loss — STOP.

## 1. Feasibility — NOW ESTABLISHED AT SCALE (this is why the lever reopened)
`dev/profiling/l3b-footprint-482.md` (worktree branch `worktree-agent-a095e0de2efe45407`): the per-clip
changed-view footprint **collapses with N** — mean **0.41 (75t) → 0.24 (240t) → 0.18 (482t)**, 84% of
clips under 0.3 at 482t, clearing the GO threshold. The earlier "L3b DEAD" verdict
(`dev/plans/2026-06-19-lever3-incremental-edgeset.md`, measured only at 37–75t where the footprint was
0.41–0.68) was **scale-limited** — cross-clip locality genuinely emerges at 482t. So incremental
maintenance is **feasible at mission scale.** The `-DTS_EDGESET_FOOTPRINT` instrumentation
(`ts_edgeset_footprint.h` + driver `l3b_footprint.R`) on that worktree branch is the validation harness;
reuse it. (Primary metric = **fp_ref** — each clip's divided tree vs the per-pass intact base tree,
order-independent = the true patch cost; fp_prev inflates ~2× under random clip order.)

## 2. The design — read `dev/plans/2026-06-19-lever3-incremental-edgeset.md` "Config C — locked kernel design"
That doc is the spec; summary:
- Visit clip edges in **DFS/Euler order** (config B = `order_clips_dfs` / `TS_TBR_DFS_CLIPORDER`, built on
  the old l3b worktree — may need reconstructing) so consecutive clips differ by ONE tree step.
- Maintain `prelim[]` (down — already incremental via `fitch_incremental_downpass`, nz→root) and `up[]`
  (view-from-above) as **persistent, incrementally-mutated** state as the clip boundary sweeps — instead
  of restore-to-full + recompute per clip.
- **Move rules:** descend p→c re-adds the chunk `subtree(p)\subtree(c)`; ascend inverts; **suppress-node**
  (clip parent `nx` removed, sibling `ns` joined to grandparent `nz`) is the single most bug-prone step —
  refresh `up[]`/`edge_set[]` at ns's new position + everything newly exposed.
- Re-materialise only the O(changed) `edge_set` entries per move; score this clip's regrafts immediately
  in place; do **NOT** cache all clips' arrays (that is O(n²) memory/copy and kills the win).
- Amortised **O(N) total** directional work per pass (vs O(N²)); each node exposed/hidden O(1)× over a
  full DFS.

## 3. Obstacles (from the design doc — plan for them)
- **Abandon restore-between-clips + the clip SHUFFLE.** `order_clips` randomises clip order to seed tight
  cutoffs early; DFS order is the opposite → cutoffs may loosen → MORE candidates scored. The B-vs-A
  kill-test measured **~1.2× candidate cost** (acceptable; inherited by both B and C).
- **Feature interactions:** gate the new path to plain EW (`use_directional && !has_na` + no
  sector_mask/constraint/tabu/pool); fall back to the current per-clip recompute whenever those are
  active. The reroot-completeness loop, collapsed skips, and TopoSnapshot all key off the
  restore-per-clip invariant — the incremental path is a NEW gated code path, not a replacement.

## 4. Correctness gate (MANDATORY — the suppress-node transition is the risk)
1. **Oracle assert (spine):** under an asserts-on build, assert the incrementally-maintained `edge_set`
   **equals the from-scratch `compute_insertion_edge_sets` result FOR EVERY CLIP, full array**, across
   all seeds/datasets, BEFORE trusting any wall number. Also `tbr_oracle.R` (0-improving) must still pass.
2. **Trajectory gate:** DFS order changes the RNG stream, so score + candidates_evaluated will NOT be
   byte-identical vs the shuffled baseline. So:
   - **C-vs-B** (config C incremental vs config B from-scratch, *identical* DFS clip order + RNG stream) ⇒
     bit-identical trajectory ⇒ **pure wall delta = the incremental saving** (the strong gate).
   - **B-vs-A** (DFS-order vs production shuffle) ⇒ the ~1.2× candidate-cost of dropping the shuffle
     (same final score, no worse candidates-per-improvement at matched wall).
   - **Ship test = C-vs-A** = (C-vs-B wall saving) − (B-vs-A candidate cost). If negative, STOP.

## 5. Measurement
C-vs-B (bit-identical) on Wortley/Zanol/Zhu locally + the **5432 subsample gradient on Hamilton** (the
~1.5× is a large-N effect; the small panel may under-show it). Re-measure the addressable ceiling at the
*shipped* `ratchetCycles` (the ceiling was originally measured at ratchet=12; it's lower now).

## 6. Build & standing constraints
Same as lever #7 (see `2026-07-14-lever7-scorer-monomorphization.md` §5): own worktree, named files only,
no push, never touch cpp-search/main, per-agent install (not load_all), stale-object ABI gotcha
(`CCACHE_DISABLE=1 --preclean`), Makevars PKG_CPPFLAGS + grep the build log, no hot-path `thread_local`,
r-conventions, verify anchors against the current tip. Note the old `claude/l3b-edgeset-timing` branch is
GONE — the a095 worktree reconstructed the footprint instrumentation from scratch; config-B
`order_clips_dfs` likely needs reconstructing too.

## 7. Deliverable
A gated incremental-TBR inner path + the oracle-asserted correctness + C-vs-B / C-vs-A wall on the panel
and 5432 subsamples + a durable `dev/profiling/l3b-land.md`. Worktree branch, NOT pushed — supervisor
merges. If the ship test (C-vs-A) is not clearly positive at 482t, record the numbers and STOP — a clean
negative is a valid result.
