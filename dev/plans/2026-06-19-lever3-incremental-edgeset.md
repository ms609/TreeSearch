# Lever 3 — incremental / amortized insertion edge-sets across clips

**Status:** BUILDING (worktree `claude/l3b-edgeset-timing`). Ceiling MEASURED —
24–33% of search compute at ratchet=12 (re-measure at shipped 6 pending) —
which OVERTURNS the lean-NO-GO below. Config B (DFS-clip-order kill-test) built;
config C (cross-clip incremental) gated on B's verdict. See "MEASURED" section.
**Owner:** this session (task #43). **Date:** 2026-06-19.
**Target:** `compute_insertion_edge_sets` (src/ts_fitch.cpp:478), the 27 % full-EW
CPU hotspot (findings T-P5a). The per-iteration deficit vs TNT.

## What is already banked / ruled out

| Lever | Outcome | Note |
|---|---|---|
| L1 skip zero-fill | **MERGED cpp-search 00d73d6a** | Zanol −16.4 %, sum −9.4 % (T-P5e/e2). Bit-identical re-gated on main checkout (EW 6-run + NA-serial 4-seed + 276 kernel tests). perclip orphaned (subset, nothing to redeem). |
| L2 vectorize combine | **AT-LIMIT** | −1.2 % in noise; n_states tiny, bandwidth-bound (T-P5f). |
| L3a fuse two sweeps | **AT-LIMIT (reverted)** | bit-identical, but +0.2/0.4 % (noise) — `up[]` (~10 KB) is cache-resident at mission sizes, so fusing saves no memory traffic (T-P5i). |
| **L3b cross-clip amortization** | **this doc — leans NO-GO** | the only remaining route; T-P5h shows TS rate flat in N ⇒ constant-factor, not asymptotic. Within-clip kernel now maxed (L1 banked, L2/L3a AT-LIMIT). |

## The proof that forces the design (why nothing lighter works)

Per clip of subtree `c` (parent `p`, sib, grandparent `gp`), we need the
divided-tree "view from edge above D", `edge_set_e[D] = combine(prelim_e[D],
up_e[D])`, for **every** regraft edge D (≈ O(N) of them).

`up_e[D]` = states of `(tree − subtree(c) − subtree_e(D))`.

- **D is an ancestor of c** (path c→gp→…→root): c ⊆ subtree(D), so c is *inside*
  D's own subtree, never in D's complement-of-view ⇒ `up_e[D] = up_full[D]`
  **unchanged**. These are O(depth) nodes.
- **D is NOT an ancestor of c** (≈ O(N) nodes): c lies in D's complement ⇒
  removing c changes `up_e[D]`. **Changed.**

So a once-per-pass full-tree directional precompute reuses only the clip's
ancestor-path; **the per-clip directional pass is irreducibly O(N)**. There is
no within-clip or full-tree-precompute shortcut (verified three ways). `prelim_e`
*is* already maintained incrementally (`fitch_incremental_downpass`, nz→root);
the cost is the up-pass + edge-set pass.

Combine (intersect-else-union) is **not invertible**, so c's contribution cannot
be "subtracted" from `up_full[D]`; it must be recomputed.

## The only route: DFS-order cross-clip incrementality (Goloboff 1996)

Visit clip edges in a **DFS / Euler order** so consecutive clips differ by one
tree step. Maintain `prelim` and `up` as a **persistent, incrementally-mutated**
state as the clip boundary sweeps the tree, instead of restore-to-full +
recompute per clip. Each subtree's directional messages are recomputed O(1)
times over a full DFS ⇒ **O(N) total directional work per pass, not O(N²)**.

### What this costs us architecturally (the obstacles)
1. **Abandon restore-between-clips.** Phase 2 (ts_tbr.cpp:1818,
   `restore_prealloc_undo`+`spr_unclip`+saved_postorder) restores to the
   pass-start tree after every clip. Cross-clip incrementality needs the clip
   point to *move* incrementally, carrying directional state with it.
2. **Abandon the clip SHUFFLE.** `order_clips` (ts_tbr.cpp:1063,1381) randomizes
   clip order — it exists to seed early candidates from across the tree so the
   bounded scorer gets a tight cutoff EARLY (ts_tbr.cpp:1489 partial-shuffle of
   main_edges serves the same end at the regraft level). DFS order is the
   opposite of "spread across the tree", so cutoffs may loosen ⇒ MORE candidates
   scored (efficiency loss) partly offsetting the directional saving.
3. **Interaction with NA / IW / sector / constraint / tabu / reroot paths.**
   The directional edge-set is used by EW and IW (`use_directional`); NA keeps
   union-of-finals. The reroot-completeness loop, sector_mask, collapsed-edge
   skips, and TopoSnapshot all key off the restore-per-clip invariant. A
   cross-clip-incremental TBR is a *new code path*, gated; it must fall back to
   the current per-clip recompute whenever those features are active.

### Why it is NOT a clear win (the honest risk)
- Saving = directional pass goes O(N²)→O(N) per pass. **But** (a) L1 already
  removed the zero-fill (the largest constant within that O(N²)); (b) dropping
  the shuffle may cost candidates (obstacle 2); (c) the combine work that
  remains is bandwidth-bound (L2 finding), so even halving node-visits may not
  halve wall. Net win is **measurement-decidable only**.

## GO / NO-GO criteria (both required before implementing)

1. **Advisor review** of this incremental-directional scheme (advisor overloaded
   on 2026-06-19 — retry). Specifically: is the DFS-incremental up-pass correct
   under suppress-node TBR, and is there a cheaper amortization I have missed?
2. **N-scaling confirmation, W-controlled — SUBSTANTIALLY ANSWERED, leans NO-GO
   (T-P5h, 2026-06-19).** Hamilton 64-bit `framing_64bit.csv` (job 17528864):
   TNT side FAILED (NA — re-run owed), but the **TS-side `ts_rate` is FLAT in N**
   (37t≈13.9, 55t≈13.4, 74t≈12.4, 75t≈20.1, 78t≈20.9 — varies with W not tips).
   TS pays no per-candidate O(N²) penalty (directional pass is O(N)÷O(N)
   cand/clip = O(1)/cand, already amortized). The local "throughput grows with N"
   was TNT-side (TNT speeds up per-candidate with N). ⇒ L3b attacks a
   **constant-factor** share, same order as L1/L3a, at high restructure risk.
   **Default decision: NO-GO.** Bank L1 (+measure L3a). REVISIT only if a
   SUCCESSFUL 64-bit TNT re-run shows tnt_rate asymptotically diverging from
   ts_rate in a way only cross-clip amortization could match.

## Fallback if NO-GO
L1 (zero-fill skip, −9.4..−16.4 %) + optional L3a (fuse) is the within-clip
ceiling. Combined with the ratchet recipe (T-P5d: `ratchetCycles` 12→6,
−20..38 % wall) these are the bankable per-iteration + recipe wins; the residual
per-candidate gap would then be declared constant-factor / near-limit and the
mission focus shifts to recipe composition (tasks #39/#40).

## Implementation sketch (only if GO) — isolated worktree, gated flag
- New `tbr_search` inner path behind a param flag (default off), active only for
  plain EW/IW search (no NA/sector/constraint/tabu/pool), large n_tip.
- Root the tree; Euler-tour the clip edges. Maintain `prelim[]` (down) and
  `up[]` (directional) as the clip boundary moves; patch only the exposed
  region per step. Score regrafts from the maintained `edge_set`.
- **Gate = bit-identical** final score AND `candidates_evaluated` vs the current
  path on {Wortley,Zhu,Zanol}×seeds (ab.R / ab_compare.R) — if the shuffle is
  dropped, candidates_evaluated WILL differ, so the gate becomes "same final
  score + no worse candidates-per-improvement at matched wall", a weaker gate
  that itself needs sign-off.

---

## MEASURED (2026-06-19) — ceiling overturns the lean-NO-GO; A/B/C build underway

Built bucketed `std::chrono` timing into `compute_insertion_edge_sets` (gated by
`-DTS_EDGESET_TIMING`, header `src/ts_edgeset_timing.h`), split by call site +
`sector_mask`, with a cumulative `driven_search` denominator, read out via
`ts_edgeset_timing_report()`.  Ran the FULL default-preset `MaximizeParsimony`
(single-thread, -O2), NOT the isolated clip-loop.  Driver: `l3b_ceiling.R`.

**Addressable ceiling = `tbr_plain` / total search compute** (the maskless
plain-clip-loop edge-set recompute — ratchet + plain TBR; the only surface L3b's
DFS restructure can touch):

| dataset | tips | ceiling (ratchet=12) |
|---|---|---|
| Wortley2006 | 37 | 24% |
| Zhu2013 | 75 | 29% |
| Zanol2014 | 74 | 33% |

- **Not low-single-digits** ⇒ the kill threshold is NOT met ⇒ the lean-NO-GO is
  overturned.  This is exactly the empirical answer the user asked for.
- The full-search figure **matches** the isolated clip-loop's 27–31% (T-P5a)
  because the default preset is clip-loop-dominated (`sector_ms`/`tbr_sectormask`
  ≈ 0 at these sizes).  **The findings.md:42 `~2%` dilution note was WRONG** for
  the default preset — superseded by this honest number.
- Ceiling **rises with tip count** (24→33%), softening T-P5h's "flat in N /
  constant-factor" framing: L3b's prize holds or grows with tree size.
- The ceiling is an **upper bound** (advisor): realisable is a fraction of it
  (L1 already banked the zero-fill; combine is bandwidth-bound (L2); see below).

### Build design — three gated configs (advisor 2026-06-19)
- **A** = production (shuffled clip order, from-scratch per clip).
- **B** = from-scratch, **DFS clip order** (`TS_TBR_DFS_CLIPORDER=1`; no rng
  consumed).  Built (`order_clips_dfs`, ts_tbr.cpp).
- **C** = incremental directional maintenance, DFS clip order (the hard build —
  not yet written).
- **C-vs-B** = identical clip order + identical RNG stream ⇒ bit-identical
  trajectory ⇒ pure wall-clock delta = the **strong-gated actual prize**.
- **B-vs-A** = the shuffle-drop trajectory cost (`order_clips` consumes rng, so
  dropping it shifts the whole downstream stream incl. the within-clip shuffle —
  hence the A/B/C split, NOT a direct C-vs-A).  **B is the cheap KILL-TEST:** if
  DFS order alone worsens score / candidates-per-rep, even a perfect C inherits
  it ⇒ STOP before building C.  Driver: `l3b_ab.R`.

### Kill-test result (B vs A, `l3b_ab.R`, ratchet=6, 8 reps, 3 datasets × 3 seeds)
**PASSED ⇒ GO on C.**  DFS clip order is NOT load-bearing:
- Final scores statistically equivalent — `dScore` ∈ {−1,0,+1}, no systematic
  loss (Zhu2013: B *better* in 2/3 seeds; Wortley: identical; Zanol: +1 in 2/3).
- candidates-per-rep ratio B/A ≈ 1.2 median (range 0.67–1.44, noisy) — a small,
  acceptable cost, NOT the ≫1 that would mean the shuffle is doing real work.
- The ~1.2× candidate cost is inherited by BOTH B and C, so C-vs-B isolates the
  incremental saving cleanly; the SHIP test is C-vs-A =
  (C-vs-B wall saving) − (B-vs-A ~20% candidate-cost).

### Honest caveats to carry into C
- Re-measure ceiling at **shipped ratchet=6** (worktree flip applied) — likely
  ~18–25%, since `tbr_plain` lives mostly in ratchet.
- **Realisable < ceiling, possibly well under:** incremental cuts node visits,
  but the combine is bandwidth-bound, so wall savings may lag the visit
  reduction.  C-vs-B measures the truth — do not promise the full ceiling.

### Config C — locked kernel design (build this; oracle-assert every clip)
Proven (3 ways) that per-clip directional work is irreducibly O(N) — combine is
non-invertible, the changed up-set is O(N), up-arrays are globally different per
clip so cannot be cached/shared. The ONLY sub-O(N²)/pass route is an **Euler-tour
clip walk with in-place incremental state** (no restore-between-clips).

**Maintain as the clip edge walks a recursive DFS Euler tour of the tree:**
- `prelim[]` (down) — already incremental via `fitch_incremental_downpass`.
- `up[]` (view-from-above) — mutate in place on each clip-edge move.
- `edge_set[D] = combine(prelim[D], up[D])` — re-materialise only the O(changed)
  entries touched by the move; score this clip's regrafts immediately (in place),
  then move on (do NOT cache all clips' arrays — that is O(N²) memory/copy and
  kills the win).

**Move rules** (clip edge = (parent, child); Euler tour enters child on descend,
returns on ascend):
- **Descend** parent p → child c: the divided tree re-adds the chunk
  `subtree(p) \ subtree(c)` (p itself + p's other child's subtree). Propagate the
  changed sibling-down-message at the divergence (nz/nx suppression point) DOWN
  the exposed path; refresh `up[]`/`edge_set[]` for the O(exposed) nodes.
- **Ascend** child c → parent p: inverse — re-remove that chunk.
- **Suppress-node** (clip parent nx removed, sibling ns joined to grandparent nz):
  the single most bug-prone step — `up[]` at ns's new position and everything it
  newly exposes must be refreshed (advisor). The amortised total over a full DFS
  is O(N) (each node exposed/hidden O(1) times).

**Gate (hard, non-negotiable):** under `-DNDEBUG`-off, assert the incrementally-
maintained `edge_set` equals the from-scratch `compute_insertion_edge_sets`
result FOR EVERY CLIP, full array, across all seeds/datasets, BEFORE trusting any
wall number. This is the correctness spine for the suppress-node transition.

**Integration:** new gated inner-loop path in `tbr_search` (env/param flag),
active only in `l3b_regime` (already wired) + DFS clip order (config B, already
wired). Replaces restore-between-clips with the Euler-tour move. Falls back to A
otherwise. Measure **C-vs-B** (bit-identical trajectory → pure wall delta) and
**C-vs-A** (net of the ~20% B-vs-A candidate cost) on Wortley/Zanol/Zhu.

---

## FINAL VERDICT (2026-06-19) — L3b is DEAD by direct measurement. STOP.

The "measure-first" reframe (advisor): the 24-33% is the *share*, not the
realizable saving. The realizable saving is bounded by the per-clip **changed-
value footprint** (Scheme 1) and the **per-descend-step delta** (Euler) — both
now measured directly on the from-scratch path (`-DTS_EDGESET_FOOTPRINT`,
`src/ts_edgeset_footprint.h`; driver `dev/profiling/drivers/l3b_footprint.R`;
rows `l3b_footprint.csv` + `l3b_euler_gate.csv`).

| dataset | tips | clips | per-clip footprint (frac of edges) | Euler descend delta / footprint |
|---|---|---|---|---|
| Wortley2006 | 37 | 11-44k | **0.68** | **1.14** |
| Zanol2014 | 74 | 22-24k | **0.46** | **1.20** |
| Zhu2013 | 75 | 19-23k | **0.41** | **1.24** |

**Both gates FAIL, decisively:**
1. **Scheme 1 (patch-from-full-tree) — DEAD.** Footprint = 41-68% of *all* edges
   change per clip (gate threshold for GO was <~0.3). Fitch combine does NOT
   saturate on these datasets. Memory-traffic floor even at the best 0.41:
   0.41×(recompute) + 0.41×(undo-save) + 0.41×(undo-restore) + boundary checks
   ≈ 0.41×3 > 1.0 — the incremental path moves MORE memory than the clean
   bandwidth-bound sweep it would replace (L2). Confirmed loss.
2. **Euler-tour (cross-clip) — DEAD.** Per-descend-step delta is **1.14-1.24× the
   per-clip footprint** — i.e. moving the clip boundary one parent→child step
   flips AS MANY OR MORE views than a from-scratch clip. The delta = tiny chunk
   (5-7 newly-exposed nodes) + LARGE common-changes (44-67 pre-existing nodes
   whose edge_set flips). **There is no cross-clip locality**: a small topological
   change near the clip point propagates view-changes across ~half the tree via
   the non-saturating intersect-else-union. This is *optimistic* for Euler (omits
   early-term boundary checks; ignores the measured ~1.2× DFS candidate cost and
   the heavy reroot/restore integration risk).

**Why the original "irreducible O(N)" proof was right after all (and the advisor's
saturation hope didn't pan out HERE):** the proof was a worst-case bound; the
advisor correctly noted saturation *could* shrink it. Measurement settles it —
on these EW morphological matrices the directional views are genuinely sensitive,
so the realizable footprint is ~half of N, not a small constant. The constant-
factor framing (T-P5h) holds.

**Disposition:** the within-clip kernel is the ceiling — **L1 (zero-fill skip,
MERGED) is the bankable per-iteration win; L2/L3a AT-LIMIT; L3b ABANDONED.** All
worktree code here is MEASUREMENT-ONLY (`-DTS_EDGESET_FOOTPRINT` /
`-DTS_EDGESET_TIMING`, gitignored Makevars) — nothing to merge to cpp-search.
Mission per-candidate throughput is at-limit for this class; focus shifts to
**recipe composition + the unexplored sectorial ~30% (#39/#40)**, where the
remaining wall actually lives. Ratchet 12→6 (T-P5d) stays the recipe win.

**Scope of the claim (advisor):** DEAD *for the mission's dataset class*
(37-75t EW morphological). `fp_frac` is data-dependent — **0.66 (Wortley) vs
0.46 (Zanol/Zhu)** — lower on the bigger/denser matrices. A much-larger-tree or
molecular-scale regime sits in UNMEASURED territory where the fraction could
cross the <0.3 GO line. **Reopen condition:** re-run `l3b_footprint.R` at large N
— the flag-gated instrumentation is KEPT in the worktree precisely so this is a
~10-min rerun, not a rebuild. Do NOT delete it.

**What this does NOT kill (advisor, so the conclusion isn't over-read):** the
measurement kills *incremental maintenance of the exact `edge_set`*. It does not
touch the **bound-then-verify / lazy-exact** route — screen candidate edges with
a cheap bound and compute the exact `edge_set` only for survivors (TNT's actual
quick-TBR; the Goloboff-93-bound vs Gladstein-exact fork from the lit sweep).
Recorded as **distinct, untested, NOT NOW** — caveated: per-candidate scoring is
already declared at-limit (T-P5c), and a cheap *admissible* (never-screens-an-
improver) lower bound on Fitch insertion cost is not established (the old
union-of-finals bound OVERcounts, the very bug the directional fix cured).

### Worktree state (claude/l3b-edgeset-timing, off 00d73d6a)
- `src/ts_edgeset_timing.h` (+ global in ts_fitch.cpp; report in ts_rcpp.cpp;
  registered in TreeSearch-init.c) — bucketed timing, `-DTS_EDGESET_TIMING`.
- `order_clips_dfs` + `l3b_regime`/`dfs_clip_order` gate + branch in ts_tbr.cpp
  (config B). Ratchet 12→6 in R/MaximizeParsimony.R default preset.
- Drivers: `dev/profiling/drivers/l3b_ceiling.R`, `l3b_ab.R`.
- `src/Makevars.win` (gitignored) carries `-DTS_EDGESET_TIMING`.
- NEXT: implement the Euler-tour kernel above + oracle micro-bench.
