# Isolated-TBR head-to-head: TreeSearch vs TNT 1.6 from identical start trees

**Date:** 2026-06-18
**Branch:** `claude/competent-chaum-6ecb56` (worktree off `cpp-search`)
**Question (project lead):** Given the *same* starting tree, how does the score
change after TBR branch-swapping? The ensemble behaviour *should* be identical in
TNT and TreeSearch. Step 1 — is there a meaningful difference? Step 2 (only if
yes) — how does TNT implement TBR, and what explains the difference?

This investigation isolates **TBR branch-swapping** from the Wagner
starting-tree confounder by feeding the **identical Newick start tree** into both
engines and running TBR to convergence. (The Wagner half is a separate task.)

## TL;DR

- **Step 1 = YES, and the difference is large.** From an *identical* poor start
  (e.g. a 1478-step Wagner tree on Zanol2014), TNT's TBR reaches ~1265 while
  TreeSearch's TBR reaches ~1300–1350. A 40–90 step gap in the hill-climb alone,
  with the starting tree held fixed.
- **Step 2 = TreeSearch's TBR neighbourhood is root-dependent.** `tbr_search`
  declares convergence at trees that are *not* unrooted-TBR local optima — TNT
  improves them, and re-rooting the same tree and re-running TS improves them too.
  Proven three independent ways (cross-feed, root-dependence test, code). The
  leading interpretation is that TS implements **rooted** TBR (a subset of TNT's
  unrooted TBR), with the fixed root blocking root-crossing rearrangements.
- **Quantified:** making TS root-invariant recovers **~half** the gap (a large,
  real effect) but leaves a **+15–36 residual** to TNT. So the root-dependent
  neighbourhood is a **proven major contributor**, not yet shown to be the whole
  cause. This still **redirects the long-standing EW-Fitch gap** away from
  "sectorial architecture" toward the TBR move set itself — a concrete, at-least-
  half-the-gap, fixable kernel deficiency.

## Method & comparability controls

- Datasets are EW-Fitch-converted (inapplicable tokens → `?`), so **both engines
  optimise the identical Fitch objective**; `TreeLength` (TS) and TNT `length`
  are directly comparable. Verified: T0 round-trips at 1271 in both.
- **TS entry point:** `TreeSearch:::ts_tbr_diagnostics(edge, ...)` — runs TBR to
  convergence from a warm-start edge matrix, returns final score + per-pass
  trajectory. `acceptEqual=FALSE, maxHits=1` = strict descent to first local
  optimum; `acceptEqual=TRUE` = single-tree plateau-walk.
- **TNT entry point:** `bbreak = tbr [no]randclip [no]mulpars;` with `tread` of
  the shared start tree and `rseed N`. `bbreak` swaps the *in-memory* tree — it
  does **not** re-randomise (verified: bbreak from T0=1271 stays 1271).
- **Two modes:**
  - **Mode A — strict single tree:** TS `acceptEqual=F`; TNT `nomulpars hold 1`.
  - **Mode B — buffer / plateau:** TS `acceptEqual=T`; TNT `mulpars hold 1000`.
    *Asymmetric by construction* (TNT swaps a buffer of equal-length trees; TS
    `tbr_search` walks a single tree — no buffer re-swap). Mode B is a mechanism
    probe, not a controlled comparison. The headline is Mode A.
- For robustness, TNT's swapped trees are saved (`tsave`) and **re-scored in R
  with `TreeLength`** — the final length is never parsed from TNT's stdout.
- Scripts: `dev/benchmarks/tbr_shared_start_lib.R` (helpers),
  `tbr_pilot.R`, `tbr_verify.R`, `tbr_crossfeed.R`, `tbr_grid.R`,
  `tbr_reroot_recovery.R`. Raw results: `dev/benchmarks/tbr_results/`.

## Pre-flight gates (all passed)

1. **Length identity:** TS tree → Newick → TNT `tread` → `length` = 1271 =
   `TreeLength`. Objective identical; round-trip faithful.
2. **Both engines seed-stochastic:** TNT `randclip`+`rseed` varies the trajectory
   (norandclip is deterministic = 1273 on the 1478 start); TS RANDOM clip order
   seeded by `set.seed`. So the "ensemble across seeds" framing is valid.
3. **TS converges genuinely:** strict descent from 1478 = 51 productive passes →
   first local optimum (not a truncated run). `bbreak` from T0 holds at 1271 in
   both engines.

## Step 1 — ensemble result (the deliverable table)

Final length over **6 seeds** per (start tree × engine), from the IDENTICAL
shared start. Six start trees per dataset spanning a quality ladder (two random
topologies, two RAS Wagner, one partially-TBR-optimised, one near-optimal
anchor). `gap` = median(TS) − median(TNT). Raw rows:
`dev/benchmarks/tbr_results/tbr_grid_raw.csv`; shared starts:
`<dataset>_starts.nwk`.

### Zanol2014 (n=74), Mode A — strict single-tree TBR

| start | start len | TNT (min/med/max) | TreeSearch (min/med/max) | gap |
|---|---|---|---|---|
| random1  | 2353 | 1264 / 1267 / 1275 | 1295 / 1318 / 1335 | **+51** |
| random2  | 2274 | 1267 / 1268 / 1272 | 1302 / 1327 / 1336 | **+59** |
| wagner1  | 1711 | 1265 / 1267 / 1270 | 1293 / 1306 / 1327 | **+40** |
| wagner2  | 1584 | 1263 / 1266 / 1274 | 1297 / 1300 / 1312 | **+34** |
| partial  | 1516 | 1262 / 1265 / 1271 | 1289 / 1296 / 1316 | **+30** |
| t0anchor | 1271 | 1271 / 1271 / 1271 | 1271 / 1271 / 1271 | +0 |

### Zanol2014 (n=74), Mode B — buffer (TNT mulpars hold 1000) / TS plateau

| start | start len | TNT (min/med/max) | TreeSearch (min/med/max) | gap |
|---|---|---|---|---|
| random1  | 2353 | 1262 / 1262 / 1271 | 1289 / 1304 / 1336 | **+42** |
| random2  | 2274 | 1262 / 1265 / 1267 | 1306 / 1328 / 1368 | **+64** |
| wagner1  | 1711 | 1261 / 1262 / 1263 | 1295 / 1304 / 1352 | **+42** |
| wagner2  | 1584 | 1262 / 1263 / 1268 | 1292 / 1308 / 1338 | **+45** |
| partial  | 1516 | 1262 / 1263 / 1266 | 1300 / 1310 / 1317 | **+46** |
| t0anchor | 1271 | 1267 / 1267 / 1267 | 1271 / 1271 / 1271 | +4 |

**Reading.** TNT lands at ~1262–1271 from *any* start (tight, low variance);
TreeSearch lands at ~1289–1368 (≈30–65 steps higher, with much wider spread).
The difference is consistent across the whole quality ladder and both modes.
Two telling cells: (i) at the near-optimal anchor both engines hold 1271 under
strict TBR — TS *can* sit at the optimum; (ii) under Mode B, TNT's buffer
*escapes* 1271 → 1267, while TS stays stuck at 1271 — TNT's neighbourhood
contains moves TS cannot see even at the optimum.

### Zhu2013 (n=75) — the gap is larger still

Mode A (strict). Same pattern, bigger magnitude (TNT target ≈ 624):

| start | start len | TNT (min/med/max) | TreeSearch (min/med/max) | gap |
|---|---|---|---|---|
| random1  | 1833 | 628 / 630 / 634 | 648 / 682 / 779 | **+52** |
| random2  | 1813 | 626 / 630 / 638 | 686 / 716 / 732 | **+86** |
| wagner1  | 1261 | 627 / 632 / 633 | 683 / 691 / 777 | **+60** |
| wagner2  | 1195 | 626 / 632 / 636 | 693 / 734 / 790 | **+102** |
| partial  | 1342 | 625 / 627 / 632 | 670 / 706 / 768 | **+79** |
| t0anchor | 631  | 631 / 631 / 631 | 631 / 631 / 631 | +0 |

Mode B (buffer). TNT reaches the project **target ≈ 624** from random starts:

| start | start len | TNT (min/med/max) | TreeSearch (min/med/max) | gap |
|---|---|---|---|---|
| random1  | 1833 | 624 / 625 / 627 | 665 / 691 / 731 | **+66** |
| random2  | 1813 | 624 / 624 / 627 | 690 / 728 / 765 | **+103** |
| wagner1  | 1261 | 624 / 626 / 627 | 676 / 688 / 718 | **+62** |
| wagner2  | 1195 | 624 / 626 / 634 | 672 / 729 / 823 | **+103** |
| partial  | 1342 | 625 / 625 / 627 | 665 / 692 / 743 | **+67** |
| t0anchor | 631  | 625 / 625 / 625 | 631 / 631 / 631 | +6 |

The effect is robust across both datasets and **larger** on Zhu (+50 to +100).
Again the buffer escapes the anchor (631 → 625, toward target 624) while TS is
stuck. Step 1 is an unambiguous YES on both datasets.

## Step 2 — mechanism: TreeSearch's TBR is rooted

The deficit is **not** "TS reaches a worse basin"; it is **TS terminates before a
true (unrooted) TBR local optimum** because its move set is root-restricted.

**(a) Reciprocal cross-feed (decisive).** Feed each engine's converged optimum
into the other:

| Fed tree | Into | Result |
|---|---|---|
| TS local optimum **1302** | TNT `bbreak` nomulpars (deterministic) | → **1267** |
| TS local optimum 1302 | TNT mulpars hold 1000 | → 1262 |
| TNT local optimum **1266** | TS strict TBR | → **1266** (holds; converged) |

TNT finds strictly-improving moves from a tree TS declared a local optimum, while
TS *holds* TNT's optimum (no wander-above ⇒ no scoring/round-trip artefact; TS
*can* represent it, it just can't path there). ⇒ neighbourhood incompleteness.

**(b) Root-dependence (engine-internal proof).** Fitch length is root-invariant,
so every re-rooting of the TS 1302 optimum is still length 1302. Re-running TS
strict TBR from those re-rootings:

| reroot at | Aciculomarphysa | Eunice_fucata | Leodice_americana | Leodice_thomasiana | Mooreonuphis | Palola_B5 |
|---|---|---|---|---|---|---|
| TS TBR final | 1296 | 1295 | **1281** | **1281** | 1291 | 1286 |

An *unrooted* TBR local optimum would hold at 1302 for every rooting. It does
not (down to 1281) ⇒ **the TS TBR neighbourhood depends on the root.**

**(c) Code (`src/ts_tbr.cpp`).** The kernel uses a rooted tree representation:
clips whose parent is the root are skipped (L804); only the *smaller* subtree of
each edge is clipped (L812); TBR rerooting is applied only to the **clipped
subtree**, never the main tree. Rearrangements that cross the fixed root are
therefore unreachable — the textbook definition of rooted (vs unrooted) TBR.

**Confirmation — emulated root-invariance recovers ~half the gap, but not all.**
Wrapping the shipping rooted kernel in an outer reroot-sweep loop (TBR → try
re-rootings → TBR, looped to convergence) over **all 74 tips** (Zanol):

| start | seed | TS rooted | reroot-invariant (strict) | reroot-inv + plateau | TNT |
|---|---|---|---|---|---|
| wagner (1711) | 1 | 1304 | 1292 | 1289 | **1265** |
| wagner | 2 | 1326 | 1304 | – | **1268** |
| random (2353) | 1 | 1330 | 1284 | 1279 | **1264** |
| random | 2 | 1295 | 1288 | – | **1264** |

Root-invariance recovers roughly **half** the strict gap (e.g. random seed 1:
1330 → 1284, recovering 46 of the 66 steps to TNT) — a large, real effect that
**proves the root-dependent neighbourhood is a first-order cause.** But a
**+15 to +36 residual to TNT remains**, and it is *not* closed by also allowing
plateau-crossing (reroot-inv + plateau ≈ strict). So root-dependence is a
**proven major contributor, not the whole story.** The residual is either (i)
*incomplete* root-crossing — this emulation reroots only the *converged* tree
between full TBR runs, whereas true unrooted TBR crosses the root *within* every
sweep, so a proper integrated implementation should beat this emulation — or
(ii) a genuine second neighbourhood/acceptance difference (candidate suspects in
`ts_tbr.cpp`: the smaller-subtree-only clipping L812, and the collapsed-edge
pruning L817/L919). Disentangling (i) from (ii) is the natural follow-up.

**Plateau-crossing is not the gap.** TS single-tree plateau-walking was tested
directly (`acceptEqual=TRUE`, `maxHits` ∈ {1, 5, 50, 500}) and does not help —
TS still lands ~1290–1350, no better than strict descent. So Mode A does not
unfairly deny TS the equal-length moves TNT's `nomulpars` takes; the deficit
survives giving TS those moves.

## Recommended fix

Make TBR root-invariant in `ts_tbr.cpp`: evaluate main-tree re-rootings *within*
the neighbourhood (true unrooted TBR — preferred, since the between-pass
emulation already recovers ~half and within-pass should do better), or as a
cheaper first cut an outer reroot-per-round loop (cf. the reroot-per-round fix
already used in tree-fusing). Expect to recover at least half the gap; then
investigate the residual (smaller-subtree clipping L812, collapsed-edge pruning
L817/L919) to close the rest. Re-run this shared-start harness to measure.

**Output caveat (per project lead):** rerooting *during* search is free —
Fitch length is root-invariant — but the search root is an internal device only.
When the final tree(s) are returned to the user they **must be re-rooted onto the
originally-specified outgroup**, so the displayed topology matches the user's
rooting. The internal reroot must not leak into the user-facing result.

## Apples-to-apples caveat resolved

`help mult` confirms plain `mult`/`mult=replic 1` is **one RAS + TBR**;
`ratchet`/`drift`/`fuse` are opt-in flags, off by default. So the prior
"1 TNT `mult` rep → 1271 vs 20 TS reps → 1287" comparison was *not* TNT secretly
running ratchet/sectorial — its only confounder was TNT's own RAS Wagner start.
This shared-start test removes even that, and the gap remains: it is the TBR
move set, not the starting tree.

---

## ADDENDUM (2026-06-18, same day): root cause is kernel move-INCOMPLETENESS, not just rootedness

The "root-dependence recovers ~half, residual unexplained" reading above is
**superseded**. A gating cross-feed + a kernel-independent neighbourhood probe
pinned the residual precisely.

### The gating cross-feed (`tbr_reroot_crossfeed.R`)

Feeding the **all-tips reroot-invariant** TS optimum (≈1284 — which *should* be a
complete unrooted-TBR optimum) into TNT `bbreak`:

| start | TS reroot-invariant opt | → TNT `nomulpars` | TNT own opt → TS reroot-inv |
|-------|------------------------:|------------------:|----------------------------:|
| wagner s1 | 1292 | **1262** | 1265 → 1265 (holds) |
| wagner s2 | 1304 | **1268** | 1268 → 1268 (holds) |
| random s1 | 1284 | **1270** | 1264 → 1264 (holds) |
| random s2 | 1288 | **1268** | 1264 → 1264 (holds) |

TNT (even single-tree `nomulpars`) **improves** the TS optimum; TS **holds** TNT's.
Asymmetric ⇒ TNT's TBR neighbourhood strictly contains moves the TS all-tips
search lacks. The residual is **neighbourhood, not basin/path**.

### The kernel-independent probe (`tbr_neighbourhood_probe.R`)

`TBRMoves`/`SPRMoves` (→ `all_tbr`/`all_spr` in `rearrange.cpp`, a separate
UNOPTIMISED enumerator — no L812, no collapsed) on the TS 1284 optimum:

- **43 improving TBR neighbours** (best 1280); **26 improving SPR neighbours**
  (best 1280). So the deficit is at the **basic clip+graft (SPR) level**.
- TNT's 1264 optimum: **0 improving** — TNT reaches genuine canonical-TBR optima;
  it does **not** exceed textbook TBR.

So the TS kernel **falsely declares convergence** while real improving moves exist.

### Mechanism: a STACK of completeness-breaking optimisations in `ts_tbr.cpp`

Validated with gated fixes behind a new opt-in `TBRParams::unrooted` (default off;
the DEFAULT search is byte-identical — confirmed 1330/1295 unchanged). Each fix
peels missed moves and lowers the all-tips optimum (Zanol2014, random start):

| kernel state | optimum | enumerator-improving |
|--------------|--------:|---------------------:|
| baseline (shipping) | 1284 | 43 |
| + collapsed pruning OFF | 1284 | 43 — **collapsed is NOT a cause; ruled OUT, keep it** |
| + L812 smaller-subtree skip relaxed | 1274 | 7–15 |
| + nz/ns graft skip fixed | 1270 | 4–6 |

(≥1 further pruning remains → 4–6; `sp==clip_node` skip and the vp-dedup were
checked and are **sound**, so the residual is something else, not yet pinned.)
Trend: TS converging toward TNT's enumerator-clean 1264.

The two confirmed bugs:

- **L812** (`clip_size > n_tip/2` skip): clipping only the smaller side is meant to
  reach larger-side moves via fragment-reroot + graft-at-original — but that graft
  is killed by the next bug, so **SPR-prune-larger-subtree** moves are lost.
- **nz/ns skip** (`above==nz && below==ns` in the *rerooting* loops): correct for
  the non-rerooted SPR loop (it's the identity) but **unsound** in the rerooting
  loops — a rerooted fragment regrafted at its original location is a *distinct*
  valid move. Fixing nz/ns alone (keeping L812) likely restores L812's soundness —
  the preferred production fix (keeps the perf optimisation).

This is a **correctness/soundness bug in the package's default TBR (and SPR)** —
every TBR search lands at non-canonical-TBR-optimal trees — materially bigger than
the EW-Fitch benchmark itself. The root-dependence finding (L804) is the *separate*
≈half, handled by the reroot mechanism. (Scripts: `tbr_reroot_crossfeed.R`,
`tbr_neighbourhood_probe.R`, `tbr_collapsed_test.R`. Logs in `tbr_results/`.)

---

## FINAL (2026-06-18, post-merge of the directional-vroot fix)

The "stacked L812/nz/ns enumeration bugs / soundness bug" framing above was
diagnosed on the **pre-fix** build and is **largely superseded**. After merging
`cpp-search` commit `2b299e4b` (the parent's EW-directional **scoring** fix —
"vroot" is a candidate-cost fix, not a root mechanism), the differential oracle
(`tbr_oracle.R`: run the in-kernel `tbr_search` to convergence, assert
`all_tbr`/`all_spr` 0-improving) shows:

| kernel state | oracle failures (random 12-tip) |
|---|---|
| pre-fix default | 23/40 |
| **post-fix default (directional scoring fix)** | **9/60** |
| post-fix + all-tips rerooting | **0/60** |

So the *scoring* under-count (union-of-finals → wrong abandonment cutoffs hiding
improving candidates) was the bulk of the apparent incompleteness; **the L812/nz/ns
move-edits were pre-fix artifacts and are not needed** (kept stashed, not applied).
The whole residual is the rooted-representation **root-edge limitation** (cannot
break the root edge; with the smaller-side clip filter also cannot clip edges whose
smaller side holds the root) — covered by re-rooting.

### What was built (opt-in, default off)

`TBRParams::unrooted` + an in-kernel **reroot-at-convergence** loop in
`tbr_search` (`ts_tbr.cpp`): after converging at one rooting, re-root at the next
tip and re-descend; stop when a full tip-sweep yields no strict improvement.
Score is root-invariant, so a re-root only changes the representation. Gated to the
plain search (no sector/constraint/tabu/pool). Exposed via `ts_tbr_diagnostics(...,
unrooted=)`. The **default path is unchanged** (`unrooted=FALSE`).

Validated: oracle in single-call mode (kernel re-roots internally) → **0/60 at
12 tips, 0/40 at 16 tips**; one real-data 74-tip Zanol check → canonical-TBR-clean.

### Cost / benefit (Zanol2014, `tbr_unrooted_validate.R`)

| start | rooted len | unrooted len | gain | time × |
|---|--:|--:|--:|--:|
| random (poor) | 1272 | 1265–1271 | 0–7 | ~3× |
| RAS-Wagner (good) | 1267–1280 | 1267–1279 | 0–1 | ~10× |

Median: rooted 1272 → unrooted 1269 (**gain ≈ 3**, **0–1 from Wagner starts**);
**median ≈ 6.5× wall-clock per `tbr_search` call**. It reaches *true* unrooted-TBR
optima but does **not** close the gap to TNT (1262–1264): 1265–1279 are clean
single-tree optima — the residual to TNT is **basin/escape (multi-tree/buffer)**, a
separate mechanism, not neighbourhood completeness.

**Recommendation:** the directional scoring fix (already merged) is the real win.
The reroot mechanism is correct but its production value is marginal — production
uses Wagner starts (gain ≈ 0–1) and it costs ~6.5×. Keep it **opt-in** (or for a
final high-effort pass), not the default, unless a cheaper variant (relax-L812 +
direct in-kernel root-edge break, single pass) is built and shown to hold quality.
(Scripts: `tbr_oracle.R`, `tbr_unrooted_validate.R`, `tbr_collapsed_test.R`.)
