# Beam sectorial: pool-aware RSS over a diverse, suboptimal-tolerant buffer

Date 2026-06-18. Worktree `C:/Users/pjjg18/GitHub/TS-selectem`, branch
`claude/selectem-diversity` (off `cpp-search`). NOT on cpp-search. Env-gated;
default path byte-identical when the flag is unset.

## The convergent diagnosis (why single-tree sectorial plateaus)
Three independent results all land on the same mechanism:

1. **Chip (TNT side):** TNT's `sectsch=rss` escape = sectorial run over a RETAINED
   diverse SET of equal-optimal trees (shared `hold` buffer). Single-tree strict
   sectorial plateaus at ~1267 forever; effort/budget cannot substitute.
2. **Budget-matched TS run (this branch):** our own single-tree sectorial, given
   TNT-like budget (20 picks x 30 rounds), reaches only ~1265-1267 (coll30 -4 on
   Zanol). Effort is not the lever.
3. **Diverse-starts test (`test_diverse_starts.R`, run bve79389o):** single-tree
   sectorial from EACH of TNT's 10 diverse 1271 trees, 30 independent lanes:
   **best 1266, 0/30 reach 1261.** Diverse *starts* alone do NOT escape.

=> The lever is neither effort, nor sector geometry, nor diverse starting points.
It is a **shared, evolving buffer**: improvements found on one tree must become
visible as starting points for later picks, AND the buffer must retain
topologically-diverse trees (including SUBOPTIMAL ones) so a sector re-solve can
reach an arrangement no single frozen tree exposes.

## The buffer-width subtlety (resolves the hold-1000-vs-10 question)
`hold 1000` is the buffer CAPACITY (explicitly set, not TNT's default); the "10
trees" is the count `mult=replic 1` deposits. The decisive point: with cap 1000
and only ~10-50 trees ever present, TNT NEVER purges -> the buffer accumulates
ALL distinct trees found across the whole length range (chip: "lengths 1261-1271
coexist at the plateau"). That suboptimal diversity fuels cross-topology sector
recombination.

Our `TreePool` defaults to `suboptimal = 0.0`: `evict()` purges everything worse
than best the instant best improves. Used as-is for a beam it collapses to
best-equal-only and discards exactly the suboptimal diversity that drives escape.
**A faithful beam needs a WIDE buffer: large `suboptimal` (or hold a length band)
+ large `max_size`, retaining diverse trees over a range of lengths.**

## Architecture (current)
- `driven_search(TreePool& pool, ...)` runs `max_replicates` reps; each rep builds
  a fresh start, runs `run_single_replicate` (single-tree sectorial/ratchet/...),
  then `pool.add_collapsed(result)`.
- `rss_search(TreeState& tree, ...)` mutates ONE tree; pool only supplies
  `split_freq` weighting. Sectorial never reads/writes the pool mid-search.
- `TreePool`: best-equal retention (`suboptimal=0`), collapsed-topology dedup,
  diversity-aware eviction when full, `best()`, `all()`, `add_collapsed()`.

## Proposed design: `beam_sectorial(TreePool& beam, DataSet&, params, cd)`
A new pool-aware sectorial driver, env-gated (`TS_BEAM`), called from the RSS
phase when enabled. Loop:
```
seed beam with the working tree (+ optional K diversification walks)
for round in 1..rss_rounds:
    T = pick_from_beam(beam)          # weighted toward better score; random among ties
    Tcopy = T
    rss_search(Tcopy, ds, sp, cd)     # one sector pass, ras_starts re-solves
    beam.add_collapsed(Tcopy, score)  # WIDE buffer: keeps diverse + suboptimal
return beam.best()
```
Beam buffer = a `TreePool` constructed with large `suboptimal` (e.g. +N steps, or
a tuned band) and large `max_size`, so it behaves like TNT `hold 1000`.

### CORRECTED design (advisor gate, 2026-06-18)
The advisor caught a real overreach: the "wide buffer" (Claim B) rests ONLY on the
chip's secondhand "1261-1271 coexist" — the same chip caught conflating capacity
with count — and the probe that would confirm it hung and never ran. The diverse-
starts test (0/30) discriminates write-back vs no-write-back; it does NOT
discriminate best-equal (Claim A) from wide (Claim B), since both have write-back.
"1261-1271 coexist" most likely = lazy-eviction ballast (hold 1000 never hit, so
intermediate trees linger), NOT evidence they are picked from.

**Ship Claim A first:** best-equal beam, `TreePool` at its default `suboptimal=0`
(its diversity-aware eviction already keeps the spread; do NOT write a new buffer
class). Minimal faithful beam:
1. Seed the beam from the working tree.
2. Each round: pick from the best-equal set (uniform among ties), copy, ONE sector
   re-solve, `add_collapsed` back.
3. **`accept_equal` ON in the sector re-solve.** This is the diversity engine:
   from a single T0 seed with accept_equal=false, the re-solve returns the tree
   unchanged when no strict improvement -> `Tcopy == T` -> add_collapsed sees a
   duplicate -> diversity never grows -> beam degenerates to single-tree. Accepting
   equal-length rearrangements and writing the DISTINCT ones back is what makes a
   beam exist at all. This also explains why accept_equal HURT before (1271): on a
   single tree it is a directionless random walk; inside a retained buffer it is
   the diversity engine. **Beam + accept_equal together is the actual test.**

**`suboptimal` is a default-0 KNOB** (already a `TreePool` ctor arg). suboptimal=0
(A) vs large (B) becomes a one-line experiment AFTER the beam exists -> **drop
probe_hold entirely**; the knob answers buffer-width in-system, faithfully.

**Budget-match:** total sector searches = rounds x picks_per_rss = ~600 (== the
coll30_20 single-tree baseline: 20 picks x 30 rounds). Count it explicitly or any
win is the budget confound, not the architecture.

**Integration:** new self-contained `beam_sectorial` with a LOCAL `TreePool`,
invoked in `run_single_replicate`'s RSS phase under `TS_BEAM`. Keeps default path
byte-identical; works under the harness's `maxReplicates=1, tree=T0` (self-seeds
from T0). Do NOT graft any fuse step into the beam loop (74-78 tips = the
[[fuse-reroot-segfault]] >64-tip zone; beam path doesn't fuse, so we are clear).

**Decision rule:** best-equal beam reaches 1261 -> done, skip the riskier buffer.
Stalls -> THEN widen `suboptimal`, and we will KNOW width is the lever, not guess.

## Target (define_target.R, ratchet-off TNT mult+sectsch, canonical hold-1000 T0)
| dataset | n | T0 | TNT target | gap |
|---|---|---|---|---|
| Zanol2014 | 74 | 1271 | 1261 | -10 |
| Wortley2006 | 37 | 485 | 480 | -5 |
| Zhu2013 | 75 | 631 | 624 | -7 |
| Giles2015 | 78 | 672 | 670 | -2 |
Single-tree TS baseline from same T0: ~0 (coll30 -4 Zanol/Zhu only). Beam must
beat this to justify the architecture change.

## RESULTS (bench_beam.R, canonical T0, budget-matched 30 rounds x 20 picks, seeds 1-3)
| configuration | Zanol (tgt 1261) | Zhu (tgt 624) |
|---|---|---|
| single-tree (baseline) | 1267 (-4) | 627 (-4) |
| beam, best-equal (Claim A) | 1266 (-5) | 627 (-4) |
| beam, wide subopt=10 + pick-all (Claim B) | 1266 (-5) | 627 (-4) |
| independent lanes, 10 diverse seeds, NO sharing (prior) | 1266 | - |
| **TNT: shared buffer over 10 diverse seeds** | **1261** | **624** |

**Both Claim A and Claim B plateau at 1266/627 == the single-tree/diverse-starts
floor.** Buffer WIDTH is moot here, and the reason is structural:

### Why single-seed beam is seed-starved (the missing ingredient)
`rss_search` NEVER returns a tree worse than its input (it reverts any sector move
that worsens the full score, and its final global TBR only improves). So in a
beam seeded from ONE T0, the only suboptimal trajectory is the T0 seed itself —
there is no source of genuinely diverse suboptimal trees for a wide buffer to hold
or pick from. Widening (Claim B) therefore changes nothing.

TNT's `mult` supplies ~10 DIVERSE 1271 seeds -> 10 independent descent
trajectories pooled in the shared buffer. The diverse-starts test already showed
those 10 seeds WITHOUT sharing = 1266. The one untested cell is **diverse seeds
WITH shared-buffer write-back** — exactly TNT's recipe. That, not buffer width, is
the next lever.

## NEXT: multi-seed beam (advisor fork)
Beam must seed from K diverse trees, not one. Implementation options under
consideration (advisor gate before more C++):
- (a) beam generates K seeds internally via RAS+TBR (faithful TNT `mult`); env
  TS_BEAM_SEEDS=K. Self-contained; abandons the fixed-T0 comparison (TNT's sectsch
  also doesn't start from a fixed T0).
- (b) seed beam from the SAME 10 TNT diverse trees (via plumbing a multiPhylo /
  file) — cleanest apples-to-apples vs the diverse-starts 1266, isolates "sharing"
  as the sole added variable.
- Budget accounting changes with K seeds (K extra TBR searches); primary question
  first ("does multi-seed beam reach 1261 at all"), wall-clock fairness second.

## MULTI-SEED RESULT — beam architecture RULED OUT as the gap-closer
First multi-seed attempt had a BUG: seeds were generated by random-addition Wagner
+ one TBR pass, which lands ~20-90 steps WORSE than T0 (1291-1358 on Zanol) —
outside the basin — so `add_collapsed` rejected all of them (subopt=10 threshold
1281). The "beamMulti = 1266" was a silent single-seed run. Fixed: seed by
plateau-collecting TBR from T0 (`accept_equal` + `collect_pool`), which gathers
distinct T0-basin (1271) trees directly into the beam.

After the fix, the beam genuinely seeds 3-4 distinct 1271 trees. Result on Zanol
(seed 2): **1267** — NO better than single-seed. Full picture:
| configuration | Zanol | Zhu |
|---|---|---|
| single-tree | 1267 | 627 |
| beam best-equal, single-seed | 1266 | 627 |
| beam wide, single-seed | 1266 | 627 |
| beam multi-seed (real diverse) + wide + pick-all | 1266-1267 | 627 |
| **TNT target** | **1261** | **624** |

**Every faithful beam variant plateaus at ~1266/627 — the single-tree + diverse-
starts floor.** What this PROVES: the shared buffer is NOT SUFFICIENT — the chip's
thesis (dev/plans/2026-06-18-tnt-sectsch-superpower.md) that the buffer is the
gap-closer is REFUTED. What this does NOT prove: that the buffer is useless. Every
beam variant calls `rss_search`, which uses the FROZEN-HTU sector re-solve — the
same one single-tree uses. So "beam plateaus where single-tree does" is exactly
what a frozen-HTU bottleneck would produce, whether or not the buffer is useful.
The experiment cannot discriminate "beam useless" from "beam capped by the
re-solve." Defensible claim: **beam-on-frozen-HTU = the floor; the re-solve is the
binding constraint.** Keep the beam behind its flag for a beam+HTU re-test once
the re-solve can float (below).

## REDIRECT: sector re-solve QUALITY (HTU floating, task #24 / D1)
The ~1266 floor sits exactly +5 above target. My own prior audit
([[sector-resolve-status]], dev/plans/2026-06-17-tnt-algorithm-audit.md, task #24)
pinned a CONFIRMED per-sector quality gap: TNT FLOATS the HTU pseudo-tip during
the sector re-solve (joint re-resolve x re-attach = a barrier crossing); we FREEZE
it. A frozen-HTU re-solve structurally cannot produce the arrangements TNT's can,
so NO buffer/beam machinery closes the gap — the moves aren't reachable. There is
already a scoring-only probe (`TS_FREE_HTU_PROBE`, ts_sector.cpp ~L978) confirming
a free-HTU re-solve finds lower reduced scores. The lever is implementing the
free-HTU re-solve + reattach, not buffer architecture.

## GATE before building free-HTU reattach (advisor)
Do NOT start the hard D1 reattach on "+5 floor ~ HTU." HTU is "+1/+3 PER SECTOR";
whether that accumulates to the +5 plateau gap AT THE PLATEAU is unverified. Cheap
check first: run a 1266-RESIDENT tree through rss_search with TS_FREE_HTU_PROBE on,
count `<<D1-CONFIRM` fires (free-HTU reduced score < anchored) AT THE PLATEAU
(not generically on a 1271 tree).
- Fires often at 1266 -> lever is live where we're stuck; build the reattach.
- Rarely fires at 1266 -> reduced-score headroom gone at plateau; float-HTU won't
  help either; floor is something else -> saved days.
Caveat: lower reduced score is necessary not sufficient — the reattach must
REALIZE it on the full tree. Positive probe => "build + verify full-tree drop."

## After HTU floats: re-test beam+HTU (do not assume)
Once the re-solve can float, run single+HTU vs beam+HTU. THAT discriminates
whether the buffer was ever load-bearing. single+HTU hits 1261 -> beam was a dead
end, delete cleanly. single+HTU stalls but beam+HTU breaks through -> buffer was
necessary, would have wrongly killed it.

## Separate flag (log, do NOT fold into HTU work)
TS's RAS+TBR-from-random lands 1291-1358 on Zanol = 20-90 steps above T0 (1271).
TS's per-replicate TBR descent is materially WEAKER than TNT's `mult`. Controlled-
for here (fixed T0) but bears on PRODUCTION where TS builds its own starts. Needs
its own investigation.

## HTU-FLOAT GATE RESULT — float-HTU also RULED OUT for the plateau (Zanol)
Ran TS_FREE_HTU_PROBE on T0 (1271) vs a plateau tree (1267), 50 sector probes each,
[31,99] coll30 sectors (S=30), 20 free RAS+TBR restarts/sector.
- **T0 (1271):** exactly ONE sector (110) fires `<<D1-CONFIRM` (free 529 < anchored
  533, -4). All others free >= anchored (cold-search weakness — free RAS+TBR from
  random underperforms the anchored search seeded from the existing good subtree).
- **Plateau (1267):** sector 110 now anchored = free = 529 — the -4 headroom is
  GONE (the anchored sectorial captured it during the 1271->1267 descent). **ZERO
  D1-CONFIRM at the plateau.**
Gate verdict (advisor): detectable free-HTU headroom is EXHAUSTED at the plateau ->
**floating the HTU will NOT break the 1267 floor.** The probe's free search is weak,
but equally so at T0 and plateau, and it DID detect the T0 sector-110 headroom, so
the relative signal (fires at T0, silent at plateau) is trustworthy. This saves the
hard D1 reattach build. NB the necessary-not-sufficient caveat cuts the other way
too: a weak cold search could under-detect, but the T0-vs-plateau contrast holds.

So BOTH the beam (buffer architecture) AND HTU-floating are ruled out as the
1267->1261 lever for Zanol. The floor is something else.

## NEW LEAD: core TBR/Wagner hill-climbing quality deficit
Surfaced incidentally: a single TS random-addition Wagner + one TBR pass lands at
**1291-1358 on Zanol = 20-90 steps above T0 (1271)**, whereas TNT's `mult=replic 1`
reaches 1271 (and deposits ~10 trees there). That is a 1.5-7% deficit in BASIC TBR
hill-climbing. If TS's core TBR descent is materially weaker than TNT's, EVERY
component (starts, sector re-solve, polish) inherits it, and no sectorial
architecture change closes the gap. This is the strongest remaining lead and is a
PRODUCTION concern (TS builds its own starts). Verify it's a real deficit, not a
measurement artifact (poor random_wagner start + too-few TBR restarts vs TNT's RAS).

VERIFIED (2026-06-18, /tmp/tbr_check.R): TS pure Wagner+TBR multistart, ratchet/
sectorial OFF, on Zanol:
| effort | TS best |
|---|---|
| 1 replicate | 1315 |
| 5 replicates | 1306 |
| 20 replicates | **1287** |
| TNT mult=replic 1 (ONE replicate) | **1271** |
**20 TS replicates (1287) can't match 1 TNT replicate (1271)** — +16 above T0,
+26 above target. NOT single-pass weakness recovered by more starts; a fundamental
core hill-climbing deficit. In production TS is doubly disadvantaged (worse starts
AND worse plateau-escape). The canonical shared-start comparison (both from 1271)
factored the start deficit OUT — which is why the EW work focused on sectorial; but
the start deficit is real and large on its own. CAVEAT: TBR-only here; confirm TNT
`mult` does not swap beyond a single TBR pass (apples-to-apples) — though 20-vs-1
magnitude makes a pure artifact unlikely. Next: profile/compare TS tbr_search
thoroughness (clip order, max_hits, convergence) vs TNT branch-swapping.

## Status / session conclusion
Beam built + wired behind TS_BEAM (knobs TS_BEAM_SUBOPT/PICKALL/MAXSIZE/SEEDS/
DEBUG); default path byte-identical. Findings, in order of confidence:
1. Beam (shared buffer over diverse seeds) does NOT close the Zanol/Zhu gap ->
   chip's "buffer is the lever" thesis REFUTED (not sufficient). Keep behind flag.
2. Float-HTU GATED OUT for the plateau (probe headroom exhausted at 1267).
3. NEW LEAD: core TBR/Wagner quality deficit (1291-1358 vs 1271) — investigate next.
All in worktree TS-selectem; nothing on cpp-search. Reported pivot to user.
