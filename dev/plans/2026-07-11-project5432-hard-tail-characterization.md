# project5432: the last hard-tail gap is COLD-START BASIN CAPTURE (rewritten 2026-07-14)

Capstone for the TNT-gap thread (`2026-06-16-closing-the-tnt-gap.md` →
`2026-06-18-gap-framing.md` → `-tnt-algorithm-audit.md`). Quality parity is
**CLOSED** on ≤88-tip matrices. This note characterizes the single matrix where a
gap survives — **project5432** (482 taxa × 189 chars, Fitch equal-weights,
gaps=missing, inapplicable=missing) — and records the corrected diagnosis.

> **Rewrite note (2026-07-14).** The 2026-07-11 version of this file concluded the
> gap was **clade generation** — "the 1943 tree's deepest 238-tip backbone
> bipartition exists in no TS tree … TS never explores deep enough to *generate*
> the deep-backbone split." **That conclusion is FALSIFIED.** The `restore238`
> experiment (job 17866684, 2026-07-12) broke exactly the 1943 tree's deepest
> 238-tip split and let TS's real main-loop TBR descend: **22/24 truly-adjacent
> starts recover the 238-split and rebuild the whole 1943 tree with all 64
> skeleton splits**, and a capture-radius sweep recovers 1943 from any start within
> ~8–12 TBR moves. TS's operators **can** generate the deep backbone. The wall is
> not generation incapacity — it is **cold-start basin capture**: cold RAS+Wagner
> lands ~109 TBR moves from 1943 and never routes within the capture radius. The
> reconciliation below keeps both facts.

> **CORRECTION (2026-07-16, user).** The best-known score for project5432 is **1939**
> (found by TNT), **not 1943**. Throughout this note "1943" is a strong, TS-HOLDABLE
> tree (TS's TBR holds it when seeded and rebuilds it from within ~8–12 TBR) but it is
> an **intermediate/suboptimal basin, not the optimum**. Corrected ladder: **1939**
> (TNT best-known) < **1943** (TS-holdable) < 1944 (TNT reliable) < 1945 (TS kick
> floor) < ~1947 (TS production). So the true cold-reach gap is TS 1945 → 1939 ≈ **6
> steps**, not the "~2 steps to 1943" stated below. The restore238 / basin-radius work
> characterised the **1943** basin; the **1939** basin is not yet characterised. The
> 1939 tree-file provenance is being confirmed on Hamilton (this session); see memory
> `project5432-best-known-1939`. Numbers below are left as originally written and should
> be read through this correction.

## TL;DR — cold-start basin capture, not a missing operator, not clade-generation

- **The gap is small and on one matrix.** TS's production multistart floor is
  **~1947**; the reweighting-kick arm reaches **1945** (see below); TNT's
  **reliable** result is **1944**; TNT reaches the **1943** optimum only ~1/3 of
  runs (its own runs are 4/12 at 1943, 8/12 at ≤1944, matched-wall). Gap
  TS→TNT-reliable ≈ **1–3 steps / ≤0.15%**.
- **1943 is a VALID, STABLE TS optimum.** Seeded at 1943, TS's own TBR HOLDS it
  (0 improving moves; `seed1943` job 17865612); MaximizeParsimony retains it. NOT
  a completeness / representation / rooting defect, NOT a retention/stopping defect.
- **TS's TBR REBUILDS the whole 1943 tree from any start within ~8–12 TBR moves of
  it** (`restore238` job 17866684; basin-radius `basin_5432` job 17863426:
  return-to-1943 40–60% at d≤12, 0/5 by d≈20). So an operator need **not** hit 1943
  exactly — landing within ~8–12 TBR suffices; normal TS descent finishes from there.
- **The wall is purely cold-start basin capture.** Cold RAS+Wagner lands ~109 TBR
  moves from 1943 (TBRDist `geom_full` job 17860581: 1943↔TS median ~109, same
  scale as TS↔TS internal ~106 — 1943 is *intermingled*, not a distant island) and
  never routes within the ~8–12 capture radius; ~2150+ cold reps → 0 hits. The
  mechanism gap is **confirmed by observation** (head-to-head jobs 17863597/8:
  TS 0/20 vs TNT 4/12 → 1943).
- **The mission (objective (i)): get a cold start — or a search trajectory — to
  land within ~8–12 TBR of the best basin on large trees.** This note is the entry
  point for that search space.

## The reconciliation — how "85/479 clades absent" and "TBR rebuilds 1943" coexist

Both facts are true and they are not in tension once the claim is stated precisely:

- **Operator capability (restore238):** started within capture radius, TS's TBR
  regenerates the 238-tip backbone split and every one of 1943's 64 skeleton splits
  and descends to 1943. TS **can** produce the deep backbone.
- **Trajectory reach (clade-check, jobs 17866068/17866225):** 85/479 of 1943's
  clades — including the deepest 238-tip split — are absent from *every* TS tree in
  both a 10 008-tree diverse pool (1953–1967) and TS's genuine best 34 (1947–1956).

The 2026-07-11 reading turned the second fact into "TS can never *generate* the
split." restore238 shows the correct reading: **the cold search TRAJECTORY never
ARRIVES within ~8–12 TBR of 1943**, so those clades never appear in trajectories
that all settle ~109 moves away. It is a navigation/basin-capture failure, not an
operator or representation limit. Every earlier "TS never generates the split /
fuse-unreachable / clade-generation gap" statement in the campaign should be read
this way.

## The honest scale

| quantity | value |
|----------|-------|
| TS production multistart floor (pre-kick, ras3+tbrMaxHits20, matched wall) | **~1947** |
| TS reweighting-kick arm floor (ratchetPerturbMaxMoves ≥40, saturated) | **1945** |
| TNT reliable score | **1944** |
| TNT optimum (lucky, ~1/3 of runs) | **1943** |
| gap TS-kick → TNT-reliable | **1 step** |
| gap TS-kick → TNT-optimum | **2 steps ≈ 0.1%** |
| TNT self-reliability at 1943 (matched-wall) | **4/12** (8/12 at ≤1944) |
| TS reach at 1943 from RAS starts (matched-wall) | **0/20** |
| capture radius of the 1943 basin | **~8–12 TBR moves** |
| cold-start distance to 1943 (TBRDist) | **~109 TBR moves** (≈ TS↔TS internal) |

## What is REFUTED — do NOT re-derive (each has a Hamilton job)

These stand, now correctly read as "the cold trajectory never arrives," not
"operator/generation incapacity":

| candidate lever | verdict | evidence |
|-----------------|---------|----------|
| **Distant-basin / diversity generation** | DEAD | 1943 intermingled (CID 0.33 ≈ TS-internal 0.37; TBRDist ~109 ≈ TS↔TS ~106); rogue-prune inert |
| **BSS** (freeze-clade block re-opt) | **structurally refuted** | ceiling probe K_ceiling=323/482, 240 singletons, largest 1943-compatible frozen block 8 tips → 1945's block internals are wrong *pervasively*; reaching 1943 needs near-per-tip re-opt = the full search (`ceiling.R`, job 17863572) |
| **Fuse / recombination** | DEAD | TS `tree_fuse` is pure clade-exchange (swaps a clade over a bipartition SHARED by two pool trees, `ts_fuse.cpp:217,414`); 85/479 of 1943's clades are in NO near-optimal TS tree → cannot be assembled, only generated; iterative/transient autoconstraint bootstraps the WRONG (1947-basin) backbone (18/64 skeleton, `job 17863085`); rigid consensus autoconstraint compat-killed (locks 4/118 splits absent from 1943) |
| **Random-TBR perturbation strength** | universal landscape property, NOT the discriminator | reach-from-near-miss-basin = 0 at every kick size on BOTH the reachable 4138 floor (`geom4138` job 17862119) and 5432 (`geom_full` job 17860581) |
| **Clip-ordering** (largefirst/antitip/invweight) | clean null | move-ordering sweep, all ~1949–1955 (job 17859608) |
| **Big-pool fuse** (TNT `hold 10000` recipe) | tested to completion, LOST | `tnt_match` (pool=10000, intraFuse, poolSuboptimal=5, sectorMaxHits=20): 8/8 seeds ~4.5h → best 1949, 0/8 ≤1944 (job 17863968) |
| **Deep-narrow** (outerCycles=30, unlimited resets) | converges ~1950 | v2 all 3 seeds stop in 32–50 min at 1950/1950/1954 — *worse* than wide (job 17865387) |
| Suboptimal buffer, addseq variety, sector-size growth, tbrMaxHits 3→20, rasStarts>3 | refuted or null | see `intensity-thorough-hardtail`, `tnt-feature-gap-audit` |

## The ONE lever that moved the floor — the reweighting kick

**`ratchetPerturbMaxMoves` (character-reweighting kick), 5 → ~40–auto, moved 1949 →
1945** (isolation sweep, job 17862445: `perturbMaxMoves` drives it, `perturbProb`
irrelevant; ladder m5→1949, m40→1946, m60(auto)→1946, m160→1945; plateaus ~1945).
This is a **different mechanism** from random-TBR: reweighting ~half the characters
distorts the landscape enough to re-optimize the whole backbone skeleton (the ≥45
splits sectorial is backbone-blind to). It coherently assembles ~49/64 skeleton
splits, then plateaus at 1945 — 2 steps from 1943.

- **Root cause (a probable GENERALIZING recipe defect):** R default
  `ratchetPerturbMaxMoves = 5L` overrides the engine auto value
  `max(20, min(200, n_tip/8))` (= 60 for 482 tips). Every stalled 1949 run kicked
  with 5 = ~12× below auto. The fixed default 5 does not scale with tree size →
  large trees are under-perturbed.
- **Reach benefit was 5432-ONLY (n=1).** Kick replication (job 17862978) on three
  other large matrices (3253/6093/6149) reached their TNT floor at both moves=5 and
  auto — no reach headroom there. So do NOT propose the default change on *reach*
  grounds.
- **Wall-clock benefit on solved matrices — NOW TESTED, and the answer is NO broad
  benefit (job 17872154; `dev/benchmarks/kick_anytime_FINDINGS.md`).** A buildless
  corpus anytime A/B (fixed5 vs auto vs tree-scaled) on the 25-matrix training
  sample (validation sequestered), per-replicate best-vs-wall trace, production
  presets. Result: **14/25 matrices tie exactly**; on the clean large/xlarge tier
  auto/scaled reach the optimum in the **same replicate** (rep2hit ratio 1.000 —
  no gain, no genuine slowdown; the ~7–12% wall2hit gap is costlier-per-replicate
  ratchet, not slower-to-optimum). The over-perturbation-on-small-trees fear is
  **refuted** (small/medium tie; auto=20 on a 25-tip tree costs nothing). A
  **marginal +reach** reappears on one hard multi-basin matrix (project2771, 94t,
  n=5) and a **−reach** on the degenerate rep-starved 4062-tip giant
  (project4284). **=> no broad benefit → do NOT change the default; the kick is an
  opt-in hard-tail intensity lever** (ties `intensity-thorough-hardtail`,
  `auto-vs-thorough-objective`).

## Mission A angles — ALL CLOSED (2026-07-16). Do NOT re-investigate.

All three angles once listed here as open are now refuted **by measurement** (Hamilton
jobs; harnesses in `dev/benchmarks/`; memory nodes cited). Cold-start CONSTRUCTION from
data structure cannot reach the 1943/1939 basin; the residual is a **search-side**
problem, not a start-generator.

1. **Cold-start CONSTRUCTOR (angle #1) — CLOSED, five independent ways.**
   - *Per-split detectability:* data-derived per-split signals (QuartetConcordance,
     PhylogeneticConcordance) do NOT rank true deep splits above wrong — they ANTI-detect
     at depth (238-split PC 9th pct). `detect_backbone.R` (job 17887714). Memory
     `detectability-proportion-signals-refuted-5432`.
   - *Absolute synapomorphy count:* does NOT discriminate true vs wrong deep splits
     (size≥60 AUC 0.51–0.55; 238-split 48th pct). `detect_synap.R`.
   - *Joint / clique recovery (count-matched):* a blind max-compatibility clique recovers
     the true deep splits at CHANCE (36th pctile). `detect_clique.R` (session 9d2149a8).
     Memory `joint-compatibility-clique-refuted-5432`.
   - *Scaffold-sufficiency:* constraining on the TRUE deep backbone (all 25 min-side≥100
     splits, held=1.0) routes in only 1/12; k5/k10/k20 flat near cold; only the FULL tree
     → 1943. `scaffold_suff.R` (job 17893429). Memory `scaffold-sufficiency-refuted-5432`.
   - *Clique-START routability (ground truth):* 12 diverse clique backbones as constraint
     starts → **0/60** hit ≤1944, best 1981 **WORSE** than RAS 1957 (compatibility cliques
     lock homoplastic-but-wrong splits; held≈1.0). `clique_start.R`/`clique_prep.R` (job
     17894813). Memory `clique-structure-5432`.
   - *Clique-proportional reweight (user variant):* ≈ uniform (effective-char ratio 0.83)
     and its mild skew favours clean-SHALLOW chars, down-weighting the ≈4/183 deep-supporting
     chars → reduces to the refuted random resample or worse. `clique_weight_probe.R`.
2. **Basin-HOPPING schedule (angle #2) — CLOSED.** Every near-optimal TS tree (all methods)
   is ≥97 TBR from 1943; ratchet is monotone ILS; no accept-worse schedule walks ~100
   rearrangements (see Route #2 section below). Memory `basin-hop-schedule-refuted`.
3. **Structured breadth (angle #3) — CLOSED.** Random char-resample START-gen lands FARTHER
   than cold (0/400 in capture radius). (job 17887750). Memory `structured-breadth-refuted-5432`.

**CAPSTONE MECHANISM — why every constructor fails.** The 1943 deep backbone has essentially
**no clean character support** (≈4/183 informative chars cleanly support any min-side≥100
split). It is an **emergent aggregate** feature — it wins only on TOTAL tree length,
integrating many individually-homoplastic, mutually-conflicting characters. So no per-character
signal, compatibility clique, or character-reweighting start-generator can detect, construct,
or bias toward it; up-weighting "clean/compatible" structure moves AWAY from it. Only the
aggregate parsimony objective over many restarts finds it (which is what TNT does, ~1/3).

**Where reach work MUST go next (NOT a constructor).** Reach is not a hard limit (TNT finds
1943/1939). The only live lever is **search-side aggregate-objective diversity**: (a) TNT-style
sectorial over a retained diverse pool — already in TS (`ts_driven.cpp`) but does not crack
5432 alone; (b) CID-tabu / diversity-gated restarts (Whelan-2007, UNTRIED; headwind: 5432's
basin is a tiny distant target — `diversity-generation-gates`). Best tree now stored =
**1942** (`floors/project5432_regen_1942.tre`); regenerate + store the 1939 tree when feasible.

## Route #2 taken (2026-07-16): the lever is the ACCEPTANCE FRAME, not the kick knobs

Mechanism read of `src/ts_ratchet.cpp:186–227`: TS's ratchet is a **monotone
strict-descent iterated local search**, not a basin-hop. Each cycle perturbs weights →
short TBR → restores weights → full TBR, then accepts the result **only if strictly
better** (`<`, line 213) and otherwise **resets to the current best tree** (line 220);
kick *strength* (`perturb_max_moves`) is constant across cycles. All three things angle
#2 lists — strength, schedule, multi-kick — operate *inside* that monotone frame, so
none can accumulate the uphill step needed to leave the 1945 local optimum. A hard
plateau at 1945 is exactly the signature of monotone ILS. **The lever is the frame
(accept-worse), not the knobs.**

Plan (advisor-reviewed, cost-gated — do NOT launch the expensive sweep before the gate):

1. **Gate (cheap, buildless — `dev/benchmarks/basin_hop_gate.R`).** Re-score the TNT
   target by label (confirm 1939 vs 1943; RenumberTips gotcha → `TreeLength` by label,
   never edge+tip_data by index) and compute **TBRDist(TS-1945-floor → target)**. Short
   bounded uphill (a few TBR, one ridge) → accept-worse can plausibly bridge it, probe/
   build worth it; 15+ TBR through rugged terrain → accept-worse from 1945 is a random
   walk and route #2 is likely dead (data to hand the user, not a route to override). TS
   *holding* 1945 under TBR already proves it is a true local optimum outside 1943's
   downhill basin, so SOME uphill is required; the distance answers only HOW MUCH.

2. **Buildless probe (only if the gate is favourable).** The whole-tree drift step
   (`src/ts_driven.cpp:452`, `drift_search` on `result.tree`) already composes with the
   ratchet in the replicate loop and already accepts worse moves within
   `driftAfdLimit`/`driftRfdLimit`. `thorough` runs it at `driftCycles=2, afd=5`;
   `default`/`sprint` run none. So the probe is an **afd-limit (uphill-tolerance) sweep**
   past the default, on top of the saturated kick, matched-wall on 5432 (+ project2771,
   the one other kick-sensitive matrix — the only shot past n=1), with every returned
   tree re-scored in R by label. This is the existence test for accept-worse, using
   validated machinery (no build). A phantom is impossible: the reported number is an
   independent R-side `TreeLength`, not the engine's internal accept score.

3. **Gated engine build (only if the probe crosses 1945).** A small, default-OFF
   accept-worse-in-ratchet change to `ts_ratchet.cpp`: in the not-improved branch, accept
   the perturbed local optimum under a schedule instead of always resetting. Keep a
   `current` tree that may wander uphill SEPARATE from a `best_tree` that tracks the true
   best (the classic basin-hop bug: conflating them returns a worse tree than was found).
   A/B on 5432 + project2771. Bar = **opt-in hard-tail lever** (no wall/reach regression
   on the thorough path), not a default — `kick-anytime-not-a-default` already settled
   that a landscape-specific kick lever is not a broad default.

### GATE RESULT — route #2 REFUTED (2026-07-16, buildless; no probe/build run)

`dev/benchmarks/basin_hop_gate.R` re-scored the near-optimal TS trees by label (ordered-char
scorer; self-check `L(1943 floor)==1943` PASSED) and measured `TBRDist(exact=FALSE)` (hard
lower bounds) to the 1943 floor. **Every generation method's near-optimal tree is `tbr_min`
≥97 TBR from 1943:**

| method | score (+gap) | tbr_min..max | CID |
|---|---|---|---|
| kick p025m160 | 1945 (+2) | **97**..279 | 0.220 |
| invweight | 1950 | 97..282 | 0.272 |
| strongperturb | 1945 (+2) | 100..286 | 0.229 |
| p025m40 | 1946 (+3) | 101..293 | 0.269 |
| cold ts5432 | 1958 | 102..299 | 0.379 |
| largefirst | 1950 | 107..310 | 0.298 |
| random | 1949 | 116..335 | 0.361 |

The **best kick tree is only 2 score-steps above 1943 yet provably ≥97 TBR away**, barely
closer than a cold tree (97 vs 102). Corroborated: basin_5432 descent returns to 1943 only
for d≤12, **0/5 by d≈20** (the floor is ≥97 out); geom_full ~109. **No accept-worse schedule
— drift OR iterated reweighting — traverses ~100 rearrangements to a narrow basin unreachable
from 20 TBR out.** `afd_limit` is a per-move score tolerance, not a distance budget. The
reweight-schedule loophole is closed empirically (all method-diverse near-optimal trees, incl.
reweight-space, are ~97–116 TBR out). The afd-sweep probe and the accept-worse build were NOT
run — the ≥97 bound guarantees their outcome.

**Kick is not useless:** CID 0.22 (kick) < 0.38 (cold) — the kick moves meaningfully closer in
information terms; it improves score by finding a **better LOCAL basin**, just bottoms out ~97
TBR short. **Scope of the refutation:** route #2 only (escape the floor via *any* accept-worse
schedule). Angle #1 (a cold-start CONSTRUCTOR that GENERATES a start inside 1943's basin) and
angle #3 (structured breadth) are UNTOUCHED — they do not walk from a floor, so the distance
wall does not apply. Memory `basin-hop-schedule-refuted`.

## Status

- **Diagnosis: cold-start basin capture** (operator-capable, trajectory-limited).
  Confirmed by observation (head-to-head), not inference.
- **Build question (freeze-a-seed operators): CLOSED.** BSS structurally refuted;
  fuse structurally cannot assemble ungenerated clades; every ported TNT feature
  built+A/B'd → refuted or matched.
- **5432 arm: cold-start CONSTRUCTION CLOSED (2026-07-16)** — all constructor / schedule /
  breadth angles refuted (see "Mission A angles — ALL CLOSED" above; capstone: the deep
  backbone is an emergent aggregate with ≈no clean character support, so no start-generator
  can bias toward it). Reach is NOT a hard limit (TNT ~1/3); the only remaining lever is
  **search-side** (sectorial-over-diverse-pool, already in TS; or UNTRIED CID-tabu restarts) —
  a different program from a start-generator. Do NOT re-open the constructor line.
- **Mission read:** not a reach *regression* (thorough's disqualifier) and not a
  wall-clock regression on anything else — one 482-tip matrix where TNT itself is
  1/3-reliable. Whether the last 1–2 steps justify a new engine is a user-level
  mission call, not an autonomous one.
- Diagnostics, all reproduced: `restore238.R` (TBR rebuilds 1943 from capture
  radius), `basin_5432.R` (capture radius ~12), `geom_full.R` / `geom4138.R`
  (landscape geometry, kicks universal-null), `ceiling.R` (BSS refuted),
  `tnt_match.R` (fuse recipe → 1949), `ts_multi` (depth → 1947), `seed1943.R`
  (HOLD = pure reachability), `clade_check.R` (85/479 clades absent from cold
  trajectories). Durable artifacts on Hamilton `/nobackup/pjjg18/`:
  `floors/project5432_tnt_floor.tre` (TS-1943), `tntwatch/`, `reeval/`. Matrix
  `.../lgsweep/matrices/project5432.nex` = TNT `m01.tnt`; floor recipe
  `sect:slack 20; xmult=replic 120 hits 20; bbreak=tbr`.
- Full diagnostic chain: memory `project5432-basin-structure` (read the whole node).
- **No code change proposed here; do not commit this note as part of a code branch
  (concurrent-session git hazard — stage the named file only if the user asks).**
