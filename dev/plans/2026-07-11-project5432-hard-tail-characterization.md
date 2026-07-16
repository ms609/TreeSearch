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

## Genuinely-open angles (Mission A's search space — none refuted)

1. **A cold-start CONSTRUCTOR that generates starts near the target basin.** BSS
   (one instantiation) is structurally dead, but the goal stands: constraint-guided
   / char-informed addition sequence; progressive assembly; any start-generator that
   lands TS within ~8–12 TBR of the best basin.
2. **Push the reweighting kick as a basin-HOPPING schedule** (the one thing that
   worked) — strength, schedule, multi-kick — beyond the single-knob ratchet
   default. Plateaus at 1945; can a schedule cross the last ~2 steps?
3. **Smarter breadth, not naive.** 2150+ uniform cold reps = 0 hits → uniform
   restarts refuted; open question is whether STRUCTURED / diversified starts raise
   the per-rep basin-capture rate.

## Status

- **Diagnosis: cold-start basin capture** (operator-capable, trajectory-limited).
  Confirmed by observation (head-to-head), not inference.
- **Build question (freeze-a-seed operators): CLOSED.** BSS structurally refuted;
  fuse structurally cannot assemble ungenerated clades; every ported TNT feature
  built+A/B'd → refuted or matched.
- **5432 arm: OPEN** only for a genuinely-new mechanism that lands cold starts (or a
  trajectory) within the ~8–12 TBR capture radius on large trees — a constructor or
  a basin-hopping schedule, not a config of the exhausted levers.
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
