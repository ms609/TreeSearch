# Closing the TNT Gap — Strategic Plan

Branch: `cpp-search` · Reference: TNT 1.6 · Comparison path: equal-weights Fitch,
apples-to-apples (`-` → `?`). Drafted 2026-06-16. Supersedes the retired
`.positai/plans` strategy thread (see `dev/plans/README.md`).

## Goal

Close the **wall-clock** gap to TNT 1.6 on equal-weights Fitch parsimony. TNT
reaches a comparable-quality tree roughly **2× faster**. The kernel is already
refined (we are competitive-to-faster *per candidate*), so the lever is search
strategy — specifically, how many candidate rearrangements we burn per unit of
score improvement.

## Reframe: three gaps, not one

| Gap | What it measures | Magnitude | Target? |
|-----|------------------|-----------|---------|
| **A. Scoring method** | Brazeau three-pass vs TNT column-Fitch on *inapplicable* data | +1 … +50 (e.g. Vinther raw 79 vs 78) | **No** — a different, arguably better objective; vanishes under `-`→`?`. |
| **B. Score quality** @ fixed time, apples-to-apples Fitch | TNT finds a shorter tree, same budget | +2/+3 steps, hardest datasets | Small; perturbation-tuning lever now **spent**. |
| **C. Wall-clock** to a comparable score | TNT is ~2× faster | ~1.5–3× | **Yes — the prize.** |

Empirical confirmation of A vs B in one row (Phase 0, Vinther2008): TreeSearch
Fitch **78** = TNT **78** (gap B = 0); TreeSearch *raw* (Brazeau) **79** (gap A = +1).

## Diagnosis: C is a candidates-per-improvement gap, concentrated in sectorial search

`wall-clock = (cost per candidate) × (candidates per unit of improvement)`.

- **Cost per candidate**: competitive. The raw Fitch kernel may lead TNT; in-search
  per-candidate cost carries StateSnapshot/rescore overhead (T-260), but this is
  not where the 2× lives.
- **Candidates per improvement**: we are far worse. **First instrumented measurement**
  (`bench_tnt_headtohead.R`, candidates-per-improvement mode): on Vinther2008, at the
  *same* score (78), TreeSearch evaluated **2.90M** candidate rearrangements vs TNT
  **0.46M** — **6.3×**. Even on a tie we burn 6× the work.
- **Where**: TNT's `xmult` is ~67% sectorial search. Our sectorial search runs, but
  with equal-score acceptance **hard-coded off** (`ts_sector.h:24 accept_equal=false`,
  never set true on the `MaximizeParsimony` path) and approximate HTU scoring.

## Plan (phased, data-gated)

**Phase 0 — instrument + baseline.**
- *0a (DONE)*: `candidates_evaluated` — total TBR/SPR-class candidates, accumulated in
  `tbr_search` into `DataSet::n_candidates_evaluated`, summed over a serial
  `driven_search`, surfaced as `attr(MaximizeParsimony(...), "candidates_evaluated")`.
  Behaviour-neutral; valid `nThreads=1` only; excludes NNI-warmup/annealing.
- *0b (DONE)*: `dev/benchmarks/bench_tnt_headtohead.R` — TreeSearch (Fitch + raw) vs
  TNT, capturing score, candidates, TNT rearrangements, wall-clock; separates gaps A/B/C.
- *0c (DONE)*: gap-panel baseline (`headtohead_phase0.csv`, 2 seeds, converged).
  **Gap B = 0..+3.5 steps** (Zanol +3.5, Wortley +3, Zhu +2.5, Giles +1.5, Dikow/Eklund 0).
  **Gap A** (raw − Fitch) = +50/+39/+12 on high-inapplicable Zanol/Giles/Zhu, 0 on
  Wortley/Eklund — pure scoring method, tracks inapplicable fraction. **Candidates-per-
  improvement ~1.3–1.9×** on real datasets (the Vinther 6.3× was a tiny-dataset outlier),
  and TNT lands a *better* score — more efficient on both axes. Wall-clock ≈1.3–2.5× vs
  *32-bit* TNT (larger vs fair 64-bit), ≈ or above the candidate ratio (per-candidate
  overhead is a co-contributor).

**Phase 1 — phase-yield diagnosis (DONE; REDIRECTS Phase 2).** `bench_phase_yield.R`
(per-phase wall-clock share + total candidates + `late_frac`):
- **Ratchet dominates wall-clock: 63–83%** (Wortley 83, Eklund 76, Dikow 68, Giles 66,
  Zanol 66, Zhu 63). **Sectorial is minor: 7–23%.** final-TBR 2%, init-TBR 2–5%, fuse 0–2%.
- **44–98% of replicates land AFTER the last improvement** (`late_frac`: Eklund 0.98,
  Giles 0.90, Wortley 0.81, Dikow 0.79, Zanol 0.61, Zhu 0.44) — large post-convergence waste.
- **Implication:** the cost centre is RATCHET, not sectorial — the *opposite* of TNT
  (~67% sectorial). We pour 66–83% of wall-clock into ratchet and still finish +2/+3 worse.
  *Caveat (to verify):* this is wall-clock share; sectorial's reduced-dataset candidates are
  cheaper per candidate, so a **per-phase CANDIDATE counter** is the next instrumentation to
  confirm ratchet also dominates candidate *count*, not just clock. `adaptive_level` likely
  scales ratchet_cycles UP on stalled (hard) datasets — pumping effort into the wasteful phase.

**Phase 2 — REDIRECTED by Phase 1 data. Experiments, cheapest first, each gated on
candidates-per-improvement + score vs baseline; default-off until validated:**
1. **Cut wasted ratchet.** 44–98% of reps are post-convergence. Test tighter stopping
   (`perturbStopFactor`, `targetHits`), fewer `ratchetCycles`, and capping/disabling the
   `adaptive_level` ratchet up-scaling on stalled datasets. Cheapest, biggest wall-clock lever.
2. **Rebalance ratchet → sectorial.** Shift budget toward sectorial (TNT's efficient phase):
   more `xss/rss` rounds, fewer ratchet cycles.
3. **Make sectorial plateau-capable** so leaning on it pays: wire + gate `accept_equal`
   (`ts_sector.h:24`, built/off) — Goloboff 2014 flat-landscape lever for high-inapp
   Zanol/Giles; + drift-done-right (large `numsub` + equal acceptance).
- *Next instrumentation:* per-phase candidate counter (mirror `PhaseTimings`/`ph_lap`) to
  confirm the clock→candidate correspondence before committing to a ratchet rewrite.

### Phase 2 results (2026-06-16) — cheap/medium levers tested, no robust global win

Via `bench_p2_levers.R` (gap panel, fixed 20 reps; the fast loop made each round ~90s).
Deltas vs the `auto` baseline (`iterate_baseline_auto.csv`):

- **Ratchet/sectorial knobs (round 1, `p2_levers.csv`):** `ratchetCycles` {3,6},
  `adaptiveLevel=off`, `xss/rss` rounds doubled, ratchet→sectorial `rebalance` — **none beats
  baseline.** Cutting ratchet saves 30–60% candidates but costs +0.5–2.5 steps on hard
  datasets (ratchet does real work); more sectorial rounds tie-or-worsen; `adaptiveLevel=off`
  is exactly neutral. The `auto` preset is near-Pareto-optimal for these knobs.
- **`accept_equal` (the #1 untried lever; hard-coded on via the fast loop, then reverted):**
  neutral-to-worse (Zanol +3, Giles +1), candidates barely move (0 to −3%). **Why it fails
  here:** sectorial is only 7–23% of our wall-clock (Phase 1), so its acceptance criterion has
  little leverage — the opposite of TNT (~67% sectorial). The built-but-off infrastructure is
  not the lever *for our pipeline shape*.
- **Fusing/ordering/starts (round 2, `p2_levers_fuse.csv` 2-seed, `p2_fuse_5seed.csv` 5-seed):**
  5-seed medians confirm a *real but per-dataset* signal: **`wagnerStarts=5` and `intraFuse`
  each robustly improve Wortley (−3/−2) and Zhu (−2/−3, → 626 vs TNT 624)** but **regress
  Zanol/Giles by +1**; Eklund/Dikow neutral. `fuseAcceptEqual` ≡ `intraFuse`. `clipOrder=2`
  saves 22–32% candidates at +1–2 steps (worse). **No feature cleanly separates helped (Wortley
  37t/8st, Zhu 75t/4st) from hurt (Zanol 74t/9st, Giles 78t/4st)** — so no safe global default.

**Conclusion — apples-to-apples Fitch gap is at the practical parameter-tuning floor.** No
single config improves all panel datasets; the only real gains (−2/−3 on Wortley/Zhu) are
dataset-specific and come with +1 regressions elsewhere, failing the "no regression on any
dataset" ship gate. `accept_equal` (the headline untried lever) has no leverage in our
ratchet-dominated pipeline. Remaining options, by cost: **(a) accept the floor** — declare the
EW-Fitch gap effectively closed (+1/+3 on the hardest datasets), redirect effort; **(b) ship an
opt-in variant** (`intraFuse`/extra Wagner starts in `thorough`) so the Wortley/Zhu wins are
available without touching `auto`; **(c) Phase 3 structural** (branch-collapsing / exact-scoring
sectorial) — weeks-scale, the only thing that could move a ratchet-dominated pipeline toward
TNT's per-candidate frugality, but hard to justify for a residual +1/+3 steps. Recommendation:
**(a)+(b)**, not (c) — the data does not justify a weeks-scale structural rewrite for this gap.

**Phase 3 — branch-collapsing search** (Goloboff 2023): search the reduced polytomy
tree space, not just skip candidates/dedup as now. Structural swing; pursue only if
Phase 1/2 data shows the candidate-frugality gap justifies it.

## Methodology guardrails

- **Optimise against candidates-per-improvement** (continuous, low-variance), not
  score-at-fixed-time (±2-step lottery on a small panel).
- **Authoritative wall-clock**: Hamilton 64-bit Linux TNT (matches the on-disk
  `t264`/`t249` reference scores). The local `tnt.exe` is **32-bit** (PE32/i386) —
  its *scores and rearrangement counts are valid* (bitness-independent), but its
  *wall-clock is not* a fair reference for our 64-bit build. (User believes a Win64
  TNT exists locally; not found at the standard path — to confirm.)
- Validate on the **MorphoBank validation split**, not just the 14 CRAN datasets.
  Report median + min over ≥5 seeds. `nThreads=1`. Everything default-off until gated.

## Artifacts

- Harness: [bench_tnt_headtohead.R](../benchmarks/bench_tnt_headtohead.R)
- Metric: `attr(MaximizeParsimony(...), "candidates_evaluated")` (serial)
- Baseline data: `dev/benchmarks/headtohead_phase0.csv`

## Decisions / dead ends (do not re-propose)

- **Perturbation escalation** (`stallEscalateFactor`, shipped 2026-06-16): score-neutral
  on Wortley/Zanol → this vein is **spent**. It ships off-by-default as a stall safety net.
- **Static perturbation re-tune** (`ratchetPerturbProb=0.15`): refuted — regresses
  Zanol/Zhu by +9/+11.
- **Drift / NNI-perturb / prune-reinsert / adaptiveStart in presets**: recorded-negative
  (T-274, T-289f, T-190); out of scope except drift's *specific* untested combination above.
- **Raw speed** (AVX2/popcnt): real but won't close the strategic gap.
