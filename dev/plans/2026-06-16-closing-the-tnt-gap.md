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
- *0c (IN PROGRESS)*: gap-panel baseline → `dev/benchmarks/headtohead_phase0.csv`.

**Phase 1 — phase-yield diagnosis.** Add a per-phase candidate counter (mirror
`PhaseTimings`/`ph_lap`) so we can attribute candidates-per-improvement to each phase
(sectorial vs ratchet vs final-TBR vs fuse). Localises *where* the 6× lives before building.

**Phase 2 — attack sectorial search**, in priority order:
1. Wire + **gate `accept_equal`** (built, off). Goloboff 2014 §5.4 flat-landscape lever
   (1881/2500 vs 1/2500 MPT on missing-data matrices). Gate to engage on plateaus.
2. **Drift done right**: large `numsub` (Goloboff Fig 8C) + equal-score acceptance —
   the untested combination (drift is off in all presets and was only tested without these).
3. Exact HTU sector scoring (kills the rescore/miss path) + adaptive sector sizing.

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
