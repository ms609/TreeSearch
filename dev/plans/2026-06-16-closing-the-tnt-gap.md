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

## Phase 3 design (2026-06-16) — structural options scoped; cheap falsifiable probe first

A 4-agent design workflow assessed three structural options to cut candidates-per-improvement;
the user opted to commit to Phase 3, so it was scoped before any code.

- **Branch-collapsing / full polytomy search (Goloboff 2023): REJECTED for this gap.** A 3–6 week
  tree-representation + Fitch-kernel + TBR-clip rewrite touching the most-optimized code in the repo
  (binary `left/right` in `ts_tree.h`, the 2-input SIMD primitives in `ts_fitch`, TBR clip/regraft).
  It attacks the wrong axis (frugality, not the escape/depth ratchet owns), and its mechanism barely
  fires: the project already measured ~0% collapsed-edge rate on near-optimal binary morphological
  trees, and the advisory collapsed-flag skip (`ts_tbr.cpp:817-820,919-921`) + `add_collapsed` pool
  dedup already bank the easy ~80%.
- **Exact-scoring sectorial (CSS): ALREADY ACTIVE on the gap datasets.** `css_search`
  (`ts_sector.cpp:1005-1073`) runs full-tree TBR restricted to a `sector_mask` — exact by
  construction, no HTU pseudo-tip, no miss-and-revert. `thorough` sets `cssRounds=2`, `large`
  `cssRounds=1`, and `auto` routes 65–119t→thorough / ≥120t→large. So "implement exact sectorial"
  is largely already done; the residual is a ratchet→CSS budget-**rebalance experiment** (days), not
  a kernel rewrite.
- **Union-based region-merging (Lever 1): cheap (days) but likely a no-op.** `compute_collapsed_regions`
  (`ts_collapsed.cpp:106-170`) is built but DEAD CODE (zero callers); wiring it merges equal-resolution
  regraft positions. The 0%-collapsed-rate finding predicts it barely fires at the optimum on the hard panel.

**Gap-closure risk: HIGH.** All three reduce candidates-per-improvement (frugality) but none finds
lower-score basins (the escape/depth axis ratchet owns at 63–83% of wall-clock). The +1/+3 most likely
remains after any rewrite — consistent with the Phase 2 floor finding.

**Decision (data-gated):** do the smallest structural slice that decides the rest — a falsifiable probe:
(1) ratchet→CSS rebalance sweep on the gap panel (`p3_rebalance.csv`); (2) a collapsed-region/edge-rate
probe at the optimum (Lever-1 go/no-go). Escalate to wiring Lever 1 ONLY if the rebalance beats baseline
(no per-dataset regression) AND regions are non-trivial; otherwise confirm the floor and rest on the
shipped (a) accept-floor + (b) opt-in `intensive` preset. Branch-collapsing is pursued only if both
probes reveal a large, real collapsed signal the advisory path leaves on the table — which the existing
data predicts they will not.

### Phase 3 outcome (2026-06-16) — structural search rewrite NOT justified; pivot to per-candidate wall-clock

The ratchet→CSS rebalance probe (`p3_rebalance.csv`, 3 seeds) is **FLAT**: every config that trades
ratchet budget for exact CSS saves 28–51% candidates but **regresses the hard datasets** (Zanol +4,
Giles +2); none beats baseline without a per-dataset score regression. With the design verdict
(branch-collapsing wrong-axis + ~0% collapsed rate; CSS already active on the gap datasets) this is the
third convergent confirmation that **the EW-Fitch score gap is at the practical floor and is
landscape/escape-bound, not frugality-bound** — no structural *frugality* lever closes it. The Lever-1
region-merging precondition (rebalance must beat baseline) failed, so it is not pursued; branch-collapsing
is rejected.

**Pivot — the movable lever for the original ~2× WALL-CLOCK concern is per-candidate COST, not count.**
The frugality analysis surfaced it as "option 4": VTune (`vtune_tbr_analysis.md`, T-260) puts StateSnapshot
save/restore at ~23% of TBR time (a full ~190 KB memcpy per accepted/rejected move) and a redundant
`reset_states` `std::fill` at ~4%. `apply_tbr_move` already knows the dirty nodes, so selective save/restore
of only those rows is est. **10–16% wall-clock**, **dataset-agnostic, no score trade-off** — it cuts the
time per candidate rather than the candidate count, which is orthogonal to the score floor and directly
targets wall-clock.

**Correction on inspection (do NOT act on the stale figures above):**
- The `reset_states` `std::fill` (design "fix #2", ~4%) was **already removed in T-261**
  (`ts_tree.cpp:265-277` — "every array entry read is written before it is read"). That win is banked.
- The cited StateSnapshot ~23% comes from a VTune doc that **predates T-261** (it *recommended* the
  fill removal T-261 then made) and likely T-300's incremental-SPR accept path — so the figure is
  **stale and the share has probably shrunk**. The remaining lever (selective `StateSnapshot`
  save/restore) is intricate, correctness-critical surgery on the most-optimised code in the repo.
- **Decision:** it must be **re-profiled in a fresh `/profile` (VTune) round** to confirm it is still a
  meaningful hotspot *before* the surgery — not done on stale data at the tail of this round. Verification
  when pursued: behaviour-neutral via **candidate-identity** (a correct timing optimisation must leave
  `candidates_evaluated` and scores bit-identical vs baseline on the iterate gate) + a wall-clock
  micro-benchmark for the delta; keep only if identical-and-faster, else revert.

**Phase 3 — branch-collapsing search** (Goloboff 2023): search the reduced polytomy
tree space, not just skip candidates/dedup as now. Structural swing; pursue only if
Phase 1/2 data shows the candidate-frugality gap justifies it.

## Challenge 2 closeout (2026-06-17) — ratchet now genuinely disableable; ratchet-OFF still trails TNT

The "ratchet is untouched / disable it to match TNT" thread (user Challenge 2) is resolved.

**Ratchet was never disableable.** `ratchetCycles = 0` still ran ratchet via three
stacked floors in `ts_driven.cpp` (ceiling-division `max(1, …)`, an unconditional call
site, and the `adaptive_level` re-floor `max(1, base * scale)`); `ratchet_search` also
runs an initial TBR pass before its cycle loop. All three are now guarded — a no-op for
every preset (all use `ratchetCycles ≥ 3`), covered by `test-ts-ratchet-disable.R`.

**With ratchet genuinely off, TreeSearch does NOT match TNT.** Patched build,
`adaptiveLevel = FALSE`, TNT-matched core, 4 datasets × 5 seeds, only `ratchetCycles`
varied (`bench_ratchet_axis.R` → `ratchet_axis.csv`). Median gap to TNT `xmult`
(arm − TNT, lower = better):

| dataset | TNT | R0 (true off) | R1 | R12 | gap R0 | gap R12 |
|---|---|---|---|---|---|---|
| Giles2015   | 670  | 675  | 675  | 672  | +5 | +2 |
| Wortley2006 | 480  | 485  | 487  | 482  | +4 | +2 |
| Zanol2014   | 1262 | 1269 | 1268 | 1267 | +8 | +5 |
| Zhu2013     | 624  | 631  | 631  | 629  | +7 | +5 |

- Ratchet-off (R0) trails TNT by **+4…+8** on every dataset; even our single best
  ratchet-off seed never reaches TNT's median. The deficit is **not** ratchet-caused
  (TNT runs no ratchet on these either) → it is structural.
- Ratchet helps **monotonically**: R12 closes the gap to +2…+5 (−2…−3 vs R0) at ~2–3×
  wall-clock. "Disable ratchet to match, then switch on to pull ahead" inverts reality —
  ratchet is *necessary to approach* TNT; it narrows but never erases the deficit.
- The residual gap is **sectorial / fusing search efficiency** — re-examined against the
  published algorithm in the Goloboff-1999 divergence analysis (2026-06-17,
  `dev/plans/2026-06-17-sectorial-divergence.md`).

Memory: `ratchet-not-disableable.md`. (Local TNT is 32-bit, so its wall-clock is not a
fair reference; scores / rearrangement counts are. R0-vs-R12 wall-clock is comparable.)

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

## Phase 4 (2026-06-18) — UPDATE: the "floor" was a cost-formula bug; score gap CLOSED

The Phase 1–3 / Challenge-2 conclusion that the EW-Fitch gap was a
**"landscape/escape-bound floor, not frugality-bound"** with a "competitive
per-candidate kernel" is **superseded**. Those phases predated finding a
correctness bug in the candidate insertion-cost function.

**Root cause (see `2026-06-18-wagner-insertion-cost-bug.md` + memory
[[wagner-insertion-cost-bug]]).** `fitch_indirect_length*` scored a candidate
insertion edge with the **union of the two endpoints' final states**
(`final[A] | final[D]`), which is not the exact edge set — it undercounts on
ambiguous trees and mis-ranks/mis-cuts moves on resolved ones. This (a) made RAS
Wagner starts **~+30% over the optimum** (near-random greedy placement) and
(b) gave TBR wrong cost magnitudes → wrong cutoffs / early abandonment. The fix
is the exact **directional** edge set `E[D]=combine(prelim[D],up[D])`
(`compute_insertion_edge_sets`, ts_fitch). Shipped to `wagner_tree`, the EW
`tbr_search` SPR scan + rerooting vroot, and `build_ras_sector` (commits
2b299e4b, 93071cae on cpp-search).

**Result — full `thorough` search now reaches the MPT across the hard panel**
(`diag_gap_panel_postfix.R`, 60s, 3 seeds):

| dataset | target | post-fix (min / median) | pre-fix floor (Phase 2 / Challenge-2 R12) |
|---|---|---|---|
| Wortley2006 | 480  | **479 / 479 (−1)** | +3 / +2 |
| Zanol2014   | 1261 | **1261 / 1261 (+0)** | +3.5 / +5 |
| Zhu2013     | 624  | **624 / 624 (+0)** | +2.5 / +5 |
| Giles2015   | 670  | **670 / 670 (+0)** | +1.5 / +2 |

So the gap the plan repeatedly called "structural / escape-bound" was, in
substantial part, this scoring-formula bug. Gap **B is now ~0** (was +1.5..+3.5).

**Candidates-per-improvement (gap C) reversed.** Vinther2008 (the canonical tie,
pre-fix "6.3× *more* candidates than TNT"): post-fix `bench_tnt_headtohead.R`
gives TS 78 = TNT 78 with **cand_ratio 0.44** — TreeSearch now examines *less
than half* TNT's rearrangements to reach the same score (counters tally slightly
different events, so indicative — but a qualitative reversal). Wall-clock tied on
this small case (0.4s = 0.4s).

**Disposition of the open threads:**
- **Core TBR/Wagner hill-climbing deficit (task #26): RESOLVED** — root-caused and
  fixed; panel now at +0/−1.
- **Drift-done-right for "+1 datasets" (task #25): MOOT** — the +1 datasets
  (Zanol/Giles) it targeted are now +0; no score gap remains for drift to close.
- **Race-to-common-target (task #22): target reached** across the panel;
  candidates-per-improvement competitive-to-better. The *only* residual is the
  authoritative **wall-clock ratio**, which is **Hamilton-gated** (local TNT is
  32-bit) — the still-live wall-clock thread, not a quality gap.

**Remaining (genuinely open):** authoritative wall-clock vs 64-bit TNT on
Hamilton (the original ~2× concern — now partly addressed by the frugality
reversal, but unconfirmed on large datasets), and the chip's TBR
move-completeness fix (L812/nz/ns; small, poor-start-only — see
[[tbr-rooted-vs-unrooted]]). The per-candidate `StateSnapshot` micro-opt
(above) is independent and still available.
