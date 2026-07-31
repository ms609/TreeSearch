# `targetHits` reach escalation — ship gate and reach study (2026-07-24 … 07-31)

Durable record for the escalation shipped in `MaximizeParsimony()`: doubling `targetHits`
or more, under `strategy = "thorough"`/`"large"`, additionally deepens the per-replicate
perturbation.  Harness: `reach_escalation_ab.R` + `reach_escalation_analyze.R` +
`reach_escalation_ab.sh` (this directory).  Raw per-cell CSVs are gitignored; the numbers
below and the committed harness are the reproducible record.

## The levers

Six `SearchControl` fields, applied as a flat bundle when
`targetHits / max(10, nTip/5) >= 2`, skipping anything the caller set:

    ratchetPerturbMaxMoves = 0   # auto/deep reweighting kick
    driftCycles            = 25
    postRatchetSectorial   = TRUE
    stallEscalateFactor    = 1.5
    intraFuse              = TRUE
    poolSuboptimal         = 3   # internal only; filtered out of the returned trees

`ratchetCycles` is deliberately NOT here — ratchet depth belongs to `.IwRatchetDepth()`,
which scales it continuously against its own 36-matrix calibration.  A flat value here
would clobber that under implied weights, and a separate 68-matrix comparison found ratchet
depth gives no equal-weights reach gain (0.970 vs 0.965), so equal weights loses nothing.

## Ship gate — general-pool A/B

**Design.**  `MBANK_FIXED_SAMPLE` (25 training matrices, 20–4062 tips) × 5 seeds = 125
cells, one SLURM task per cell with both arms on the same node.  Validation split
sequestered (asserted per matrix).  Regime EW Fitch, gaps→missing.  `strategy = "auto"`,
so the escalation layers on whatever preset auto picks — the faithful test.

**The confound this avoids.**  Comparing default-`targetHits` against raised-`targetHits`
would confound the levers with "raised `targetHits` searches longer anyway".  So *both*
arms run at the SAME raised `targetHits` (exactly the trigger threshold) and differ ONLY in
the six levers.  Run on a gate-free engine with the levers passed as dots, so the result
does not depend on the gate implementation being correct.

**Target.**  Union-best final score across arms within each cell (the established mbank
convention — there is no canonical best-known table for these matrices).  Note this is
self-referential: if one arm alone attains a score, the other "misses" by construction, so
the reach fractions restate the paired counts rather than measuring absolute optimality.
The paired win counts below are the honest statistic.

**Result.**  No cell in either arm hit its wall cap, so nothing is a budget artefact.

Strict paired final-score wins (deltas vs base):

| tier | n | deltas better | base better | tie |
|------|---|---------------|-------------|-----|
| small (≤30 tips) | 35 | 0 | 0 | 35 |
| medium (31–60)   | 35 | 0 | 0 | 35 |
| large (61–120)   | 35 | 1 | 1 | 33 |
| xlarge (≥121)    | 20 | 5 | 0 | 15 |

**Read the xlarge row carefully: all five wins are the same matrix**, `project4284` at 4062
tips, which improved on every one of its five seeds (by 1–9 steps).  The other three xlarge
matrices (125, 131, 173 tips) all tied.  So the demonstrated benefit is *datasets far too
large to converge within an ordinary budget*, **not** a property of "over 120 tips".

Cost, median over cells: total wall ×3.56, time-to-best ×2.49, but **replicates-to-best
×1.00**.  The wall gap is entirely work per replicate (~15× candidates evaluated against a
`sprint` baseline; ~1.5× against `thorough`), not slower convergence.  Reporting wall alone
would read as a 2.5× regression; `rep2hit` is what shows it is not.  (Same lesson as
`kick_anytime_FINDINGS.md`.)

The one loss: `project2771` (large, seed 5821), base 911 → deltas 912.  Against 6 better /
1 worse overall and a tied `large` tier, read as stochastic — but recorded, not swept away.

**Pre-registered rule** (fixed before results): ship iff reach does not regress and no tier
regresses; wall cost is accepted by construction, since the escalation only fires when the
user has asked for it.  → **SHIP.**

## Why it is gated, and gated to `thorough`/`large`

Below ~200 tips the levers buy no reach and cost ~3.5× the wall, so they must not be
defaults.  Gating on `targetHits` puts the cost only where the user asked for it.  Scoping
to `thorough`/`large` matches `.IwRatchetDepth`, keeps `sprint`/`default` at their own
measured implied-weights operating point (`.iwStopPackage`), and costs nothing measured —
the small and medium tiers were 0 better / 0 worse across 70 cells.

## Hard-tail reach study (project5432, 482 tips, EW; the levers' origin)

Best known 1943 (shared TS/TNT floor; a 1942 tree exists from a rare TNT run).  All arms
`targetHits = 999`, `maxReplicates = 500`, matched seeds, current engine.

| arm | configuration | reach 1943 |
|-----|---------------|------------|
| A | stock `thorough`, single call | 0/3 (floor 1944) |
| B | + the levers, single call | 0/3 (floor 1944) |
| C | + the levers, R block loop (`tree = best` carry-forward, ~48 re-entries) | **5/16** |
| D | + the levers, `TS_POOL_RESEED=0.5`, single call — **time-truncated** | 2/16 |
| E | as D, un-starved (`enumTimeFraction = 0`, 68 h) | **4/16** |

**Mechanism.**  A supplied `tree =` seeds only replicate 0 (`src/ts_driven.cpp`), and
nothing inside a single call ever re-optimises the global best through the deep pipeline
again — the restart strategies are all from-scratch.  So a single call re-solves the
incumbent once; the block loop does so once per block.  Reaching the floor is a
restart-volume phenomenon, and `outerCycles`/`maxOuterResets` are not a substitute (they
were identical in every arm).

**Vehicle verdict.**  Arm C 5/16 vs arm E 4/16 is statistically indistinguishable (Fisher
p ≈ 1.0), so the in-engine `TS_POOL_RESEED` matches an external block loop while keeping
the pool, conflict table and fuse donors that each re-entry discards, and paying the
MPT-enumeration reserve once instead of per block.  **Reach is ~25–31% on any vehicle: the
recipe improves the odds, it does not guarantee the floor.**

## Methodological traps worth not re-learning

- **`maxSeconds` is not the search budget.**  `main_deadline = maxSeconds × (1 −
  enumTimeFraction)`, default 10% reserve.  Identical elapsed times across seeds is the
  time-truncation signature; arm B and arm D were truncated exactly there while appearing
  replicate-capped.  For pure reach runs set `enumTimeFraction = 0`.
- **Whole-suite test runs need `test_local()`/`devtools::test()`**, not `test_dir()` +
  `library()`: several test files call internals unqualified and error otherwise.
- **A top-level `skip_on_cran()` makes a file report "0 pass 0 fail"** — that is *not* a
  pass, it never ran.  Use `NOT_CRAN=true`.
- **To attribute a suite failure, use a pristine detached worktree at the parent commit.**
  Swapping a single source file is invalid once you have added a test file, because your own
  tests then error on the missing functions and the totals stop being comparable.

## Not measured

- The shipped configuration was never A/B'd *as such*: the run above tested seven levers
  ungated; six ship, gated to `thorough`/`large`.  The inference is a tier decomposition of
  that run (the xlarge win is under `large` and preserved; the excluded tiers measured 0
  better / 0 worse, so the gate should remove only cost), not a direct measurement.
- Equal weights only.  Use under implied weights and profile parsimony is a deliberate but
  unmeasured extrapolation; ratchet depth there is governed separately.
