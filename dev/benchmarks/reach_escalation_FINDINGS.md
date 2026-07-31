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

**Result.**  Budget regime (re-derived 2026-07-31 with the *fixed* deadline detector — the
original write-up said "no cell in either arm hit its wall cap", which was **wrong**, an
artefact of the `0.95 × cap_s` bug): 36 base and 40 deltas cells stopped at the deadline, 85
cells converged in both arms, and **4 cells are asymmetric** (project2184 seeds 5821/5823/
5824/5825 — deltas at the deadline, base not).  All four tied at 563, so no score comparison
here rests on an asymmetric budget, but the clean-sweep claim does not stand.

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

Per-matrix concentration for this run (added 2026-07-31): only **2 of the 25 matrices changed
at all** — `project4284` 5 win / 0 loss / 0 tie, and `project2771` **1 win / 1 loss / 3 tie**.
So the `large` row's "1 better, 1 worse" is *both* project2771, i.e. the one matrix that is
demonstrably high-variance in both runs.  **project4284 is the only matrix with a clean win in
either A/B.**

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

## Confirmation A/B of the shipped form (job 18127149, 2026-07-31)

The gap left by the ship gate — the shipped configuration is *six* levers, and was inferred
from a seven-lever run — is now closed on score.  `reach_ab6.R`: the same design, restricted
to the 11 `large`+`xlarge` members of `MBANK_FIXED_SAMPLE` (on small/medium the gate is inert,
so those cells carry no information), 5 seeds, 55 cells, both arms at the same raised
`targetHits`, levers as dots on a gate-free engine.

**Verdict: SHIP — confirmed.**  Paired final score **8 better, 1 worse, 46 tie**; no tier
regression.  Budget regime clean: 40 cells both-arms-at-deadline (equal-wall), 15 cells both
converged, **0 asymmetric**.

**The whole effect is 2 of the 11 matrices.**  Read this before quoting any tier number:

| matrix | tips | win | loss | tie |
|--------|------|-----|------|-----|
| project4284 | 4062 | **5** | 0 | 0 |
| project2771 | 94   | 3 | 1 | 1 |
| the other nine (63–173 t) | | 0 | 0 | 45 |

So the analyzer's `xlarge = 20/20 vs 15/20` is *one matrix*, exactly as in the ship gate.
Do not restate reach fractions as evidence: the union-best target is self-referential, so
base "misses" on 8 cells only because deltas6 got there.  **The paired counts are the
statistic.**

**`project2771` is noise, not a second winning matrix.**  Pooled across both A/Bs it is
**4 win / 2 loss / 4 tie** over 10 seeds — a coin flip, and it supplied the single loss in each
run.  Read the two-row table above as *one* winning matrix plus one high-variance matrix, and
treat any single 2771 cell as noise.

**The claim that survives both runs:** `project4284` is **10 win / 0 loss over 10 seeds**, and
**no matrix regresses net** in either A/B.  That is the whole of the positive evidence, and it
is enough — but it is one matrix, and it is the largest one in the battery.

**project4284 won with ZERO completed replicates.**  base completed 1–2 replicates in ~1330 s;
deltas6 completed **0** and still returned a tree 2–9 steps better on every seed.  At 4062
tips the escalation's benefit is *spending the whole budget deepening a single replicate*, not
more restarts.  Note this sits crosswise to the hard-tail study below, which concluded reach
is a **restart-volume** phenomenon — different regimes (482 vs 4062 tips), not a
contradiction, but the two mechanisms are opposite and neither generalises to the other.

**Cost.**  In this run `wall_total` ratio is 1.000 — a tautology of the shared deadline, not a
measurement; the honest figure is that at equal wall the escalation completes **2.24× fewer
replicates**.  On the 30 tied cells where both arms were deadline-bound (7 matrices, 86–173
tips) that cost nothing in final score.  The 16 tied cells that converged early (project4286
at 5–9 s of a 720 s cap, project4359 stopping on `targetHits` at 28 replicates) carry **no**
cost information and must not be counted as evidence of cost-neutrality.  The ×3.56 wall
figure belongs to the ship-gate run, which was not deadline-bound.

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
- **…and writing that lesson down did not stop me re-committing it.**  `reach_escalation_
  analyze.R` flagged truncation at `wall >= 0.95 × cap_s`, but cells stop at `0.90 × cap_s`,
  so it reported **1 of 110** deadline-bound ab6 cells when the true count was 59 — three
  sections below the paragraph above.  Both the harness (records `enum_time_fraction`) and the
  analyzer (derives the deadline from it, defaulting to 0.1) are fixed.  A threshold that
  encodes a magic number the engine owns will rot: derive it from recorded data.
- **Deadline-bound is not the same as spoiled.**  When *both* arms stop at the same deadline
  the cell is a valid equal-wall comparison — the stronger test, since the cheaper-per-replicate
  arm gets ~2.2× more replicates and still has to win.  What invalidates a cell is *asymmetry*
  (one arm converged, one cut off).  The analyzer now classifies cells three ways rather than
  discarding every deadline-bound cell, which would have thrown away 40 of 55 informative
  cells and printed "INCONCLUSIVE" over a clean result.
- **A tier win is not a tier property until you count matrices.**  Twice now an `xlarge`
  reach jump has been one matrix repeated across five seeds.  The analyzer prints the
  per-matrix win/loss/tie table and the "n of N matrices changed" line for this reason.
- **Union-best-across-arms targets make reach self-referential.**  If one arm alone attains a
  score the other misses *by construction*, so reach fractions inflate the loss count.
  Report paired win/loss/tie; use reach only as a secondary description.
- **You cannot `TreeLength()` a tree that `MaximizeParsimony()` returned by default.**
  `collapse = TRUE` (the default) contracts zero-length branches, and `TreeLength` errors with
  "`tree` must be binary".  To re-score a returned tree, ask for `collapse = FALSE`.  This
  killed the first tree-recovery run (job 18128376) after 22 minutes of search, in the
  *verification* step — so **persist the deliverable before verifying it**: write the trees
  out first, then check them, or a failed check throws away the search too.
- **A `collapse = FALSE` return is the pool verbatim.**  That makes it the natural place to
  test whether `poolSuboptimal` leaks suboptimal trees to the caller: re-score every returned
  tree and compare max against min.  Under a gate-free harness `escalatedPool` is FALSE, so
  the filter is inert and the leak (if real) is visible.
- **Whole-suite test runs need `test_local()`/`devtools::test()`**, not `test_dir()` +
  `library()`: several test files call internals unqualified and error otherwise.
- **A top-level `skip_on_cran()` makes a file report "0 pass 0 fail"** — that is *not* a
  pass, it never ran.  Use `NOT_CRAN=true`.
- **To attribute a suite failure, use a pristine detached worktree at the parent commit.**
  Swapping a single source file is invalid once you have added a test file, because your own
  tests then error on the missing functions and the totals stop being comparable.

## Not measured

Measured as of job 18127149: the **six levers' effect on score**, under `auto`→`thorough`/
`large`, equal weights, 61–4062 tips.

Still **unit-tested only** — no benchmark evidence, because that run passed the levers as dots
on a gate-free engine and therefore never executed the gate:

- the trigger threshold (`targetHits / defaultHits >= 2`);
- the `thorough`/`large` scoping;
- the `userSet` skip (caller-supplied fields must survive);
- the `escalatedPool` return filter.  `escalatedPool` was FALSE throughout both benchmark
  runs, so the filter that stops `poolSuboptimal = 3` leaking suboptimal trees into a
  `collapse = FALSE` result has never run outside `test-reach-escalation.R`.

Also not measured:

- **Equal weights only.**  Implied weights and profile parsimony are a deliberate but
  unmeasured extrapolation; ratchet depth there is governed separately by `.IwRatchetDepth`.
- **No tree was retained for the project4284 result** the confirmation rests on — the harness
  recorded `attr(res, "score")` only, so ~353 could not be independently re-scored at the
  time of writing.  A degenerate partial tree would score *worse*, so the win is very
  unlikely to be an artefact, but "unlikely" is what the record says.  Recovery run
  `reach_recover4284.R` re-runs the identical config, writes the trees and re-scores by label
  with `TreeLength` (job 18128376 died in the check — see the `collapse` trap above; job
  **18128526** is the corrected `collapse = FALSE` run).
- **The `escalatedPool` filter is only needed on the `collapse = FALSE` path** — verified by
  reading, not benchmark: the `collapse = TRUE` branch already restricts to
  `scores == best_score` before collapsing (`R/MaximizeParsimony.R`), which the `collapse`
  roxygen also documents.  So the default path never leaked; the fix sits on the one branch
  that did.
- **`poolReseed` (v2) cannot help where v1 helped most.**  It reseeds *replicates* from a
  pool needing `size >= 2`; project4284's winning arm completed zero replicates.  Any v2
  validation must therefore use matrices on which enough replicates complete for a pool to
  form — 4062 tips is the wrong test bed for it.
