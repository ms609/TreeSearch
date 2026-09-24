# Route-3: does cheap exact-drift earn a place in the recipe? — NO-GO

Follow-on to `drift-exactness-gate.md`. The drift matched-wall gate showed
exact-drift wins the `auto`/wall-race regime; route-3 asked whether that's worth
**re-tuning recipes** to install/expand drift in wall-limited presets. Answer:
**no** — the marginal-value win does not survive into the full search ensemble.

## Experiments (all Hamilton, EW, reuse lib-driftexact)

1. **Marginal value** (`drift-exactness-route3-bench.R` + `-replicate-bench.R`):
   from a shared per-seed local optimum, best marginal spend of a wall-second —
   arms `ratchet` / `drift_exact` / `drift_union` / `replicate`. 3 gate datasets.
2. **Recipe whole-search** (`drift-exactness-recipe-bench.R`): whole `thorough`
   search on **8 training-split MorphoBank keys** (65–135t, EW-recoded, sector
   drift OFF to isolate main-tree drift), arms `none` / `union` / `exact`,
   `maxReplicates` knob. TRAINING split only (sequestered; harness refuses
   non-training keys). Raw CSVs in `route3-data/`; reanalyse `route3-analyze.R`.

## Result 1 — marginal value: exact-drift is the BEST single continuation
At matched wall on all 3 datasets, `drift_exact` beats `ratchet`, `drift_union`
AND `replicate` (a fresh restart from a local-opt is the *worst* — it discards
progress). E.g. Dikow2009 @0.44s: exact 1607.85 < ratchet 1608.30 < union 1608.20
< replicate 1610.00. So *in isolation*, a wall-second is best spent on exact-drift.

## Result 2 — recipe whole-search: drift is REDUNDANT in the ensemble (both ends)
The isolated win does **not** translate. Across the 8 training keys (4 were
trivially converged = zero signal; effective N=4):

- **Saturated (reps=8):** 4/8 exact-tie; non-ties sub-step noise split across
  arms; **`none` ties on score and is faster** (project2124: none 33.2s vs union
  41.9s vs exact 35.6s). Adding drift costs ~10–25% wall for ~nil.
- **Tight budget (the regime route-3 is *about*; reps=2 / low W, non-flat keys):**
  still no systematic exact win — project3617 `none` best, project4614 `union`
  best (exact ≈ none), project2124 `exact` best, project3807 wash. exact-drift
  beats no-drift on **1/4** keys; where drift helps it is as often union as exact.

So the multistart + ratchet + sectorial + TBR ensemble already captures the score
that exact-drift finds as an isolated continuation; layering main-tree drift on
top is redundant at every budget and only adds wall. (Consistent with
`diversity-generation-gates`: main-tree drift OFF on the hard path.)

## Verdict & what lands
- **NO recipe surgery.** Do not install/expand main-tree drift in any preset;
  per `completeness-secondary-to-wallclock`, adding wall for no reach gain is a
  disqualifier. exact-drift stays **opt-in** (`TS_DRIFT_EXACT`).
- **exact-drift's only residual value** is as the *cheaper* scorer *where drift
  already runs* (thorough's sector godrift): ~equal reach, less wall — a minor
  opt-in, not a default flip. Sector-godrift scorer itself was a wash/crossover
  (`drift-exactness-gate.md`), reinforcing opt-in.
- **The genuine win from this line of work is `spr_search`** — exact promoted to
  default (a real premature-convergence fix), see `drift-exactness-gate.md`.
- **D preset-refactor: not worth it.** No preset default should flip to
  exact-drift, so the ABI/`SearchControl` plumbing has no consumer; the existing
  `TS_DRIFT_EXACT` env flag suffices as the opt-in.
- **mbank 180t marginal-value curve: cancelled** (both original tasks TIMED OUT
  at the 8h cap; the resubmit was cancelled once recipe-B settled the decision) —
  it would only have re-confirmed the isolated marginal win, which is moot given
  the ensemble redundancy.
