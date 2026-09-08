# TNT-style sector-size scaling for RSS (cpp-search)

Status: **implemented, flag-gated, default-OFF; awaiting matched-wall Hamilton A/B.**
Closes the human-raised thread *"TNT sets sector size based on nTip (nTip/2?) once
nTip exceeds a threshold (~88?)"* for `large` recipes. From the TNT feature-gap
audit; a candidate reach lever for large trees.

## What TNT has (`help sectsch`)

- `minsize N` / `maxsize N` — bounds for random sector selections.
- `increase N` — grow size when enough selections of the current size complete:
  `S := S + (S*N)/100`.
- `selfact N` — max selections of size `S` for `T` active taxa: `M = (T*100)/(N*S)`.
- `moveon N` — if `N` selections fail to improve, move on.
- `xmult: level 0-10` sets intensity (incl. sector size) scaled by taxon count.

Today TreeSearch's RSS uses a **fixed** `sectorMaxSize` per preset (50 default /
80 thorough & large). There is no nTip-scaled cap and no adaptive growth. (The
long-rumoured hard "45 cap" does not exist — it is the tunable `sectorMaxSize`.)

## Two additive levers (both in `rss_search`, `src/ts_sector.cpp`)

### (a) nTip-scaled cap
Effective max sector size = `min(cap, ceil(frac * nTip))` once `nTip >= threshold`,
where `cap` is the absolute ceiling. **Design note (critical):** with a small
preset ceiling (50/80) the `min()` only ever *shrinks* sectors or is inert on
large trees — so the ceiling must be raised for the fraction to bind. `TS_SECT_MAXSIZE`
raises it. (`min` is deliberate and TNT-faithful: `maxsize` is a hard cap; the
taxon-scaled size ramps *up toward* it, never past it, so a 3000-tip tree can't
get 1500-tip "sectors".)

### (b) adaptive growth (increase / selfact / moveon)
Ramp the sector-size window's upper bound `S` from a small start up to the
effective cap: do `M = (T*100)/(selfact*S)` selections at each `S` (or the
existing `2T/S` auto-heuristic when `selfact` unset), grow `S` by `increase`
percent per band, and stop early after `moveon` consecutive non-improving
selections **once the ramp has reached the cap** (`S == eff_max`). Gating
`moveon` on the ceiling — and resetting the fail counter on each growth — is
deliberate: an ungated `moveon` would trip inside the first band and leave the
growth lever inert, so the A/B would silently measure "no growth". Selection is
uniform within `[minSize, S]` (conflict-weighting is disabled under growth to
avoid desync on the per-band rebuild).

## Knob reference

| Concern | `SectorParams` field | Env override | Default (OFF) |
|---|---|---|---|
| absolute ceiling | `max_sector_size` | `TS_SECT_MAXSIZE` | preset (50/80) |
| nTip fraction | `max_sector_size_frac` | `TS_SECT_MAXFRAC` | `0.0` (off) |
| fraction threshold | `sector_size_threshold` | `TS_SECT_THRESHOLD` | `88` |
| growth on/off | `sector_grow_increase>0` | `TS_SECT_GROW` (`0`=off, else on) | off |
| increase % | `sector_grow_increase` | `TS_SECT_GROW_INC` | `0` |
| selfact | `sector_grow_selfact` | `TS_SECT_GROW_SELFACT` | `0` (auto `2T/S`) |
| moveon | `sector_grow_moveon` | `TS_SECT_GROW_MOVEON` | `0` (no early stop) |
| start size | `sector_grow_start` | `TS_SECT_GROW_START` | `0` (→ minSize) |

Env wins over params, so a matched-wall A/B runs from a **single isolated-lib
build** with no R change. `SectorParams` fields are declared for a follow-up
R-wiring PR (SearchControl + `large` preset default), which lands only if the
A/B wins.

## Off-path is byte-identical to baseline (proven)

Every off-path change is either (i) a variable that *equals* the old constant
when off (`eff_max == max_sector_size`, `size_upper == eff_max`, `loop_cap ==
n_picks`), (ii) a pure read (`best_before`), or (iii) a block guarded by
`if (grow_enabled)`. So with no env override and default params the executed
statements and RNG-draw sequence are unchanged.

Verified empirically: a fixed-seed `ts_rss_search` run against a pristine-HEAD
build and against this build gave identical `score`, edge checksum, and sector
counts. End-to-end via `MaximizeParsimony` (through `ts_driven`), the same-seed
baseline score is identical before and after the env block, and the env arm
diverges (engagement). Unit tests: `tests/testthat/test-ts-sector-scaling.R`.

## Validation plan (Hamilton only — no local heavy compute)

1. **Isolated-lib build** (skills/hamilton/r-infrastructure.md). Env
   `TERM=dumb OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Build for trajectory
   validity with `CCACHE_DISABLE=1 --preclean` (stale-object ABI gotcha).
2. **Engagement diff** first: same-seed trajectory diff on one large matrix to
   prove the knob engages. Feature (b) diverges on its own; feature (a) needs
   the raised ceiling (`TS_SECT_MAXSIZE`) to diverge — without it the `min()`
   clamps to the preset cap and the A/B would be a false negative.
3. **Matched-wall A/B** on large matrices (nTip >= 120) and the uncracked hard
   floors, TRAINING split only, `rasStarts >= 2`. Anytime-curve / time-to-optimum
   reporting. Harness: `dev/benchmarks/bench_sector_scaling_ab.R` (arms:
   `baseline` / `capraise` / `frac` / `grow` / `both`).

Adopt (preset-default-ON for `large`) **only** on a matched-wall win with no
regression; the knob values in the harness are the tuning surface and must be
selected on the training split. The validation / project175 split stays
sequestered.

**Interpretation caveat.** These levers touch **RSS only** (the correct and only
home for TNT `minsize/maxsize/increase/selfact/moveon`), yet the `large` recipe
also runs `xssRounds=5 / rssRounds=3 / cssRounds=2`. If XSS/CSS dominate the
sector gain on a given matrix, a null A/B may reflect RSS's small share of total
sector effort rather than the lever failing. When reading a null result, check
the RSS phase's contribution (verbosity timings) before concluding the scaling
lever is inert; consider an RSS-weighted recipe (raise `rssRounds`, lower
`xssRounds`) as a follow-up arm.
