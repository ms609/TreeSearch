# Sectorial profiling — production-ready levers (Round 6, 2026-06-20)

This separates the **shippable, byte-identical** wins from the env-gated
measurement scaffolding, per the /profile round on the isolated sectorial
component (findings.md T-S6a–e).

## What to land: `sector-levers.patch`

`sector-levers.patch` applies **4 byte-identical micro-levers to `src/ts_sector.cpp`
only**, against pristine cpp-search (`da0f203f`):

1. **`compute_from_above_for_sector`** — hoist the per-step `new_from_above`
   allocation out of the path loop (allocate once + `std::swap`).
2. **`search_sector` `ras_starts==1` fast path** (the default) — skip the
   provable no-op `best_*` snapshot + post-loop restore + `build_postorder`
   round-trip (the single start's result already sits in `rd.subtree`, and
   `reinsert_sector` never reads postorder).
3. **`search_sector` getenv hoist** — `TS_FREE_HTU_PROBE` was read 2–3× per
   sector pick; cache it in one `static const bool`.
4. **`rss_search` getenv hoist** — same for the per-accept `TS_SECT_DEBUG`.

Apply with: `git apply --directory=src dev/profiling/sector-levers.patch`
(the patch paths are file-local; adjust `-p`/`--directory` to your layout).

### Evidence (byte-identical + faster)
- **Byte-identical:** the CLEAN patched `ts_sector.cpp` (no measurement code)
  produces identical per-call score + n_sectors_searched/improved vs pristine
  across {Zanol,Zhu,Wortley} × {Wagner,TBR} starts (6/6). The
  instrumentation-included build additionally passed 12/12 across 2 seeds and
  8 search test files.
- **Final-build score-identity end to end:** the mission A/B (full
  `MaximizeParsimony`, `dev/profiling/drivers/mission-getenv-ab.R`) returned
  identical scores on 4 datasets × 3 seeds (625/624/625, 1261×3, 479×3, 272×3).
- **Wall:** ~2.8 % of isolated-sectorial wall (Zanol 48×80 base 3.08→2.99 s,
  8/8 rounds faster); breakdown ras_starts fast path ~0.05 s, getenv hoist
  ~0.028 s, from_above swap ~0.006 s — sums exactly to the A/B delta.

## NOT in this patch (deliberately)

- **`src/ts_tbr.cpp` is entirely excluded.** My ts_tbr.cpp changes are
  (a) env-gated measurement timing and (b) the per-clip `TS_REVERT_CHECK` /
  `TS_IW_SCANCHK` / `TS_PHYS_REROOT` getenv hoists — the latter are the
  **TBR-agent's production domain** and worth ~20–26 % MISSION-WIDE wall on
  their own (findings T-S6d, memory `getenv-ucrt-cost`, spawn_task
  `task_2f451c4f`). Land those via a clean compile-out / centralized
  debug-flags mechanism on cpp-search, NOT by copying my instrumentation.
- **Measurement-only scaffolding** in the working-tree `ts_sector.cpp`
  (`TS_SECT_TIMING` chrono timers + counters, `TS_SECT_NOREROOT`,
  `TS_SKIP_RSS_GTBR`, Probe-A counters) — env-gated, zero behaviour change when
  off, but should be stripped before merge. The patch above already excludes
  all of it.

## Efficiency axis (T-S6e) — no clean behaviour-neutral win
Probed and recorded; the top idea (suppress the redundant trailing global TBR)
is a speed/quality tradeoff, not safe. See findings.md T-S6e.
