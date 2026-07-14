# S7-LAND — lever #7: plain-EW scorer monomorphization (template-by-flavour)

**Date:** 2026-07-14 · **Base:** `cpp-search` `f5894837` · **Branch:** worktree only, NOT pushed
(supervisor merges) · **Spec:** `dev/plans/2026-07-14-lever7-scorer-monomorphization.md` ·
**Sizing:** `dev/profiling/s7-fastpath-sizing.md` · **Memory:** `tnt-per-move-kernel-gap`.

## Headline

Landed the EXACT (byte-identical) per-candidate speedup: the plain-EW SPR candidate scan is now a
**compile-time monomorphized template** (`spr_scan_plain_ew<UseFlat, Wrong>` in `src/ts_tbr.cpp`),
dispatched ONCE per weight-class at the scan site. The dead-in-plain-EW per-candidate branches
(NA / IW / sector_mask / constrained / collapsed / b2_ceiling) compile away, and — the actual payoff —
the `UseFlat=true` instantiation uses the weight-blind **flat kernel** (`fitch_indirect_cached_flat`),
which is UNSAFE as a global toggle but SAFE as an instantiation because dispatch selects it ONLY for
`all_weight_one` data.

- **Byte-identity gate: PASSED** (score + `candidates_evaluated` + per-pass `n_candidates_evaluated`),
  via BOTH `ts_tbr_diagnostics` AND full `MaximizeParsimony`, plus a gold-standard clean-HEAD
  cross-check and a positive control. **The MaximizeParsimony gate caught a real bug** (below).
- **Per-candidate SPR win (min-of-runs, direct):** flat kernel **+12–19%** of the SPR scan on
  unit-weight data; dead-branch strip alone **+6%** on weighted data. Matches the sizing prediction
  (~17% flat / ~5% strip).
- **Whole-search:** modest and dataset-dependent — Zhu2013 whole-call wall **+2.7%**; honest ceiling
  ~5% on unit-weight, ~1% on weighted, shrinking at 5432 scale (fixed cost / growing gather). **We take
  the 5%.** Orthogonal to 5432 reach — does NOT change which trees are found.

Default ON (it is exact); kill-switch `TS_EW_MONO_OFF` reverts to the general per-candidate-dispatch
loop.

## What changed (src/ts_tbr.cpp only; +155/−2)

1. `template<bool UseFlat, bool Wrong> spr_scan_plain_ew(...)` — file-scope, before `tbr_search`.
   Identity-skip + scorer + strict-`<` accept only. `if constexpr (UseFlat)` picks
   `fitch_indirect_cached_flat` (unit-weight) vs `fitch_indirect_length_cached` (general). `Wrong` is a
   positive-control-only corruption.
2. Dispatch at the SPR scan site: `if (ew_mono_plain) { use_flat ? <flat> : <cached> } else { <verbatim
   general loop> }`. `ew_mono_plain = ew_mono && !has_na && !use_iw && sector_mask==nullptr &&
   !constrained && collapsed_all_zero && !b2_ceiling`. The general loop is textually unchanged (baseline
   integrity — verified by git diff AND a clean-HEAD binary cross-check).
3. `collapsed_all_zero` gate + `refresh_collapsed_all_zero()` lambda, refreshed at EVERY
   `compute_collapsed_flags` recompute (see the bug below).
4. Positive-control machinery: `ew_mono` / `ew_mono_wrong` flags, split fire counters
   (`spr_mono_flat_fired` / `spr_mono_cached_fired`) surfaced in the `TS_IW_TIMING` line.

**Reroot path: no change, by design.** For unit-weight EW the reroot is ALREADY the branch-free x4
flat kernel (`fitch_indirect_cached_flat_x4`) — the sizing doc measured "nothing to strip" (control
delta ±2.5% noise). Its residual `has_na` branch is per-batch and well-predicted (~0%). Extracting the
250-line reroot loop into a template risks byte-identity for ~0% payoff, so it was NOT done; "cover both
paths" is satisfied honestly — the reroot is already monomorphic for the case that matters.

## Correctness gate (MANDATORY — passed before any timing)

Scripts (session scratchpad): `s7_gate.R` (ts_tbr_diagnostics + positive control),
`s7_gate_mp.R` (MaximizeParsimony), `s7_xcheck.R` / `s7_xcheck_mp.R` (clean-HEAD cross-checks),
`s7_determinism.R`. Datasets Wortley2006 / Zanol2014 / Zhu2013, gaps `-`→`?` (EW), ≥2 seeds.

- **`ts_tbr_diagnostics` (kernel-level), 6/6 PASS:** `score`, `n_evaluated`, and per-pass
  `n_candidates_evaluated` all exactly `identical()` between BASE (`TS_EW_MONO_OFF=1`) and NEW.
- **`MaximizeParsimony` (ratchet `use_flat=false` + sector regimes), 6/6 PASS:** `candidates_evaluated`,
  per-replicate scores, best score, and `n_topologies` all identical.
- **Gold-standard clean-HEAD cross-check:** a SEPARATE binary built from clean `f5894837` (`.agent-base`,
  no toggle) vs the templated binary, both at their DEFAULTS — byte-identical for BOTH
  `ts_tbr_diagnostics` and `MaximizeParsimony`. Rules out a symmetric perturbation of the general loop.
- **Positive control (replicated per spec, strengthened):** split fire counters prove the FLAT
  instantiation actually executed — `flat=450/450` (Wortley), `flat=803/803` (Zhu), `cached=1220/1220`
  (Zanol); i.e. Wortley2006 AND Zhu2013 are `all_weight_one` and drive `UseFlat=true`. A deliberately-
  wrong instantiation (`TS_EW_MONO_WRONG`, data-dependent always-≥1 corruption) DIVERGES on all three
  (flat and cached), proving the specialised path's output drives the search — not a no-op that passes
  the gate falsely (the false-0% trap from the sizing doc).
- **Determinism control:** `MaximizeParsimony` at `nThreads=1` + fixed seed is reproducible run-to-run
  (BASE==BASE, NEW==NEW) — so a BASE-vs-NEW difference is a real signal, not noise.
- **Targeted testthat (extra net, default-ON build):** 383 assertions across `test-ts-tbr-search`,
  `-symmetry`, `-collapsed`, `-tbr-dirty-rescore`, `MaximizeParsimony-features` — 0 failures.

### Bug found + fixed (why the MaximizeParsimony gate is mandatory)

The first MaximizeParsimony gate FAILED on Zhu2013/seed202: `candidates_evaluated` 262689699 (NEW) vs
262689659 (BASE) — +40 of 262M, identical scores/topologies. The determinism control proved it was
real, not noise. Root cause: `collapsed_all_zero` was computed ONCE at `tbr_search` entry, but
`compute_collapsed_flags(tree, ds, collapsed)` is **re-run after every accepted move** (and after each
root-edge re-descent) — a new tree can introduce zero-length branches that set collapsed flags. When it
did, the stale gate let the template skip the collapsed check while the general loop's
`if (use_collapsed && collapsed[below]) continue;` still fired → the template evaluated edges the
general loop skipped (NEW's +40, "evaluated more"). Fix: `refresh_collapsed_all_zero()` at every
recompute site. After the fix, all 6 MaximizeParsimony runs are byte-identical including Zhu2013/seed202
(262689659==262689659). `ts_tbr_diagnostics` never hit this (single descent, few/no zero-length
branches) — exactly the regime the spec warns the kernel oracle can't reach.

## Measurement (per-candidate ns, TS_IW_TIMING chrono, single-thread, min-of-runs, 20 reps)

Panel recoded `-`→`?` (EW). BASE = cached general loop; NEW = template. REROOT is the untouched control.
Script: `s7_timing.R` + `parse_timing2.R`.

| Dataset (n_tip) | weight class | BASE SPR ns | NEW SPR ns | SPR-specific | REROOT ctrl | SPR share of scan |
|---|---|---|---|---|---|---|
| Zhu2013 (75)    | all_weight_one → **flat**  | 14.35 | 12.08 | **+18.6%** | −2.7% | 20% |
| Wortley2006 (74)| all_weight_one → **flat**  | 14.80 | 12.94 | **+12.5%** | +0.1% | 37% |
| Zanol2014 (74)  | weighted → cached (strip)  | 19.35 | 18.10 | **+6.0%**  | +0.5% | 18% |
| Dikow2009 (88)  | weighted → cached (strip)  | 23.07 | 21.93 | **+5.8%**  | −0.9% | 21% |

- SPR-specific = SPR raw delta − REROOT control delta. The flat kernel roughly triples the strip-only
  win on unit-weight data (+18.6% vs +6.0%), exactly the sizing prediction (~17% flat vs ~5% strip).
- **Whole-call wall-clock** (a full `ts_tbr_diagnostics` descent, min of 12×40 loops, `s7_wall.R`):
  Zhu2013 **+2.7%**; Wortley2006 / Zanol2014 within timer noise (both fast; the weighted strip on an
  18%-share SPR is ~1%, sub-resolution). The whole call dilutes the SPR win with non-scan overhead
  (full_rescore, per-clip precompute), so ~2–4% whole-search on unit-weight, ~1% weighted — consistent
  with the spec's ~5% ceiling.

## Honest framing

- The lever is EXACT and byte-identical (gate passed on both harnesses + a clean-HEAD cross-check).
- The win is a modest, general wall-clock lever: biggest where the flat kernel fires (all_weight_one)
  AND the SPR share is high (Wortley 37%). On weighted data only the ~5%-of-SPR dead-branch strip
  applies. Per the sizing doc the fixed branch/bookkeeping saving is a SMALLER fraction at 5432 scale
  (gather grows with n_tips), so the panel is the conservative upper bound.
- **Orthogonal to 5432 reach** — byte-identical means it does NOT change which trees are found; do NOT
  sell it as a 5432 fix. "We'll take the 5% when we find it."

## Build / provenance

- Full `R CMD INSTALL --preclean` under `CCACHE_DISABLE=1` for the first build; `dev/build-fast.R`
  (ccache incremental — safe, no header touched) for iteration. Per-agent libs in the session
  scratchpad (`.agent-s7` = templated, `.agent-base` = clean HEAD, `.agent-assert` = flat==cached
  assertion build).
- The `-DTS_EW_MONO_ASSERT` diagnostic (flat vs cached accept-decision check) confirmed ZERO mismatches
  across a full MaximizeParsimony run, empirically proving flat≡cached on all_weight_one; the temporary
  assertion was then removed from the deliverable. `-D` went in `PKG_CPPFLAGS` (the machine's
  `~/.R/Makevars.win` zeroes `PKG_CXXFLAGS`), confirmed live by grepping the compile line.
- Deliverable: `src/ts_tbr.cpp` (+155/−2) + this doc. Not pushed.
