# S7-LAND — lever #7: plain-EW scorer monomorphization (template-by-flavour)

**Date:** 2026-07-14 · **Base:** `cpp-search` `f5894837` · **Branch:** worktree only, NOT pushed
(supervisor merges) · **Spec:** `dev/plans/2026-07-14-lever7-scorer-monomorphization.md` ·
**Sizing:** `dev/profiling/s7-fastpath-sizing.md` · **Memory:** `tnt-per-move-kernel-gap`.

## Headline

Landed the EXACT (byte-identical) per-candidate speedup: the plain SPR candidate scan is now a
**compile-time monomorphized template**, dispatched ONCE per criterion flavour at the scan site —
`spr_scan_plain_ew<UseFlat, Wrong>` for equal-weights and `spr_scan_plain_iw<Wrong>` for implied
weights (both in `src/ts_tbr.cpp`). The dead-in-that-regime per-candidate branches
(NA / the other criterion / sector_mask / constrained / collapsed / b2_ceiling) compile away, and — the
EW payoff — the `UseFlat=true` instantiation uses the weight-blind **flat kernel**
(`fitch_indirect_cached_flat`), which is UNSAFE as a global toggle but SAFE as an instantiation because
dispatch selects it ONLY for `all_weight_one` data. IW has no flat analogue (per-pattern weights), so
its win is the dead-branch strip only.

- **Byte-identity gate: PASSED for BOTH EW and IW** (score + `candidates_evaluated` + per-pass
  `n_candidates_evaluated`), via `ts_tbr_diagnostics` AND full `MaximizeParsimony`, plus a gold-standard
  clean-HEAD cross-check and a positive control. **The MaximizeParsimony gate caught a real bug** (below).
- **Per-candidate SPR win (min-of-runs, direct):** EW flat kernel **+12–19%** of the SPR scan on
  unit-weight data; EW dead-branch strip **+6%** on weighted; **IW strip +2.7–4.0%** (heavier
  double-accumulator scorer → smaller branch fraction, no flat bonus). Matches the sizing prediction.
- **Whole-search:** modest and dataset-dependent — the one resolvable unit-weight run,
  full-`MaximizeParsimony` on Zhu2013, is **+1.3%** (single-descent proxy +2.7%); weighted Zanol2014
  +0.1% (noise). The ~5% is a scan-level ceiling, NOT a whole-search floor, and shrinks at 5432 scale
  (fixed cost / growing gather). **We take the 5% where we find it.** Orthogonal to 5432 reach — does
  NOT change which trees are found.

Default ON (it is exact); kill-switch `TS_EW_MONO_OFF` reverts to the general per-candidate-dispatch
loop.

## What changed (src/ts_tbr.cpp only; +232/−2)

1. `template<bool UseFlat, bool Wrong> spr_scan_plain_ew(...)` — file-scope, before `tbr_search`.
   Identity-skip + scorer + strict-`<` accept only. `if constexpr (UseFlat)` picks
   `fitch_indirect_cached_flat` (unit-weight) vs `fitch_indirect_length_cached` (general). `Wrong` is a
   positive-control-only corruption.
2. `template<bool Wrong> spr_scan_plain_iw(...)` — the IW/XPIWE analogue. Same dead-branch strip, but
   NO flat variant: implied weights must apply the per-pattern `iw_delta`, so there is no weight-blind
   kernel to swap in — the win is the ~5%-of-SPR strip only. Scorer + accept reproduced verbatim from
   the general loop's IW branch (incl. the dead-on-IW `cutoff` update, kept so byte-identity does not
   rest on a "cutoff is IW-dead" argument).
3. 3-way dispatch at the SPR scan site: `if (ew_mono_plain) {flat|cached} else if (iw_mono_plain) {iw}
   else { <verbatim general loop> }`, gated on a shared `mono_plain = ew_mono && !has_na &&
   sector_mask==nullptr && !constrained && collapsed_all_zero && !b2_ceiling` (then `!use_iw` for EW,
   `use_iw` for IW). The general loop is textually unchanged (baseline integrity — verified by git diff
   AND a clean-HEAD binary cross-check).
4. `collapsed_all_zero` gate + `refresh_collapsed_all_zero()` lambda, refreshed at EVERY
   `compute_collapsed_flags` recompute (see the bug below).
5. Positive-control machinery: `ew_mono` flag (default ON; `TS_EW_MONO_OFF` reverts BOTH EW+IW to the
   general loop), the compile-gated `ew_mono_wrong`, and split fire counters
   (`spr_mono_flat_fired` / `spr_mono_cached_fired` / `spr_mono_iw_fired`) in the `TS_IW_TIMING` line.

**Reroot path: no change, by design (EW and IW).** The reroot is ALREADY x4-batched and branch-free for
the case that matters: `fitch_indirect_cached_flat_x4` for unit-weight EW, `indirect_iw_cached_flat_x4`
for IW/XPIWE. The sizing doc measured "nothing to strip" on the EW reroot (control delta ±2.5% noise);
the reroot control deltas in the IW timing table above (±0.4%, and the −3.0% Dikow outlier) confirm the
same for IW. The residual per-batch `has_na` branch is well-predicted (~0%). Extracting the ~250-line
reroot loop into a template risks byte-identity for ~0% payoff, so it was NOT done; "cover both paths"
is satisfied honestly — the reroot is already monomorphic where it counts, for both criteria.

## Correctness gate (MANDATORY — passed before any timing)

Scripts (session scratchpad): `s7_gate.R` (EW ts_tbr_diagnostics + positive control),
`s7_gate_mp.R` (EW MaximizeParsimony), `s7_gate_iw.R` / `s7_gate_mp_iw.R` (IW analogues, concavity=10),
`s7_xcheck.R` / `s7_xcheck_mp.R` (clean-HEAD cross-checks), `s7_determinism.R`. Datasets
Wortley2006 / Zanol2014 / Zhu2013, gaps `-`→`?`, ≥2 seeds (EW) / 3 seeds (IW MP).

- **`ts_tbr_diagnostics` (kernel-level), 6/6 PASS:** `score`, `n_evaluated`, and per-pass
  `n_candidates_evaluated` all exactly `identical()` between BASE (`TS_EW_MONO_OFF=1`) and NEW.
- **`MaximizeParsimony` (ratchet `use_flat=false` + sector regimes), 18/18 PASS (6 seeds × 3 sets):**
  `candidates_evaluated`, per-replicate scores, best score, and `n_topologies` all identical. Extended
  from 2 to 6 seeds specifically to confirm the collapsed fix is not seed-luck (the bug below first
  surfaced on 1 of 2 seeds).
- **Gold-standard clean-HEAD cross-check:** a SEPARATE binary built from clean `f5894837` (`.agent-base`,
  no toggle) vs the templated binary, both at their DEFAULTS — byte-identical for BOTH
  `ts_tbr_diagnostics` and `MaximizeParsimony`. Rules out a symmetric perturbation of the general loop.
- **IW/XPIWE gate (`s7_gate_iw.R`, concavity=10), byte-identical:** `ts_tbr_diagnostics` **6/6** (IW
  scores are doubles — identical to full precision — plus `n_evaluated` + per-pass) and
  `MaximizeParsimony` (XPIWE) **9/9** (3 seeds × 3 sets: `candidates_evaluated`, per-replicate scores,
  best, `n_topologies` all identical). The IW template fired 100% of clips (`iw=764/1563/1173`,
  `flat=0 cached=0`, `regime=IW`), and the compile-gated wrong IW instantiation DIVERGES on the flagged
  build (Wortley 30.04→29.84, Zanol 64.54→64.92, Zhu 47.22→47.72). **EW regression re-checked 6/6 PASS**
  after the 3-way dispatch refactor.
- **Positive control (replicated per spec, strengthened):** split fire counters prove the FLAT
  instantiation actually executed — `flat=450/450` (Wortley), `flat=803/803` (Zhu), `cached=1220/1220`
  (Zanol); i.e. Wortley2006 AND Zhu2013 are `all_weight_one` and drive `UseFlat=true`. A deliberately-
  wrong instantiation (data-dependent always-≥1 corruption) DIVERGES on all three (flat and cached),
  proving the specialised path's output drives the search — not a no-op that passes the gate falsely
  (the false-0% trap from the sizing doc). The wrong instantiations are emitted ONLY under a compile
  flag (`-DTS_EW_MONO_WRONG`), so a PRODUCTION binary has no env-var path to a corrupted scorer:
  verified that `TS_EW_MONO_WRONG=1` is INERT on the default build, and the divergence reproduces only
  on a flag-built binary (`.agent-wrong`).
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
recompute site. After the fix, all 18 MaximizeParsimony runs (6 seeds) are byte-identical including
Zhu2013/seed202 (262689659==262689659). `ts_tbr_diagnostics` never hit this (single descent, few/no
zero-length branches) — exactly the regime the spec warns the kernel oracle can't reach.

The fix is structurally complete, not a patch of one instance: `collapsed_all_zero` is the ONLY
per-clip-mutable input to the `ew_mono_plain` gate. Every other input (`ew_mono`, `has_na`, `use_iw`,
`sector_mask`, `constrained`, `b2_ceiling`, and the `use_flat` weight-class sub-branch) is computed once
and constant for the whole `tbr_search` call — the reroot path reads the SAME once-computed `use_flat`,
so SPR and reroot cannot diverge on weight class. `collapsed` was uniquely dangerous because it is a
vector recomputed mid-search; refreshing the bool at every recompute closes the entire staleness class.

## Measurement (per-candidate ns, TS_IW_TIMING chrono, single-thread, min-of-runs, 20 reps)

Panel recoded `-`→`?` (EW). BASE = cached general loop; NEW = template. REROOT is the untouched control.
Script: `s7_timing.R` + `parse_timing2.R`.

| Dataset (n_tip) | weight class | BASE SPR ns | NEW SPR ns | SPR-specific | REROOT ctrl | SPR share of scan |
|---|---|---|---|---|---|---|
| Zhu2013 (75)    | all_weight_one → **flat**  | 14.35 | 12.08 | **+18.6%** | −2.7% | 20% |
| Wortley2006 (74)| all_weight_one → **flat**  | 14.80 | 12.94 | **+12.5%** | +0.1% | 37% |
| Zanol2014 (74)  | weighted → cached (strip)  | 19.35 | 18.10 | **+6.0%**  | +0.5% | 18% |
| Dikow2009 (88)  | weighted → cached (strip)  | 23.07 | 21.93 | **+5.8%**  | −0.9% | 21% |

IW (concavity=10, `s7_timing_iw.R` — all datasets go through `spr_scan_plain_iw`; strip-only, no flat):

| Dataset (n_tip) | BASE SPR ns | NEW SPR ns | SPR-specific | REROOT ctrl | SPR share of scan |
|---|---|---|---|---|---|
| Zhu2013 (75)    | 38.38 | 36.87 | **+4.0%** | −0.0% | 15% |
| Wortley2006 (74)| 36.75 | 35.30 | **+3.5%** | +0.4% | 35% |
| Zanol2014 (74)  | 50.57 | 49.19 | **+2.7%** | +0.1% | 21% |
| Dikow2009 (88)  | 54.27 | 54.90 | +1.8% (within −3.0% control noise) | −3.0% | 24% |

- SPR-specific = SPR raw delta − REROOT control delta. For **EW** the two weight classes match the sizing
  prediction: the dead-branch strip alone recovers ~5% of SPR (weighted: Zanol +6.0%, Dikow +5.8%),
  and adding the flat kernel roughly triples that to ~13–19% (unit-weight: Zhu +18.6%, Wortley +12.5%).
  The clean WITHIN-dataset flat-vs-strip comparison is the sizing doc's Zhu2013 figures (−6.5% strip →
  −17.4% flat); the unit-vs-weighted rows above differ by dataset too, so read them as "flat-eligible
  vs strip-only," not a controlled flat/strip contrast.
- **IW** SPR-specific is +2.7–4.0% (Dikow +1.8% is within its own control noise) — the strip only, as
  predicted (no flat kernel possible). It is SMALLER than EW's ~5–6% strip because the IW
  double-accumulator scorer is ~2.5–4× heavier per candidate (BASE SPR 38–54 ns vs EW's 14–23 ns), so
  the FIXED branch cost is a smaller fraction of it. IW/XPIWE is the production default, so this is a
  real (if small, exact) win on the common path; whole-scan(timed) ~0.6–1.7%.
- **Single-descent whole-call** (a full `ts_tbr_diagnostics` descent, min of 12×40 loops, `s7_wall.R`):
  Zhu2013 **+2.7%** — an UPPER proxy (one plain-EW TBR descent, the regime where the flat kernel fires
  most).
- **Faithful whole-search** (a full `MaximizeParsimony` run — ratchet + sector + drift + fuse +
  enumeration — min-of-runs, `nThreads=1`, `s7_wall_mp.R`): **Zhu2013 +1.3%** (18.8s → 18.5s). Lower
  than the single-descent proxy, as expected: the ratchet (upweight → cached, not flat), sectorial, and
  enumeration phases don't hit the flat SPR path, diluting the win. Zanol2014 (weighted, cached path
  only): **+0.1%** (22.83s → 22.80s, within run-to-run noise, as expected for a strip-only 18%-share
  SPR). Wortley2006's full-MP converges in ~0.27s (too short to resolve above timer noise). So the one
  meaningfully-resolvable whole-search number is **+1.3% (Zhu, unit-weight)**, with weighted sub-1% —
  a modest, EXACT lever; the ~5% is a scan-level ceiling, not a whole-search floor, and it shrinks
  further at 5432 scale (fixed cost / growing gather).

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
