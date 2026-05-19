# Profiling rounds — log

One entry per `/profile` invocation. Most recent at the top.

Append a new round when you finish step 6 of the round. Update `last_focus:`
at the bottom of the file before saving.

## Round template

```
### Round N — YYYY-MM-DD — area #K (<name>)
- Driver:        dev/profiling/drivers/<area>.R   (bare wall: X.X s)
- Build:         .vtune-lib mtime YYYY-MM-DD HH:MM (vs src/ HH:MM)
- profvis:       <2 % R overhead | top R line / [Port] finding>
- VTune top 3:   <fn1 X %>, <fn2 Y %>, <fn3 Z %>   (module=TreeSearch.dll)
- Finding:       [Port|Optimise|AT-LIMIT] short — verified Δ via micro-bench
- Filed:         T-NNN row(s) in findings.md
- Cleanup:       result_<area>_<date> removed; .vtune-lib <kept|deleted>
- Next reviewer: <what to look at next time on this area>
```

---

### Round 0 — 2026-05-18 — scaffold
- Scaffolded `dev/profiling/` per `/profile init`.
- Built focus-areas.md from the phase distribution in
  `.positai/expertise/profiling.md` (Zhu2013 thorough, 2026-03-27) and the
  active profiling tasks in `to-do.md` (T-274, T-298, T-300, S-PROF round 7).
- Three live targets ranked first: NNI-perturb (#1), Ratchet (#2),
  RSS/sector (#3).
- AT-LIMIT recorded for Wagner (#10), per-candidate indirect scoring (#11),
  R-loop search (#12) — rotation will skip these.
- No profiling run on this turn (per skill: do not profile on same turn as
  scaffold; user should review the ranking first).

---

### Round 1 — 2026-05-18 — area #2 (Ratchet inner loop)
- Driver:        dev/profiling/drivers/ratchet.R   (bare wall: 2.80 s median; Zhu2013 thorough ×1 rep)
- Build:         .vtune-lib rebuilt 2026-05-18 from TreeSearch_2.0.0.tar.gz (CXXFLAGS=-O2 -g -fno-omit-frame-pointer via dev/profiling/Makevars.vtune; R CMD INSTALL into .vtune-lib)
- profvis:       ~3 % R overhead (`MaximizeParsimony` wrapper); `ts_driven_search` dominates → no [Port] finding
- VTune top 3:   NOT COLLECTED — VTune not installed on this machine; WPR (`wpr -start CPU`) requires admin. Phase-level data from verbosity=2 used instead: Ratchet 62 % of inner-loop search time (802 ms / 1 301 ms); XSS 3 %, RSS 11 %, CSS 9 %, Wagner+TBR 15 %
- Finding:       [Profiled — unverified] Ratchet is a TBR wrapper. Perturbation save/restore overhead is O(n_chars) ≪ TBR time. T-300 (`full_rescore` after every accepted TBR move, `ts_tbr.cpp:1136–1137`) is the actionable target; implementation plan in `.AGENTS/memory/t300-lazy-tbr-rescore.md`. Cannot verify Δ without per-function hotspot data.
- Filed:         No new row in findings.md (unverified). T-300 already tracked. Note: also fixed `flat_blocks`/`all_weight_one` missing from `build_reduced_dataset` in ts_sector.cpp (S-RED area-5 finding, done inline).
- Cleanup:       No VTune result dirs; tarball removed (TreeSearch_2.0.0.tar.gz); .vtune-lib kept (needed for next round)
- Next reviewer: Install VTune (or run gprof build with `-pg`) to measure `full_rescore` share inside `tbr_search` → then implement T-300 → re-profile to verify speedup. Baseline: 2.80 s/rep (see baselines.md).

---

### Round 2 — 2026-05-19 — area #4 (TBR full-rescore at acceptance)
- Driver:        dev/profiling/drivers/tbr-rescore.R   (bare wall: 3.9 s; 12 ratchet reps × nCycles=12; Zhu2013 75t nThreads=1)
- Build:         dev/profiling/.vtune-lib-20260519061049 (built 2026-05-19 06:10:49 from HEAD c504ea87, src/ clean; CXXFLAGS=-O2 -g -fno-omit-frame-pointer; debug symbols confirmed via objdump .debug_info/.debug_line)
- profvis:       <2 % R overhead (258/287 samples in ts_ratchet_search; remaining = loadNamespace startup, not hot path)
- VTune top 5 (TreeSearch.dll, 3.211 s total DLL CPU):
    1. ts::fitch_na_score       0.585 s  18.2 %  (full-tree Fitch pass — full_rescore path confirmed via callstack)
    2. ts::simd::any_hit_reduce_avx2  0.309 s  9.6 %  (SIMD candidate hit reduction, inner evaluation)
    3. ts::tbr_search (residual)  0.297 s  9.3 %  (control-flow overhead outside child callees)
    4. ts::fitch_na_pass3_score  0.281 s  8.8 %  (incremental uppass, candidate evaluation)
    5. ts::fitch_na_incremental_uppass  0.110 s  3.4 %
- full_rescore attribution:
    - ts::fitch_na_score (self) + load_tip_states = 0.617 s = **19.2 % of DLL time (self-time lower bound)**
    - Attribution method: callstack report confirms fitch_na_score → fitch_score_ew → full_rescore → tbr_search → ratchet_search
    - ⚠ Caveat: 19.2 % is self-time only. fitch_na_score has SIMD callees (any_hit_reduce_avx2 0.309s, 9.6 %) that are shared with the incremental evaluation path — unknown fraction comes from full_rescore vs incremental. [Unknown source file] 2.076 s (39 %) includes inlined code from both paths. Full_rescore **inclusive time** is plausibly 22–30 %. The prior S-PROF round 7 estimate of 28 % was likely inclusive time and is not contradicted by the 19.2 % self-time measurement — they measure different things.
    - full_rescore at line 1138 (acceptance) >> line 563 (entry): ratchet-driven TBR accepts ~100–200 moves per call from perturbed trees vs 1 entry call, so ~99% of full_rescore time is the acceptance-path T-300 target
    - Source-line attribution for lines 1138/1283 not available via software sampling (inlined into [Unknown]).
- Finding:       [Optimise] T-300 is confirmed: full_rescore after accepted move ≥ 19.2 % of DLL time (inclusive estimate 22–30 %). Incremental path (fitch_na_pass3_score + incr_uppass + incr_downpass = 12.2 % self) already costs less per call. T-300 (in-flight by parallel agent) is justified — predicted gain 15–30 % of DLL time.
- Filed:         T-300 row in findings.md (unverified — micro-bench pending T-300 implementation)
- Cleanup:       result_tbr-rescore_20260519/ removed; .vtune-lib-20260519061049 deleted
- Next reviewer: After T-300 lands — re-run this driver to verify fitch_na_score drops from 18.2 % toward the incremental path baseline. Also look at ts::simd::any_hit_reduce_avx2 (9.6 %) as next T-300-independent target.

---

last_focus: 4
