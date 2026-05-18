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

last_focus: 2
