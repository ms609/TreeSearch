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

last_focus: 0
