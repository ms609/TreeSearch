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

### Round 3 — 2026-06-16 — area #13 (standard-Fitch TNT-parity path) — NEW AREA
First profile of the **standard-Fitch** path (inapplicable `-`→`?`, so
`has_na=FALSE`, flat/x4 kernels). Rounds 1-2 profiled the NA three-pass path
on raw `inapplicable.phyData`; that path is ~20× slower per replicate with an
entirely different hotspot mix (`fitch_na_*` dominate). Standard Fitch is the
path the TNT benchmark actually compares against.

- Driver:        dev/profiling/drivers/fitch-tnt.R   (bare: 5.57 s / 8 reps = 0.56 s/rep; Zhu2013 75t, auto→thorough, nThreads=1, score 627 vs TNT 624)
- Build:         dev/profiling/.vtune-lib-20260616052323 (HEAD 841eead3, -O2 -g)
    ⚠ GOTCHA: default Windows R build STRIPS the DLL (`DLLFLAGS=-s` in Makeconf)
    → VTune shows `func@0x…`/`[Unknown]`. Override `MAKEFLAGS="DLLFLAGS=-static-libgcc"`
    to drop `-s`; verify `objdump -h DLL | grep debug_info` + `nm DLL` (23089 syms).
    ⚠ GOTCHA2: even symboled, VTune's CSV reporter emits `func@0x…` (MinGW DWARF
    unparsed). Resolve via `nm -C DLL` — image base 0x2cc1a0000 is stable across
    builds, so VTune addresses map 1:1 to nm addresses.
- profvis/Rprof:  99.5 % self-time in `ts_driven_search` (single .Call); R <0.5 %; no [Port].
- Phase dist (attr "timings", 3 reps): ratchet 63.0 %, rss 9.2 %, xss 9.2 %,
    css 6.8 %, wagner 5.5 %, tbr 4.0 %, final_tbr 2.3 %, drift/nni/anneal 0 %.
- VTune top fns (TreeSearch.dll self, total 2.70 s; names via nm):
    1. ts::tbr_search (orchestration, 2 ranges)   25.1 %  — candidate-loop control + collapsed/sector vector<bool> bit-tests + inlined scoring
    2. ts::simd::any_hit_reduce_avx2              14.5 %  — core 2-op Fitch reduce
    3. ts::uppass_node                            13.2 %  — incremental uppass; SCALAR state-update loop (cf. vectorised fitch_combine)
    4. ts::simd::any_hit_reduce3_avx2              6.3 %  — 3-op reduce (SPR bounded)
    5. ts::TreeState::build_postorder_prealloc     5.2 %  — O(n) postorder rebuild, per clip AND per accept
    6. ts::fitch_incremental_downpass              4.1 %
    7. ts::fitch_indirect_bounded_flat             4.0 %
    8. hash_tree / fitch_indirect_length_cached (scalar) / validate_topology  ~2.9 % each
   (Scoring SIMD+wrappers ≈ 50 %; per-clip bookkeeping ≈ 18 %; orchestration 25 %.)
- Findings:
    [AT-LIMIT] SIMD `any_hit_reduce` (21 %): disasm of `hor_or256` shows GCC
      already elides the store-reload (vextracti128/vpsrldq/vpor/vmovq, register-only)
      → compiler-optimal. No win.
    [AT-LIMIT] `uppass_node` vectorisation (13 %): micro-bench
      dev/profiling/microbench/bench_uppass_combine.cpp — AVX2 update loop
      bit-identical (value+changed flag) but only **1.22×** at n_states=4, and the
      4-wide path does NOT trigger for 2-state (binary) morph chars → ~1 % wall,
      not worth the incremental-uppass correctness risk.
    [Optimise, modest] Per-clip/accept allocation churn (~3-4 %): compute_from_above,
      collect_main_edges/collect_subtree_edges, validate_topology heap-alloc
      std::vector scratch per call (_M_realloc_append 1.1 %). Extend existing
      prealloc pattern (work_stack/saved_postorder/clip_actives_buf). Low risk.
    [Strategic] Standard-Fitch is **bookkeeping- + strategy-bound**, not
      scoring-bound. Per-candidate scoring is at the AVX2/compiler limit; remaining
      levers are (a) reduce per-clip O(n) bookkeeping (postorder rebuild +
      incremental passes ≈ 18 %), (b) ratchet evaluation economy (63 %). Both align
      with the TNT-outperformance analysis (strategy > code). Score near parity.
- Filed:         findings.md row T-S3a (allocation churn) + AT-LIMIT rows.
- Cleanup:       result_fitch_tnt_* + result_fitch_sym_* removed; stripped lib
    .vtune-lib-20260616051420 removed; symboled lib kept pending follow-up. microbench kept.
- Next reviewer: code lever for parity is per-clip bookkeeping (incremental
    postorder across clip/unclip), NOT the scoring kernel. Strategy lever: ratchet
    eval economy (time-adjusted expected-best).

---

### Round 4 — 2026-06-17 — area #13 (standard-Fitch) — StateSnapshot re-profile + build-protocol hardening
- Trigger: re-confirm the stale "StateSnapshot ~23%" before the deferred
  selective save/restore surgery (task #10), on a FRESH symboled build.
- Driver:  dev/profiling/drivers/fitch-tnt.R (Zhu2013, 12 reps; bare ~3.5 s)
- Build:   **build-symboled-lib.ps1 (NEW, in /profile skill dir)** — isolated
  tarball (src/ untouched, concurrent-safe) + PKG_CXXFLAGS `-g -fno-omit-frame-pointer`
  + DLLFLAGS=-static-libgcc; HARD-FAILS if no .debug_info. 23,221 syms.
  ⚠ GOTCHA caught: the prior `Makevars.vtune` set `CXXFLAGS`, which R SILENTLY
  BYPASSES for C++17 (uses CXX17FLAGS) → only ~214 KB .debug_info (partial; -g
  on cache-hit TUs only, -fno-omit-frame-pointer absent). PKG_CXXFLAGS fixes it
  → 19 MB .debug_info (all TUs). resolve_syms.R maps VTune `func@0x` via `nm`.
- VTune top (resolved, % of total CPU): ts::tbr_search 12.8 %, simd any_hit_reduce
  ×2 = 14.2 % (AT-LIMIT), uppass_node 7.7 % (AT-LIMIT), memcpy 5.0 %, malloc+free
  6.4 %, fitch_indirect* ~17 %, build_postorder 3.1 %, save_node_state 2.5 %,
  hash_tree 2.3 %.
- Finding: **[AT-LIMIT] selective StateSnapshot (task #10) — NOT worth it.**
  save/restore (save_node_state 2.5 % + its memcpy share) ≈ 3-5 % ceiling;
  selective restore reclaims ~half → ~2 % for risky surgery on the most
  correctness-critical code. The cited "23 %" was stale NA-path (pre-T-261/T-300).
  Per-candidate cost is at floor; ~19 % is alloc/copy churn. The 2× gap is
  candidates-per-improvement (SEARCH STRATEGY) → see TNT sectorial reverse-
  engineering: [[tnt-sectorial-recipe]] memory + dev/benchmarks/tnt_sector_defaults.csv.
- Cleanup: result_statesnap_* removed; old partial-symbol libs removed; kept
  .vtune-lib-20260617081344 (validated symboled). build-symboled-lib.ps1 +
  resolve_syms.R retained.
- Next reviewer: profile the NEW multi-start sectorial once built (regenerate
  the symboled lib with build-symboled-lib.ps1 first — never reuse a stale one).

---

last_focus: 13
