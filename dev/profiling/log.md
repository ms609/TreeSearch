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

### Round 5 — 2026-06-18 — area #13 — FRAMING the gap (efficiency × throughput) on a FRESH post-fix build

**Why:** user directive — before VTune, test "is TNT really running 10× more
iter/sec, and is each iteration similar complexity?"  ALL prior rounds (1–4)
are stale: they predate the Wagner directional cost fix (2b299e4b) AND the
`unrooted=TRUE`-by-default flip (25e35be7).  Fresh `-O2` instrumented build
`.agent-p0` at **HEAD 25e35be7**.

- Drivers: `dev/profiling/drivers/framing.R` (TS converge vs local 32-bit TNT,
  xmult hits10/replic50, 4 datasets × seeds 1–2) + `dev/profiling/drivers/perclip.R`
  (per-clip economy via `ts_tbr_diagnostics`).  CSVs: `framing_latest.csv`,
  `perclip_latest.csv`.

**Unit-commensurability GATE (advisor) — PASSED.** `++n_evaluated` counts per
(reroot,regraft) candidate; TNT "rearrangements examined" is the same
granularity → cand_ratio is O(1) (buggy-era prior 1.2–1.9, headtohead_phase0.csv),
not O(n).  The efficiency factor is meaningful.

**Decomposition  wall_ratio(TS/TNT) = efficiency(ts_cand/tnt_rearr) × throughput(tnt_rate/ts_rate)**
(local same-machine identity; 32-bit TNT → throughput is a LOWER bound on the
true per-second gap; counts are bitness-free so efficiency is exact):

| dataset | tips | gapB | efficiency | throughput(32-bit) | local wall ratio |
|---------|------|------|-----------|--------------------|------------------|
| Wortley2006 | 37 | **0** | 0.49 | 1.36 | **0.68** |
| Zanol2014   | 74 | **0** | 0.95 | 2.30 | 2.17 |
| Zhu2013     | 75 | **0** | 1.08 | 1.84 | 1.98 |
| Giles2015   | 78 | **0** | 0.42 | 1.81 | **0.75** |

**Findings.**
1. **gapB = 0 on all four** — fresh build reaches TNT's exact optimum (quality
   parity holds; matched stopping rule satisfied — both converge to same score).
2. **Efficiency 0.42–1.08 — TS examines comparable-or-FEWER candidates than
   TNT.**  Full reversal from buggy-era (1.4–1.9).  The cost-bug + unrooted-TBR
   fixes made TS candidate-efficient.  Strategy is NOT the residual lever.
3. **Throughput 1.36–2.30× (32-bit, same machine) — NOT 10×.**  TNT examines
   rearrangements ~1.4–2.3× faster/sec; 64-bit TNT would widen this to perhaps
   2–4× (uncertain — needs a 64-bit local TNT to measure cleanly).  The historic
   "10–16×" gap was TNT-64bit-on-EPYC-wall vs TS-full-thorough-budget — NOT a
   like-for-like throughput comparison.  Like-for-like (matched 50-start effort,
   same machine) the gap is throughput-only and modest.
4. On Wortley/Giles TS is **faster than local 32-bit TNT** (wall 0.68/0.75).

**Per-clip economy (perclip.R, RandomTree start, TBR-to-convergence):**

| dataset | cand/clip | ns/cand | µs/clip |
|---------|-----------|---------|---------|
| Wortley2006 | 198 | 102 | 20 |
| Zhu2013 | 590 | 46 | 28 |
| Giles2015 | 763 | 42 | 31 |

cand/clip is LARGE (200–760) → per-clip overhead is amortized, BUT ns/cand
(42–102) materially EXCEEDS the bare bounded scan (~18–23 ns, stale Round 3) →
per-clip/per-reroot overhead is **roughly half** the per-candidate time on the
75–78-tip sets.  This overhead = the per-clip state TNT maintains incrementally
(Goloboff 1996) but TS RECOMPUTES every clip: `compute_insertion_edge_sets`
(NEW, O(n·chars) + 3 heap allocs/clip), `compute_from_above` (O(n·chars)),
`vroot_cache` rebuild, plus the unrooted reroot machinery (`fitch_join_states`
+ `fast_hash` + dedup per reroot).  **This IS the high-order "buffering/caching"
the user asked about, now quantified — it's the entire ~1.8× throughput gap.**

**Answer to the framing question:** NO, TNT is not running 10× more iter/sec
(~1.4–2.3× same-machine); each iteration IS of similar complexity (within ~2×);
the residual is throughput, and the high-order lever is per-clip state
recomputation that TNT amortizes via incremental maintenance — NOT line-level
scoring (the AVX2 scan is at the compiler limit per Round 3, re-confirmed by the
~18–23 ns floor here).

**VTune attribution (collected — `tbr-clip.R`, Zanol2014, 30 TBR-to-convergence;
symboled lib `.vtune-lib-20260618212528`; names via `nm`; TreeSearch.dll = 63%
of 6.5 s total CPU):**

| % total CPU | function | note |
|------------:|----------|------|
| **31.1** | **`ts::compute_insertion_edge_sets`** | **#1 — the NEW Wagner-fix code; 2.5× the next; ~49% of TreeSearch self-time. Inner `combine` is SCALAR.** |
| 12.4 | `simd::any_hit_reduce_avx2` | bounded scan — AT-LIMIT (Round 3) |
| 11.4 | `ucrtbase func@0x180020b2c` | heap/zero-fill — consistent with the per-clip `up`/`edge_set` `vector(N,0)`+`.assign(0)`+free churn |
| 5.9 | `fitch_indirect_length_cached` | scalar cached scorer |
| 4.5 | `tbr_search` (self) | candidate-loop control |
| 1.9 | `uppass_node` | AT-LIMIT |
| 1.2 | `build_postorder_prealloc` | per clip |
| 0.9 | `fitch_join_states` | per reroot (cheap) |
| — | `malloc_base` 2.4 + `memcpy` 1.0 | per-clip alloc/copy |

**The single dominant cost is the recently-added `compute_insertion_edge_sets`
(31% CPU + ~11–14% allocation/zero-fill).**  Three behaviour-neutral levers
(verify score-identical — this is the quality-fix function):
1. **Vectorize the combine** — its inner combine is scalar; the IDENTICAL
   intersect-else-union already exists as `simd::fitch_combine` (ts_simd.h:353,
   AVX2).  Swap it in.  Attacks the 31%.
2. **Hoist `up`/`pre` to caller-owned buffers** (emutls lesson T-S3d: plain
   locals, NOT thread_local) — kills the per-clip malloc.
3. **Skip the zero-fill** — every `up[D]`/`edge_set[D]` slot is written before
   read in preorder; the `vector(N,0)`/`.assign(0)` memset is unnecessary.
   Attacks the ucrtbase 11.4%.

Next: A/B all three in an isolated worktree; gate on score-identical across
datasets × {EW, IW} and wall delta on `framing.R`/`perclip.R`.

---

last_focus: 13

---

last_focus: 13

---

## Round 7 (2026-07-02, Opus) — Ratchet phase internals (T-P5c CLOSED)

Focus: the still-open T-P5c (ratchet internals). Lead: the ratchet
perturbed-weight search runs on the scalar scorer path (ts_tbr.cpp:2129) because
MIXED-mode upweight_mask disables use_flat+x4 — the one hot path the whole
flat/x4/edge-set/getenv program never optimised. thorough/intensive/large use
ratchetPerturbMode=2 (MIXED); .AutoStrategy routes the >=65t,>=100char roster there.

Method: env-gated std::chrono (TS_RATCHET_TIMING) around ratchet_search's 3
tbr_search calls + reject rebuild. Worktree TS-ratchet-prof off HEAD 3273d10b.
Driver dev/profiling/drivers/ratchet-internals.R. thorough, nThreads=1.

Result (Zanol2014 reps3, score 1261 = best-known): ratchet 1.322s / 2.33s wall = 57%.
Split: s1_init 4% / s2_perturb 10% / s3_restore 86% / reject ~0%. Zhu2013 same shape.

Verdict: CLOSE T-P5c as at-limit-by-inheritance.
- reject-branch per-cycle rebuild = 0.000s/120 cycles -> T-P5c's redundant-recompute
  concern CLEARED by measurement.
- s3 (86%) is full tbr_search-to-convergence on the FAST flat+x4 path = at-limit kernel.
- Fable lead (perturbed scalar path) = 5-7% of wall; weighted-flat kernel caps <2% wall
  and per-candidate delta argued-wash -> NOT built. This VINDICATES the prior
  kernel-reading dismissal of batch-scalar-x4-for-ratchet. getenv (T-P5n) was the
  exception: uncharacterised/profiler-hidden cost, not a characterised per-candidate delta.
- Only remaining wall lever = RECIPE (cycle count x restore depth); ~49% of wall is
  repeated full reconvergence after 5-move perturbations. Handed to recipe owners (#40).
  Byte-identical localized reconvergence DEAD (L3b non-locality + lever-c).

last_focus: ratchet-internals (T-P5c)
