# Profiling Expertise — TreeSearch

## Purpose

Profile the C++ search engine to identify bottlenecks. Produce specific,
actionable optimization tasks for `to-do.md`.

## Tools

### 1. Built-in Phase Timing (Quick)

The driven search already has `std::chrono` phase timing at `verbosity >= 2`.
Use the R-level interface:

```r
library(TreeSearch)
library(TreeTools)
dataset <- TreeSearch::inapplicable.datasets[["Vinther2008"]]
result <- MaximizeParsimony(dataset, maxReplicates = 3, verbosity = 2L)
```

This prints per-phase timing. For programmatic access, use the
`ts_bench_tbr_phases` diagnostic function (7 args, registered in
TreeSearch-init.c).

### 2. std::chrono Micro-Benchmarks (Medium)

For fine-grained timing of specific functions, add `steady_clock` timing
around the code path of interest. See `inst/benchmarks/bench_memory.R`
and `inst/benchmarks/bench_simd.R` for examples.

Key metrics to measure:
- Per-candidate indirect scoring cost (ns)
- Clip+incremental phase time (μs per TBR pass)
- Full rescore time (μs)
- Snapshot save/restore time (μs)

### 3. VTune (Thorough)

For instruction-level hotspot analysis, use the `r-package-profiling`
skill (load via the skill tool). Key steps:

1. Build with debug symbols: set `DLLFLAGS` via `MAKEFLAGS` env var
2. Run a representative workload under VTune
3. Analyze hotspots in the VTune GUI

See `.positai/skills/r-package-profiling/references/` for detailed
VTune workflow on Windows.

**Current version: VTune 2025.10** (updated 2026-03-19). Requires Ice Lake
or newer CPU (10th gen Intel Core / 3rd gen Xeon Scalable+). VS 2019
integration and Eclipse integration are removed in 2025.x. Command-line
workflow (`vtune -collect hotspots`) is unchanged.

### 4. R-Level Profiling

For R overhead identification:

```r
Rprof("profile.out")
result <- MaximizeParsimony(dataset, maxReplicates = 5)
Rprof(NULL)
summaryRprof("profile.out")
```

## Known Baselines

### Latest run: 2026-03-27 by Agent A (round 6: post-T-261/T-262/T-263 phase distribution)

See "Phase distribution: current thorough preset" section below for updated numbers.
The 2026-03-18 baselines used strategy='none' (TBR-only); the thorough preset
now dominates medium-scale search, making direct comparison impractical.

### Previous run: 2026-03-18 16:00 by Agent A (v2.0.0, single-agent, quiet machine)

Previous baselines (2026-03-17) were inflated ~30–40% by multi-agent machine
contention. Scores are identical. Timings below are authoritative.

### End-to-end benchmarks (3-run medians, 5 reps, strategy='none', EW):

| Dataset | Tips | Chars | Median (s) | Score |
|---------|------|-------|------------|-------|
| Vinther2008 | 23 | 57 | 0.390 | 79 |
| Agnarsson2004 | 62 | 242 | 1.860 | 778 |
| Zhu2013 | 75 | 253 | 2.720 | 655 |
| Dikow2009 | 88 | 220 | 3.860 | 1614 |

### Per-phase breakdown (Zhu2013, 5 reps, two runs averaged):

| Phase | % of time | Avg ms/rep |
|-------|-----------|------------|
| Wagner | <0.1% | <1 |
| TBR | 24–37% | 110–160 |
| XSS | 10% | 35–55 |
| RSS | 2% | 9–13 |
| Ratchet | 24–28% | 90–155 |
| Drift | 25–33% | 90–200 |
| Final TBR | 2% | 7–10 |

Ratchet (24-28%) and drift (25-33%) dominate. TBR (24-37%) varies
substantially by run. XSS ~10%, RSS ~2%, both stable.

### Wagner tree construction: Negligible (<0.1% of search time)

| Dataset | Tips | µs/tree | % of replicate |
|---------|------|---------|----------------|
| Vinther2008 | 23 | 300 | <0.1% |
| Agnarsson2004 | 62 | 1000 | 0.3% |
| Zhu2013 | 75 | 600 | 0.1% |
| Dikow2009 | 88 | 1400 | 0.2% |

Not a bottleneck at any dataset size. No optimization needed.

### Parallel scaling (2 threads)

| Dataset | Reps | 1T (s) | 2T (s) | Speedup | Efficiency |
|---------|------|--------|--------|---------|------------|
| Zhu2013 | 5 | 2.53 | 1.59 | 1.59× | 80% |
| Zhu2013 | 10 | 5.16 | 3.29 | 1.57× | 78% |
| Zhu2013 | 20 | 10.70 | 5.20 | 2.06× | 103%* |
| Zhu2013 | 40 | 18.63 | 11.35 | 1.64× | 82% |
| Dikow2009 | 10 | 7.76 | 5.11 | 1.52× | 76% |

*Superlinear at 20 reps is stochastic noise (different search paths).

**Finding:** Typical 2-thread efficiency is 78–82%. The old 1.24× measurement
was a multi-agent machine contention artifact. The implementation (dynamic
work-stealing via `atomic::fetch_add`, mutex-guarded pool) is sound.
Main loss is stochastic load imbalance between replicate times.

### XSS/RSS effectiveness (5 reps per dataset)

| Dataset | Tips | XSS hits | XSS avg Δ | XSS avg ms | RSS hits | RSS avg Δ | RSS avg ms |
|---------|------|----------|-----------|------------|----------|-----------|------------|
| Agnarsson2004 | 62 | 3/5 | 3.8 steps | 59 | 0/5 | 0 | 14 |
| Zhu2013 | 75 | 5/5 | 26.6 steps | 43 | 2/5 | 1.0 | 11 |
| Dikow2009 | 88 | 0/5 | 0 | 93 | 1/5 | 3.2 | 29 |

**Finding:** XSS effectiveness is highly dataset-dependent — from zero
improvement (Dikow2009) to 27-step average improvement (Zhu2013). No obvious
predictor from simple nTip/nChar statistics. XSS cost is ~10% of replicate
time; acceptable when effective but wasted when not.

RSS is marginal across all datasets (0–3 steps, 2% of time). One exception:
Dikow2009 where RSS found 16 steps while XSS found 0 — suggests they
explore different neighbourhoods.

### Auto strategy (reference — unchanged from T-066/T-068 study)

Threshold: ≥75 tips AND nChar < 100 triggers "thorough". Signal-density gate
prevents unnecessary thorough runs on character-rich datasets.

### R overhead: <0.5% of wall time (confirmed via Rprof, unchanged)

### Scaling exponent: ~2.82 (TBR pass time vs tips, unchanged)

### Drift/ratchet cycle tuning (reference — unchanged from T-029 study)

| Config | Med score | Min score | Med time | Speedup |
|--------|-----------|-----------|----------|---------|
| d5_r5 (default) | 656 | 648 | 5.7s | — |
| d2_r5 | 660 | 646 | 4.1s | 28% |
| d2_r2 | 662 | 656 | 3.8s | 33% |
| d0_r5 | 658 | 650 | 2.8s | 51% |
| d5_r0 | 662 | 660 | 4.8s | 16% |

Lower score = better. Current defaults: d2_r5.

### CSS effectiveness: Marginal (adds 2-6% time, no consistent improvement)
Disabled by default (cssRounds=0).

### Latest EW regression check: 2026-03-19 by Agent A (v2.0.0, post T-115–T-124)

All datasets pass regression benchmark. EW baselines updated with 7-run medians:

| Dataset | Tips | Chars | Median (s) | Score (range) | Notes |
|---------|------|-------|------------|---------------|-------|
| Vinther2008 | 23 | 57 | 0.420 | 79 | stable |
| Agnarsson2004 | 62 | 242 | 1.790 | 778 | stable |
| Zhu2013 | 75 | 253 | 3.170 | 648–666 | high variance (2.5–7.6s range) |
| Dikow2009 | 88 | 220 | 4.900 | 1612–1614 | high variance (4.0–12.4s range) |

Zhu2013/Dikow2009 appear slightly slower than 2026-03-18 baselines (~17–27%) but
within stochastic noise. Phase breakdown unchanged. No regression in C++ engine.
The recent DataSet changes (inapp_state field, HSJ/XFORM modes) have no measurable
effect on EW search paths.

### HSJ and XFORM scoring baselines: 2026-03-19 by Agent A

Synthetic hierarchical datasets (valid hierarchy structure: primary + secondary chars,
secondaries are inapplicable when primary absent). 3-run medians, 5 reps per run.

| Config | Tips | Chars | Blocks | EW (s) | HSJ (s) | XFORM (s) | HSJ/EW | XFORM/EW |
|--------|------|-------|--------|--------|---------|-----------|--------|----------|
| small | 20 | 19 | 3 | 0.020 | 0.010 | 0.020 | 0.5× | 1.0× |
| medium | 40 | 50 | 5 | 0.170 | 0.100 | 0.280 | 0.6× | 1.6× |
| large | 60 | 82 | 8 | 0.610 | 0.360 | 1.330 | 0.6× | 2.2× |
| xlarge | 80 | 120 | 10 | 5.920 | 3.560 | 9.460 | 0.6× | 1.6× |

**HSJ is faster than EW** (~0.6× at medium/large sizes) because:
1. Fitch candidate screening guards expensive full HSJ rescore — most candidates
   are rejected by Fitch before HSJ is called.
2. Hierarchy datasets have a simpler parsimony landscape (secondaries add signal
   only when primary is present), leading to faster search convergence.

**XFORM is slower than EW** (~1.6–2.2× at medium/large sizes) due to Sankoff
cost per candidate. Phase breakdown (large config, 5 reps):

| Phase | EW avg ms/rep | HSJ avg ms/rep | XFORM avg ms/rep |
|-------|---------------|----------------|------------------|
| TBR | 25 | 23 | 29 |
| XSS | 14 | 7 | 14 |
| RSS | 4 | 2 | 5 |
| Ratchet | 51 | 28 | 86 |
| Drift | 22 | 13 | 36 |
| Final TBR | 2 | 1 | 4 |
| **Total** | **117** | **74** | **174** |

XFORM overhead concentrated in Ratchet (+69%) and Drift (+64%), which perform
more scoring iterations than TBR. XSS/RSS overhead is negligible.

**Conclusion:** Both modes are acceptable. XFORM at ~1.7× overhead for real
workflows is reasonable given the algorithmic complexity (Sankoff vs Fitch).
No optimization tasks raised — XFORM at this cost is expected behavior.

### Hierarchical resampling: 2026-03-19 by Agent A

Medium config (40 tips, 50 chars, 5 blocks), jackknife, 20 reps:

| Mode | 1 thread (s) | 2 threads (s) | Speedup |
|------|-------------|--------------|---------|
| Brazeau (C++ parallel) | 5.19 | 2.05 | 2.5× |
| HSJ hierarchical (serial R loop) | 1.76 | 1.64 | 1.1× |
| XFORM hierarchical (serial R loop) | measured via 10-rep: ~1.58 | — | — |

**Finding 1 (positive):** HSJ/XFORM hierarchical resampling is faster than Brazeau
per-replicate because the block-level resampling units (35 vs 50 units) produce
simpler per-replicate datasets. No performance concern here.

**Finding 2 (known limitation):** Hierarchical resampling uses a serial R loop
across replicates — `nThreads` only applies within each replicate's internal search.
Brazeau gets full 2.5× at 2 threads; HSJ/XFORM get only ~1.1×. For users running
50–100 jackknife replicates with large HSJ/XFORM datasets, wall time will be ~2×
longer than equivalent Brazeau. This is documented in AGENTS.md as a known future
optimization (C++-level inter-replicate parallelism for hierarchical resampling).
No new task filed — already on the roadmap.

### Preset tuning benchmark: 2026-03-22 by Agent A

Compared updated presets (wagnerStarts=3, sprFirst=TRUE, adaptiveLevel=TRUE
for default; wagnerStarts=3, sprFirst=TRUE for thorough) against old presets
(wagnerStarts=1, sprFirst=FALSE, adaptiveLevel=FALSE). 7-run medians via
`MaximizeParsimony()`, strategy=auto, 10 reps, 1 thread.

| Dataset | Tips | Preset | Old time (s) | New time (s) | Δ time | Old score | New score |
|---------|------|--------|-------------|-------------|--------|-----------|-----------|
| Vinther2008 | 23 | sprint | 0.76 | 0.65 | –14% (noise) | 79 | 79 |
| Agnarsson2004 | 62 | default | 3.59 | 2.41 | **–33%** | 778 | 778 |
| Zhu2013 | 75 | thorough | 23.65 | 24.83 | +5% (noise) | 647 | 648 |
| Dikow2009 | 88 | thorough | 49.19 | 39.24 | **–20%** | 1611 | 1612 |

**Findings:**
- `adaptiveLevel` in `default` preset: consensus-stability triggers early exit
  on easy landscapes (Agnarsson2004), saving 33%. No score regression.
- `sprFirst + wagnerStarts=3` in `thorough`: 20% faster on Dikow2009 (better
  starting tree reduces initial TBR descent). Neutral on Zhu2013.
- **Do not enable `adaptiveLevel` in `thorough`**: with 20 ratchet + 12 drift
  base, 1.5× scaling creates 30 ratchet + 18 drift per hard replicate,
  causing 3–4× slowdowns for only 2–3 step improvement (benchmarked separately).

### 180-tip large-preset baselines: 2026-03-26 by Agent E (Hamilton HPC, EPYC 7702)

Dataset: mbank_X30754 (180 taxa, 425 chars, 418 patterns, 40% missing, 20% inapplicable).
Strategy: auto → "large" preset. 5 seeds per budget, single-threaded.

**Score quality by budget (median, 5 seeds):**

| Budget | Median score | Range | Reps/seed |
|--------|:-----------:|:-----:|:---------:|
| 30s | 1202 | 1189–1214 | ~1.5 |
| 60s | 1190 | 1190–1202 | ~3 |
| 120s | 1185 | 1171–1189 | ~6 |

Per-replicate time: median 17.3s (range 13.7–21.2s). MPT enumeration adds
0–2 steps beyond best single-replicate score.

**Phase distribution (rep 1, 30s budget, 5-seed averages):**

| Phase | % time | Mean ms | Steps/s | Hit rate |
|-------|:------:|--------:|:-------:|:--------:|
| TBR | 43.6% | 7313 | 91.4 | 5/5 (661 steps avg) |
| Ratchet | 32.2% | 5390 | 4.5 | 5/5 (26.6 steps avg) |
| SA (anneal) | 7.4% | 1241 | 0.8 | 7/50 (14%, 1.3 steps) |
| XSS | 5.4% | 897 | 13.8 | 4/5 |
| Wagner+NNI | 4.7% | 790 | — | starting point |
| RSS | 3.2% | 530 | 4.8 | 3/5 |
| CSS | 2.5% | 424 | 11.2 | 2/5 |
| Final TBR | 1.0% | 174 | 5.2 | 1/5 |

**SA (simulated annealing) phase is the least productive:** 7.4% of time,
14% hit rate (7/50 reps improved by 1.3 steps on average). Efficiency =
0.8 steps/s, far below ratchet (4.5) or XSS (13.8). annealCycles=3,
annealPhases=5 may be overtuned. Reducing could save ~1.2s/rep → 1 extra
replicate per ~17s saved.

**Comparison with earlier Intel desktop baselines (T-179, pre-T-206):**

| Budget | Intel (pre-T-206) | EPYC (post-T-206) | Delta |
|--------|:-:|:-:|:-:|
| 30s | 1276 | 1202 | −74 |
| 60s | 1255 | 1190 | −65 |
| 120s | 1250 | 1185 | −65 |

The 65–74 step gap is **primarily due to T-206** (outer cycle reset cap),
not hardware. T-206 was merged 2026-03-24 19:27; the Intel baselines were
recorded at 12:56 the same day (pre-T-206). Without the reset cap, each
replicate performed 3–5 pipeline cycles (~51–85s) vs ~17s with cap=0.
At 120s budget: ~2 replicates pre-T-206 vs ~6 post-T-206. Hardware
differences (Intel desktop vs EPYC 7702) are a secondary factor.

### Phase distribution: current thorough preset (2026-03-27, Agent A, round 6)

Dataset: Zhu2013 (75t, 253 chars). Strategy: auto → thorough.
3 reps, single-threaded, post-T-261+T-262+T-263. Total: 33.7 s = ~11.2 s/rep.

| Phase | Calls | Total ms | Mean ms | % |
|-------|:-----:|:--------:|:-------:|:---:|
| Ratchet | 14 | 15617 | 1116 | 46.3% |
| NNI-perturb | 14 | 11565 | 826 | **34.3%** |
| RSS | 14 | 2488 | 178 | 7.4% |
| CSS | 14 | 1477 | 106 | 4.4% |
| XSS | 14 | 1079 | 77 | 3.2% |
| TBR (post-phase) | 14 | 622 | 44 | 1.8% |
| Initial TBR | 3 | 468 | 156 | 1.4% |
| wag+NNI | 2 | 427 | 214 | 1.3% |

**Key findings vs 2026-03-18 baselines:**

1. **TBR is no longer a bottleneck** (1.4% + 1.8% = 3.2%). T-261+T-262+T-263
   combined are working — TBR has become fast enough that other phases dominate.
   Drift was 25–33% before T-255; its removal freed that budget to more ratchet.

2. **NNI-perturb at 34.3% with poor efficiency:**
   - Hit rate: 14% (2/14 calls improved score)
   - Mean improvement when hit: 1 step
   - Efficiency: 0.17 steps/s vs ratchet's ~4–8 steps/call at comparable cost
   - Cost grows within a replicate (early calls ~300ms, late calls ~1300ms)
   - This phase likely over-tuned for 75-tip datasets. Filed **T-274** (P2).

3. **RSS at 7.4%** — higher than old 2% baseline. With conflict-guided RSS and
   outerCycles/reset mechanism creating ~4.7 RSS calls per replicate at ~178ms each
   (~837ms/rep). Old uniform RSS: ~11ms/rep. 16× overhead increase. Most of this
   is the actual sector TBR cost (more calls × similar per-sector time), not conflict
   computation overhead. The reset mechanism is the multiplier.

4. **wag+NNI at 1.3%**: biased Wagner + 3 starts + NNI warmup adds ~214ms per
   replicate start. Negligible at this scale; confirms T-246/NNI-warmup tuning is fine.

### Round 7 — 2026-06-16 — STANDARD-FITCH path (TNT-parity), first profile

**Critical distinction:** all rounds above profiled the **NA three-pass path**
(raw `inapplicable.phyData`). The TNT benchmark compares the **standard-Fitch**
path: inapplicable `-` replaced with `?` so `has_na=FALSE` and the flat / 4-wide
(T-245) kernels run — NOT `fitch_na_*`. This path is **much cheaper per replicate**
(0.56 s/rep here) with a *completely different* hotspot mix — the NA three-pass adds
~3-4× scoring work plus subtree_actives bookkeeping. (The round-6 NA ~11 s/rep figure
is NOT a clean comparison — it also ran now-disabled phases, NNI-perturb + outer resets.)
Driver:
`dev/profiling/drivers/fitch-tnt.R` (Zhu2013 `-`→`?`, auto→thorough, nThreads=1).
Score 627 vs TNT 1.6 = 624 (was +11 in March; now +3 — near parity).

Phase distribution (standard Fitch, attr "timings"): **ratchet 63%**, rss 9.2%,
xss 9.2%, css 6.8%, wagner 5.5%, tbr 4%, final_tbr 2.3%. (Contrast the NA/thorough
round-6 table above — there ratchet 46% + NNI-perturb 34%; NNI-perturb now off.)

VTune top functions (TreeSearch.dll self, total 2.70 s; names via `nm`):
| % DLL | function | note |
|------:|----------|------|
| 25.1 | `tbr_search` (orchestration self) | candidate-loop control + `vector<bool>` collapsed/sector bit-tests + inlined scoring |
| 14.5 | `simd::any_hit_reduce_avx2` | **AT-LIMIT** — disasm: GCC already elides the `hor_or256` store-reload |
| 13.2 | `uppass_node` | scalar update loop; **AT-LIMIT** — vectorising = 1.22× only, nil for 2-state (micro-bench) |
|  6.3 | `any_hit_reduce3_avx2` | SPR-bounded reduce |
|  5.2 | `build_postorder_prealloc` | O(n) rebuild **per clip AND per accept** |
|  4.1 | `fitch_incremental_downpass` | per clip |
|  4.0 | `fitch_indirect_bounded_flat` | SPR candidate scoring |
| ~2.9 | `hash_tree` / `fitch_indirect_length_cached` (scalar) / `validate_topology` | |

**Conclusion:** per-candidate scoring is at the AVX2/compiler limit. Standard-Fitch
is **bookkeeping- + strategy-bound**: per-clip O(n) work (postorder rebuild +
incremental passes + edge/from_above/vroot construction ≈ 18%) and ratchet (63%)
are the levers — exactly what TNT minimises (Goloboff 1996). This corroborates
`.positai/plans/2026-03-21-tnt-outperformance-analysis.md` (strategy > code).

**Two build gotchas for VTune on this MinGW/Windows toolchain** (cost an hour):
1. The default R Windows build **strips** the DLL (`DLLFLAGS=-s` in Makeconf) →
   VTune shows `func@0x…`/`[Unknown]`. Override with
   `MAKEFLAGS="DLLFLAGS=-static-libgcc"` (drops `-s`); the personal `~/.R/Makevars`
   supplies `-g` via `CXXFLAGS` (and zeroes `PKG_CXXFLAGS`, so a package
   `src/Makevars.win PKG_CXXFLAGS` is ignored). Verify: `objdump -h DLL | grep debug_info`.
2. Even with symbols, VTune's CSV reporter still prints `func@0x…` (MinGW DWARF
   unparsed). Resolve names with `nm -C DLL`; the image base (0x2cc1a0000) is stable
   across rebuilds, so VTune addresses map 1:1 to `nm` addresses. `addr2line` gives
   the source FILE but not always the line (DWARF5 line tables).

Full round notes + ranked candidates: `dev/profiling/log.md` (Round 3),
`dev/profiling/findings.md` (T-S3a/b/c/d), `dev/profiling/baselines.md`.

**Wins shipped this round (per-clip allocation churn).** Tier 1 (T-S3a, committed):
hoist per-clip scratch (`dirty`, DFS stacks, preorder) to reusable buffers ≈ **4%**.
Tier 2 (T-S3d, uncommitted): replace the per-clip `unordered_set<uint64_t>`
rerooting-dedup with a reusable open-addressed table (generation-stamped O(1) reset)
≈ **3%**. Both behaviour-neutral (score 627 unchanged; scores identical on 6 datasets
× {NA, standard} via `dev/profiling/bench_equiv.R`). On the fast standard-Fitch path,
per-clip *allocation* is a bigger proportional share than on the NA path (T-260 saw
the dedup destructor at ~0.4%, but the full alloc footprint — bucket array + a node
malloc per insert + teardown, every internal-node clip — is several %).

**emutls gotcha (cost two rebuilds — record for next time).** Tier 2 was *neutral*
when first written as `static thread_local`: MinGW resolves a `thread_local` living
in a **dynamically-loaded DLL** via **emutls** = a function call (`__emutls_get_address`)
on *every access*. For a buffer touched once per clip (Tier 1) that's invisible, but
the dedup table's `insert()` runs *per reroot candidate*, so the per-access TLS call
cancelled the allocation saving. Fix: a **plain local declared once before the loop**
— in a per-thread-entered function (each search thread calls `tbr_search` with its own
`TreeState`) a stack local is already per-thread-safe and has *zero* TLS cost. Rule:
prefer hoisted plain locals over `static thread_local` on per-inner-iteration hot paths;
reserve `thread_local` for buffers inside helpers you can't hoist (Tier 1's case).

**Verifying from a temp lib (two traps).** (1) `test_dir()` switches CWD to
`tests/testthat/`, so a **relative** `lib.loc` breaks package-data lazy-load
(`cannot open .../data/Rdata.rdb`) even though `library()` loaded fine — pass an
**absolute** `lib.loc`. (2) There is a **flaky `nThreads=2` crash** (exit 127, no R
error) in the parallel test files under `Rscript`; it reproduces on the *unmodified*
build too (not change-specific) and passes on re-run, so the full `ts-` suite is not a
reliable gate from a temp lib — verify behaviour-neutrality with a single-threaded
score-equivalence harness instead (`bench_equiv.R`). [Flaky parallel crash: open
question — real concurrency bug vs Rscript/sandbox thread-resource artifact.]

## What to Profile

Status key: ✅ resolved, ⚠ partially explored, ❌ not yet investigated

1. ✅ **Drift + ratchet inner loops** (50–60% of C++ time combined). Both use
   TBR internally. Per-candidate indirect evaluation at memory-throughput
   limit (~23 ns at 75 tips per T-075). Cycle counts tuned (d2_r5).
   **Drift threshold sensitivity (2026-03-18 Agent E):** AFD={1,3,5,8} ×
   RFD={0.05,0.1,0.2} on Zhu2013 (75 tips, 15 runs each): no significant
   score difference between any config (Wilcoxon p=0.60–1.00). Permissive
   thresholds (AFD=8, RFD=0.2) waste time; tight vs default indistinguishable.
   On Dikow2009 (88 tips), d2 drift provides no benefit over ratchet alone
   (p=0.54); d6 gives 2-step improvement (p=0.006) at 2× time cost.
   **Conclusion:** Current defaults (AFD=3, RFD=0.1) are fine. Cycle count
   matters more than threshold values. No optimization task raised.

2. ✅ **Sectorial search effectiveness** (12% of time). XSS effectiveness is
   dataset-dependent (0–27 steps). RSS is marginal (0–3 steps). No clear
   predictor from simple dataset statistics. Could make XSS adaptive (skip
   after N unproductive reps) but time savings would be <10%.

3. ✅ **Wagner tree construction**: <0.1% of search time. Not a bottleneck.

4. ✅ **R overhead**: <0.5% of wall time. Not a bottleneck.

5. ✅ **Parallel scaling**: 78–82% efficiency at 2 threads. Implementation is
   sound (dynamic work-stealing, low-contention pool). Main loss is stochastic
   load imbalance. No obvious improvement without algorithmic changes.

6. ✅ **IW scoring overhead** (2026-03-18 Agent E). Compared EW vs IW (k=10,
   k=3) on three datasets (5 runs each, d2_r5, 5 reps, serial):
   - Vinther2008 (23 tips): IW 64% *faster* (landscape converges quicker)
   - Agnarsson2004 (62 tips): IW 26–39% slower
   - Zhu2013 (75 tips): IW 40–57% slower
   IW overhead scales with dataset size due to per-character weighted delta
   computation in indirect scoring. No optimization opportunity — the delta
   lookup is already O(n_blocks) per candidate, same as EW Fitch.

7. ✅ **Fuse effectiveness** (2026-03-18 Agent E). Compared fuseInterval=0 vs
   3 on three datasets (8 runs each, 10 reps):
   - Agnarsson2004: identical scores/time (pool deduplicates to 1 tree)
   - Zhu2013: identical scores/time
   - Dikow2009: negligible overhead (13.65s vs 13.78s with poolSuboptimal=5)
   Fuse is cheap when pool is small, free when pool=1. Current default
   (fuseInterval=3) is appropriate. No optimization task raised.

## Comparing Search Strategies: Time-Adjusted Expected Best

When comparing strategies that differ in per-replicate cost (e.g. NNI→TBR
vs TBR alone), the **median per-replicate score is the wrong metric**.
Multi-start search keeps the best tree across all replicates, so what
matters is the expected minimum from k independent draws, where
k = budget / time_per_replicate.

A strategy with high variance but occasional excellent scores can dominate
a consistent-but-mediocre one — if it's fast enough to get more draws.

**Bootstrap estimation:**
```r
expected_best <- function(scores, k, n_boot = 5000) {
  mean(replicate(n_boot, min(sample(scores, k, replace = TRUE))))
}

# k = budget / median_time_per_rep for each strategy
k <- floor(budget / median_time)
exp_best <- expected_best(observed_scores, k)
```

Compare `exp_best` across strategies at fixed budget (e.g. 20s, 60s, 120s).
This naturally trades off per-replicate quality against replicate throughput.

**When median IS acceptable:** comparing parameter changes on a fixed pipeline
(same time-per-rep), e.g. ratchet perturbation probability. All runs take
roughly the same time, so k is constant and the median is a reasonable proxy.

See AGENTS.md "NNI in the driven pipeline" for the reference application of
this metric (NNI→TBR vs TBR at 88 and 180 tips).

## Reporting Format

For each finding, add to `to-do.md`:

```
| T-NNN | P2 | OPEN | — | [Profile] Brief description | X% of time. Potential Y% improvement via Z approach. |
```

Include the measurement methodology and baseline numbers so the implementer
can verify the improvement.

8. ✅ **HSJ scoring overhead** (2026-03-19 Agent A). HSJ is ~0.6× EW wall time
   (faster) on synthetic hierarchical data. Fitch screening gates full HSJ rescore
   effectively. No optimization needed.

9. ✅ **XFORM (Sankoff) scoring overhead** (2026-03-19 Agent A). XFORM is ~1.6–2.2×
   EW wall time. Overhead concentrated in Ratchet (+69%) and Drift (+64%). This
   is expected Sankoff vs Fitch arithmetic cost — no obvious optimization target.

10. ✅ **Hierarchical resampling parallelism** (2026-03-19 Agent A). Serial R loop
    means `nThreads` only applies within each replicate. Brazeau 2T = 2.5× speedup;
    HSJ/XFORM hierarchical 2T = 1.1× only. Known limitation, future optimization
    (C++-level inter-replicate parallelism for hierarchical resampling).

11. ✅ **MaddisonSlatkin internal bottlenecks** (2026-03-19 Agent A, T-149).
    VTune hotspot collection (software sampling, `-g -fno-omit-frame-pointer`
    symbols build) on 57 calls at boundary cases: k=3/n=20–25, k=4/n=14–18,
    k=5/n=9–12. Total ~23 s CPU time; 63% in `TreeSearch.dll`.

    **CPU time breakdown within TreeSearch.dll (14.1 s):**

    | Category | CPU (s) | % DLL |
    |----------|---------|-------|
    | `logB_cache::find` (k=3,4,5) | 2.72 | 19% |
    | `SolverT<N>::LogB` compute | 1.88 | 13% |
    | `logPVec_cache::find` (k=3,4,5) | 1.91 | 14% |
    | `SolverT<N>::LogPVec` compute | 1.24 | 9% |
    | `LogPVecKey::operator==` | 1.11 | 8% |
    | `StateKeyT<N>::operator==` | 1.01 | 7% |
    | `expl`/`_expl_internal` (LogB LSE) | 0.91 | 6% |
    | `logRD_cache::find` | 0.74 | 5% |
    | `std::isfinite` (all sites) | 0.70 | 5% |
    | `vector<double>::~vector` (eviction) | 0.60 | 4% |
    | `logconv` actual convolution | 0.20 | 1% |

    **Key findings:**
    - `logconv` is only **1%** of DLL time — the Phase 2 vectorization worked
      perfectly; the algorithm itself is no longer the bottleneck.
    - **Hash map infrastructure dominates** (53% of DLL time): `unordered_map::find`
      + key equality checks across the three caches (logB, logPVec, logRD).
      Switching to a flat/open-addressing map would help but adds complexity.
    - **`expl()` in `LSEAccumulator`** (6%) uses long-double arithmetic. Switching
      to `double`/`exp()` would save ~0.7s at negligible precision cost. → **T-151**
    - **`std::isfinite`** (5%) routes through `_fpclassify` on MinGW/Windows.
      Replacing with `x != NEG_INF` saves the function-call overhead. → **T-152**
    - `memcmp` in ucrtbase.dll (1.6 s / 7% of total) is the `StateKeyT::operator==`
      fall-through when `cached_hash` and `cached_sum` both match — unavoidable
      with the current key design.

    **Estimated combined T-151 + T-152 saving: ~1.4 s (6%) per cold-cache run.**

## Build and Test (Reminder)

Always use isolated library:
```bash
R CMD build --no-build-vignettes --no-manual . && R CMD INSTALL --library=.agent-X TreeSearch_*.tar.gz && rm -f TreeSearch_*.tar.gz
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='ts-')"
```

Max 2 CPU cores. Use `nThreads = 2L` at most in benchmarks.
