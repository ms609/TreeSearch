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

### 4. R-Level Profiling

For R overhead identification:

```r
Rprof("profile.out")
result <- MaximizeParsimony(dataset, maxReplicates = 5)
Rprof(NULL)
summaryRprof("profile.out")
```

## Known Baselines

### Latest run: 2026-03-18 16:00 by Agent A (v2.0.0, single-agent, quiet machine)

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

## What to Profile

Status key: ✅ resolved, ⚠ partially explored, ❌ not yet investigated

1. ⚠ **Drift + ratchet inner loops** (50–60% of C++ time combined). Both use
   TBR internally. Per-candidate indirect evaluation at memory-throughput
   limit (~23 ns at 75 tips per T-075). Cycle counts tuned (d2_r5).
   Remaining question: is drift acceptance threshold optimal?

2. ✅ **Sectorial search effectiveness** (12% of time). XSS effectiveness is
   dataset-dependent (0–27 steps). RSS is marginal (0–3 steps). No clear
   predictor from simple dataset statistics. Could make XSS adaptive (skip
   after N unproductive reps) but time savings would be <10%.

3. ✅ **Wagner tree construction**: <0.1% of search time. Not a bottleneck.

4. ✅ **R overhead**: <0.5% of wall time. Not a bottleneck.

5. ✅ **Parallel scaling**: 78–82% efficiency at 2 threads. Implementation is
   sound (dynamic work-stealing, low-contention pool). Main loss is stochastic
   load imbalance. No obvious improvement without algorithmic changes.

## Reporting Format

For each finding, add to `to-do.md`:

```
| T-NNN | P2 | OPEN | — | [Profile] Brief description | X% of time. Potential Y% improvement via Z approach. |
```

Include the measurement methodology and baseline numbers so the implementer
can verify the improvement.

## Build and Test (Reminder)

Always use isolated library:
```bash
R CMD INSTALL --library=.agent-X .
Rscript -e "library(TreeSearch, lib.loc='.agent-X'); testthat::test_dir('tests/testthat', filter='ts-')"
```

Max 2 CPU cores. Use `nThreads = 2L` at most in benchmarks.
