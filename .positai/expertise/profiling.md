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

## Known Baselines (as of 2026-03-17, post T-025 fix + d2_r5 tuning + CSS disabled)

### Per-phase breakdown (EW, strategy='none', 5 replicates, 3-run medians):

| Dataset | Tips | TBR% | XSS% | RSS% | CSS% | Ratch% | Drift% | Other% | Total ms |
|---------|------|------|------|------|------|--------|--------|--------|----------|
| Vinther2008 | 23 | 11.0 | 21.9 | 4.3 | 0.0 | 39.8 | 19.6 | 3.4 | 233 |
| Agnarsson2004 | 62 | 18.0 | 13.5 | 4.2 | 0.0 | 37.9 | 24.0 | 2.4 | 3141 |
| Zhu2013 | 75 | 33.2 | 9.8 | 2.8 | 0.0 | 24.7 | 27.8 | 1.7 | 3995 |
| Dikow2009 | 88 | 20.5 | 12.4 | 2.6 | 0.0 | 37.5 | 24.7 | 2.3 | 6666 |

Note: Current defaults are driftCycles=2, ratchetCycles=5, cssRounds=0.
Previous baselines (d6_r10) showed drift dominating at 40-50%; with d2_r5,
ratchet (25-40%) and TBR (11-33%) are now the largest phases.

### End-to-end benchmarks (3-run medians, 5 reps, strategy='none'):

| Dataset | Tips | EW (s) | IW k=10 (s) |
|---------|------|--------|-------------|
| Vinther2008 | 23 | 0.230 | 0.360 |
| Agnarsson2004 | 62 | 3.250 | 4.730 |
| Zhu2013 | 75 | 4.080 | 5.920 |
| Dikow2009 | 88 | 6.630 | — |

Speedup vs old defaults (d6_r10): 41-52% faster with equivalent score quality.

### Auto strategy benchmarks (3-run medians, 5 reps):

| Dataset | Tips | Preset | EW (s) | Score |
|---------|------|--------|--------|-------|
| Vinther2008 | 23 | sprint | 0.110 | 79-80 |
| Agnarsson2004 | 62 | thorough | 13.670 | 778 |
| Zhu2013 | 75 | thorough | 16.390 | 647-648 |

Note: "thorough" preset is 4× slower than raw defaults for 62-tip datasets
with zero score improvement. Threshold at 61 tips may be too aggressive.

### R overhead: <0.5% of wall time (confirmed via Rprof)

### Parallel scaling (Zhu2013, 5 reps):
- 2 threads: 1.24× speedup (62% efficiency)
- Note: degraded from previous 1.86× due to shorter per-replicate time with
  d2_r5 defaults. Thread overhead is proportionally larger.

### Scaling (TBR pass time vs tips):

Exponent ~2.82 (23→88 tips). Consistent with prior measurement of 2.78.
Super-quadratic growth is from candidate count growth, not per-candidate
degradation.

### Drift/ratchet cycle tuning (Zhu2013, 10 seeds, 5 reps):

| Config | Med score | Min score | Med time | Speedup |
|--------|-----------|-----------|----------|---------|
| d5_r5 (default) | 656 | 648 | 5.7s | — |
| d2_r5 | 660 | 646 | 4.1s | 28% |
| d2_r2 | 662 | 656 | 3.8s | 33% |
| d0_r5 | 658 | 650 | 2.8s | 51% |
| d5_r0 | 662 | 660 | 4.8s | 16% |

Lower score = better. Drift has diminishing returns beyond 2 cycles.
Ratchet 5 is more valuable than drift 5 for quality.

### CSS effectiveness: Marginal (adds 2-6% time, no consistent score improvement)

## What to Profile

In order of likely impact:

1. **Drift + ratchet inner loops** (68% of C++ time combined). Both use TBR
   internally. The per-candidate indirect evaluation is near-optimal. Look for:
   - Number of TBR passes per drift/ratchet cycle (are we doing too many?)
   - Cycle counts (are defaults too high for the dataset?)
   - Acceptance criteria (is drift accepting too many/few suboptimal moves?)

2. **Sectorial search effectiveness** (8% of time). Questions:
   - Is XSS finding improvements? (Count accepted vs rejected sectors)
   - Are CSS rounds redundant after XSS+RSS?
   - Sector size distribution — are we getting good partitions?

3. **Wagner tree construction** — now O(n × depth × C) per insertion.
   Profile whether the incremental scoring early termination is triggering
   effectively.

4. **R overhead** — the `AdditionTree()` → `RandomTree()` fix eliminated the
   biggest R bottleneck. Look for remaining overhead in:
   - `PrepareDataIW()` / `PrepareDataProfile()`
   - `Renumber()` on output trees
   - Pool → multiPhylo conversion

5. **Parallel scaling** — with `nThreads > 1`, measure:
   - Thread spawn/join overhead
   - Pool mutex contention (should be negligible)
   - Per-thread memory allocation cost

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
