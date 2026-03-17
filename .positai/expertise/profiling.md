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

## Known Baselines (as of 2026-03-17, post all Phase 2-3 optimizations)

### Per-phase breakdown (EW, default params, 5 replicates, 3-run medians):

| Dataset | Tips | Wagner% | TBR% | Sect% | Ratch% | Drift% | Other% | Total ms |
|---------|------|---------|------|-------|--------|--------|--------|----------|
| Vinther2008 | 23 | 0.2 | 9.4 | 14.0 | 33.5 | 39.8 | 3.2 | 255 |
| Agnarsson2004 | 62 | 0.1 | 14.4 | 10.2 | 27.9 | 45.6 | 1.9 | 2482 |
| Zhu2013 | 75 | 0.0 | 25.0 | 6.8 | 17.8 | 49.4 | 1.1 | 5100 |
| Dikow2009 | 88 | 0.1 | 16.5 | 8.7 | 27.4 | 45.7 | 1.7 | 8141 |

### Per-replicate cost (ms/rep, default params):

| Dataset | Tips | Wagner | TBR | Sect | Ratch | Drift | Total |
|---------|------|--------|-----|------|-------|-------|-------|
| Vinther2008 | 23 | 0.1 | 9.3 | 12.9 | 33.1 | 39.6 | 98 |
| Agnarsson2004 | 62 | 0.8 | 127 | 96 | 238 | 472 | 957 |
| Zhu2013 | 75 | 0.5 | 275 | 69 | 180 | 503 | 1040 |
| Dikow2009 | 88 | 0.9 | 292 | 147 | 449 | 785 | 1703 |

### End-to-end benchmarks (3-run medians, 5 reps, default params):

| Dataset | Tips | EW (s) | IW k=10 (s) |
|---------|------|--------|-------------|
| Vinther2008 | 23 | 0.250 | 0.310 |
| Agnarsson2004 | 62 | 2.510 | — |
| Zhu2013 | 75 | 5.320 | 7.930 |
| Dikow2009 | 88 | 8.270 | — |

### R overhead: <0.5% of wall time (confirmed via Rprof)

### Parallel scaling:
- 2 threads: 1.86× speedup (93% efficiency) on Zhu2013

### Wagner + TBR only (single replicate, no perturbation):

| Dataset | Tips | Wagner (ms) | TBR (ms) |
|---------|------|-------------|----------|
| Vinther2008 | 23 | 0.1 | 6.8 |
| Agnarsson2004 | 62 | 0.7 | 92.0 |
| Zhu2013 | 75 | 0.3 | 226.9 |
| Dikow2009 | 88 | 1.0 | 300.3 |

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
