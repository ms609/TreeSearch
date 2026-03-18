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

## Known Baselines (as of 2026-03-17 18:30, v2.0.0 verification run by Agent A)

### Per-phase breakdown (EW, strategy='none', 5 replicates, verbose run):

| Dataset | Tips | TBR% | XSS% | RSS% | CSS% | Ratch% | Drift% | Med ms |
|---------|------|------|------|------|------|--------|--------|--------|
| Vinther2008 | 23 | 13.6 | 18.6 | 4.5 | 0.0 | 41.2 | 22.1 | 550* |
| Agnarsson2004 | 62 | 21.5 | 14.8 | 2.8 | 0.0 | 38.7 | 22.3 | 3420 |
| Zhu2013 | 75 | 39.6 | 10.1 | 2.1 | 0.0 | 24.1 | 24.2 | 4930 |
| Dikow2009 | 88 | 28.8 | 11.2 | 2.7 | 0.0 | 35.5 | 21.7 | 6490 |

*Vinther2008 median inflated by MPT enumeration (49 trees in pool).

Note: Current defaults are driftCycles=2, ratchetCycles=5, cssRounds=0.
Phase distribution stable vs prior baselines (stochastic variation only).
Ratchet (24-41%) and TBR (14-40%) remain the largest phases.

### End-to-end benchmarks (3-run medians, 5 reps, strategy='none'):

| Dataset | Tips | EW (s) | IW k=10 (s) |
|---------|------|--------|-------------|
| Vinther2008 | 23 | 0.550* | 0.170 |
| Agnarsson2004 | 62 | 3.420 | 2.890 |
| Zhu2013 | 75 | 4.930 | 4.610 |
| Dikow2009 | 88 | 6.490 | — |

IW consistently faster than EW (fewer MPTs → smaller pool → less enumeration).

### Auto strategy: default vs thorough comparison (8 seeds, 5 reps each)

| Dataset | Tips | default med | thorough med | Slowdown | Score improvement |
|---------|------|-------------|-------------|----------|-------------------|
| Wilson2003 | 61 | 889 (3.8s) | 887 (10.2s) | 2.7× | 2 (noisy) |
| Agnarsson2004 | 62 | 778 (2.2s) | 778 (5.4s) | 2.5× | 0 |
| Zhu2013 | 75 | 658 (4.6s) | 649 (12.2s) | 2.6× | **9** |
| Dikow2009 | 88 | 1614 (5.4s) | 1612 (17.9s) | 3.3× | 2 |

**Finding:** "thorough" benefit is dataset-dependent, not purely size-dependent.
At 61-62 tips: 0-2 step improvement (noisy). At 75 tips (Zhu2013): 9-step
improvement justifies the cost. At 88 tips: only 2-step improvement.

**Recommendation:** Consider raising the threshold from 61 to ~75 tips, or
introducing an intermediate preset (e.g. 10 ratchet + 5 drift) for 61-75.
The current threshold wastes 2.5-3× compute on many 61-tip datasets.

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
