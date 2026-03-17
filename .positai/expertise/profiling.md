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

## Known Baselines (as of 2026-03-17)

### Per-phase breakdown (Zhu2013 EW, 75 tips, 3 replicates):

| Phase | Avg ms/rep | % of C++ time |
|-------|-----------|---------------|
| Drift | 365 | 39% |
| Ratchet | 274 | 29% |
| Initial TBR | 210 | 23% |
| Sectorial (XSS+RSS+CSS) | 75 | 8% |
| Final TBR + Wagner | 9 | 1% |

### End-to-end benchmarks (3-run medians, post R-overhead fix):

| Dataset | Tips | EW (s) | IW (s) |
|---------|------|--------|--------|
| Vinther2008 | 23 | 0.160 | 0.180 |
| Zhu2013 | 75 | 2.730 | 3.560 |
| Agnarsson2004 | 62 | 1.520 | — |

### Per-candidate indirect scoring:

| Tips | Words | ns/candidate |
|------|-------|-------------|
| 20 | 12 | ~33 |
| 75 | 20 | ~46 |
| 62 | 59 | ~110 |

Per-candidate cost scales linearly with `total_words`. Near memory-throughput
limit — hard to improve without algorithmic changes.

### Scaling (TBR pass time vs tips):

Exponent ~2.78 (vs expected 2.0). Super-quadratic growth is from candidate
count growth (exponent 2.66), not per-candidate degradation.

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
