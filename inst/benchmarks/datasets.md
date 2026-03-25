# Benchmark Dataset Suite

Selected from the 30 `inapplicable.phyData` datasets bundled with TreeSearch.
Criteria: cover small → large tip counts, varying inapplicable proportions,
varying state counts, and varying matrix densities (% missing data).

## Dataset Selection

| # | Dataset | Tips | Chars | Patterns | %Inapp | States | %Missing | Category |
|---|---------|------|-------|----------|--------|--------|----------|----------|
| 1 | Longrich2010 | 20 | 93 | 80 | 4.2 | 3 | 45.3 | Small, high missing |
| 2 | Vinther2008 | 23 | 57 | 50 | 6.1 | 4 | 21.0 | Small, moderate |
| 3 | Sansom2010 | 23 | 109 | 97 | 6.1 | 4 | 40.0 | Small, high missing |
| 4 | DeAssis2011 | 33 | 50 | 36 | 21.4 | 3 | 0.2 | Medium-small, high inapp |
| 5 | Aria2015 | 35 | 50 | 50 | 6.7 | 6 | 12.7 | Medium-small, multi-state |
| 6 | Wortley2006 | 37 | 105 | 105 | 2.7 | 8 | 31.4 | Medium, many states |
| 7 | Griswold1999 | 43 | 137 | 118 | 6.2 | 6 | 5.6 | Medium, dense matrix |
| 8 | Schulze2007 | 52 | 58 | 57 | 16.7 | 3 | 2.4 | Medium, high inapp, dense |
| 9 | Eklund2004 | 54 | 131 | 131 | 7.8 | 6 | 29.8 | Medium, moderate |
| 10 | Agnarsson2004 | 62 | 242 | 225 | 6.9 | 7 | 6.1 | Large, many chars, dense |
| 11 | Zanol2014 | 74 | 213 | 210 | 16.8 | 9 | 11.9 | Large, high inapp, many states |
| 12 | Zhu2013 | 75 | 253 | 253 | 12.4 | 4 | 42.6 | Large, high missing |
| 13 | Giles2015 | 78 | 236 | 236 | 11.8 | 4 | 41.5 | Large, high missing+inapp |
| 14 | Dikow2009 | 88 | 220 | 204 | 1.2 | 9 | 0.4 | Largest, dense, many states |

## Selection Rationale

- **Size range**: 20 → 88 tips (5× range). Covers small (exhaustive-feasible)
  through large (heuristic-only).
- **Inapplicable variation**: 1.2% (Dikow) → 21.4% (DeAssis). Tests the
  NA three-pass scoring path under varying load.
- **State count variation**: 3–9 applicable states. Affects `total_words`
  (state word count per block) and thus inner-loop iteration count.
- **Missing data variation**: 0.2% (DeAssis) → 45.3% (Longrich). High missing
  data creates more ambiguous tokens, affecting scoring and simplification.
- **Dense vs sparse**: DeAssis (0.2% missing) and Dikow (0.4% missing) are
  nearly complete matrices; Longrich (45.3%) and Zhu (42.6%) are sparse.

## Best-Known EW Scores (standard Fitch, inapplicable → missing)

Best scores from TreeSearch (cpp-search, tuned "default" strategy) and
TNT 1.6, across benchmark rounds 1 and 2 (2026-03-25). All inapplicable
tokens replaced with missing data so both engines optimize the same Fitch
objective. TNT scores are authoritative lower bounds.

| Dataset | TNT Best | TS Best | Δ | Published Tree | Notes |
|---------|----------|---------|---|---------------|-------|
| Longrich2010 | 131 | 131 | 0 | 167 | |
| Vinther2008 | 78 | 78 | 0 | 93 | |
| Sansom2010 | 188 | 188 | 0 | — | |
| DeAssis2011 | 64 | 64 | 0 | 89 | |
| Aria2015 | 142 | 142 | 0 | 185 | |
| Wortley2006 | 479 | 482 | +3 | 518 | |
| Griswold1999 | 394 | 394 | 0 | 511 | |
| Schulze2007 | 155 | 155 | 0 | 212 | |
| Eklund2004 | 440 | 440 | 0 | 496 | |
| Agnarsson2004 | 765 | 765 | 0 | 1035 | |
| Shultz2007 | 431 | 431 | 0 | — | Round 2 only |
| Rougier2012 | 1147 | 1148 | +1 | — | Round 2 only |
| Wilson2003 | 860 | 860 | 0 | — | Round 2 only |
| OMeara2014 | 1208 | 1209 | +1 | — | Round 2 only |
| Wetterer2000 | 549 | 549 | 0 | — | Round 2 only |
| Conrad2008 | 1725 | 1726 | +1 | — | Round 2 only |
| Capa2011 | 381 | 381 | 0 | — | Round 2 only |
| Geisler2001 | 1293 | 1298 | +5 | — | Round 2 only |
| Liljeblad2008 | 2840 | 2840 | 0 | — | Round 2 only |
| Zanol2014 | 1261 | 1263 | +2 | 1802 | |
| Zhu2013 | 624 | 627 | +3 | 810 | |
| Aguado2009 | 575 | 575 | 0 | — | Round 2 only |
| Giles2015 | 670 | 671 | +1 | 1005 | |
| Dikow2009 | 1606 | 1606 | 0 | 1646 | |

Note: Published tree scores are higher because they were optimized for
inapplicable-aware scoring (Brazeau et al. 2019), not standard Fitch.

## Usage

### Running benchmarks

Benchmarks require a built copy of TreeSearch from `TreeSearch-a` (cpp-search):

```bash
# From the GitHub/ root
bash build-agent.sh TreeSearch-a F          # builds TreeSearch.F
cd TS-TNT-bench/inst/benchmarks
Rscript --no-save bench_round1_cppsearch.R  # 14 datasets, 10s timeout
Rscript --no-save bench_round2_hard.R       # 16 datasets, 120s timeout
```

Both scripts use tuned "default" strategy parameters from `MaximizeParsimony`
(25% ratchet perturbation, adaptive level, outer-cycle resets). Results are
saved to `round1_cppsearch.csv` and `round2_hard.csv`.

**Important:** Do not use the stale compat-layer defaults (`ratchetPerturbProb
= 0.04`, `driftCycles = 6`). These were the old `ts_driven_search()` defaults
and produce substantially worse scores. Always use the tuned parameters from
`MaximizeParsimony`'s "default" strategy.

### TNT requirement

TNT 1.6 must be installed at `C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe`.
See `.positai/expertise/tnt.md` for invocation notes.

### Dataset utilities

```r
source("inst/benchmarks/bench_datasets.R")
datasets <- load_benchmark_datasets()
score_dataset("Vinther2008", maxSeconds = 10)
run_benchmark_suite(maxSeconds = 30, replicates = 5)
```
