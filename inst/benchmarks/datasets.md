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

## Best-Known EW Scores

Scores from the C++ driven search engine (5 replicates, 5s timeout per
dataset, `set.seed(42)`). These are the standard Fitch parsimony scores
(not inapplicable-aware). Published tree scores from `inapplicable.trees`
are generally higher because they may not be optimized for standard Fitch.

| Dataset | C++ Best | Published Tree | Notes |
|---------|----------|---------------|-------|
| Longrich2010 | 131 | 167 | |
| Vinther2008 | 79 | 93 | |
| Sansom2010 | 189 | — | |
| DeAssis2011 | 64 | 89 | |
| Aria2015 | 145 | 185 | |
| Wortley2006 | 496 | 518 | |
| Griswold1999 | 409 | 511 | |
| Schulze2007 | 167 | 212 | |
| Eklund2004 | 445 | 496 | |
| Agnarsson2004 | 778 | 1035 | |
| Zanol2014 | 1338 | 1802 | |
| Zhu2013 | 649 | 810 | |
| Giles2015 | 720 | 1005 | |
| Dikow2009 | 1614 | 1646 | |

Note: C++ scores are lower than published because (a) the published trees
were optimized for a different scoring method (inapplicable-aware), and
(b) our driven search may find better trees. These scores were obtained
with `set.seed(42)`, 10s timeout, 10 replicates. Use `bench_datasets.R`
with longer search times for authoritative best-known scores.

## Large-Tree Benchmark Datasets

Separate tier for datasets >= 100 tips, loaded from `inst/benchmarks/`.
These have fundamentally different search dynamics: single TBR convergence
takes seconds to minutes, replicates take minutes rather than sub-second.

| # | Dataset | Tips | Chars | Patterns | %Missing | %Inapp | Source |
|---|---------|------|-------|----------|----------|--------|--------|
| L1 | mbank_X30754 | 180 | 425 | 418 | 40% | 20.5% | MorphoBank P30754 |

### mbank_X30754

MorphoBank project X30754 (downloaded 2025-06-16). 180 taxa, 425 characters
with ~40% missing data and ~20% inapplicable entries. This is a realistic
large morphological matrix that exposes scaling issues in the search engine:
NNI warmup is essential, single TBR convergence takes ~13s, and the standard
strategy presets (calibrated for ≤88 tips) are poorly suited.

Best-known EW score: TBD (to be established after systematic benchmarking).

## Usage

```r
source("inst/benchmarks/bench_datasets.R")

# Load standard benchmark datasets (14 datasets, ≤88 tips)
datasets <- load_benchmark_datasets()

# Load large-tree benchmark datasets (≥100 tips)
large <- load_large_benchmark_datasets()

# Load all (standard + large)
all_ds <- load_all_benchmark_datasets()

# Score a single dataset
score_dataset("Vinther2008", maxSeconds = 10)

# Run standard benchmark suite
run_benchmark_suite(maxSeconds = 30, replicates = 5)

# Run large-tree benchmark (from bench_framework.R)
# source("inst/benchmarks/bench_framework.R")
# benchmark_large(maxSeconds = 120)
```
