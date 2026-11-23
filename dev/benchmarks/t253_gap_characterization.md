# T-253: Gap Characterization by Dataset Features

**Date:** 2026-03-27  
**Agent:** F  
**Data sources:**
- `t265_results/t265_phase1_20260326_1617.csv` — 8 named datasets, fitch_mode EW, 120s (TNT vs TreeSearch, apples-to-apples)
- `t252_mbank_*` CSVs — 25 MorphoBank training matrices, TreeSearch 30/60/120s (convergence proxy)

---

## Summary

**ntax is the primary predictor of search difficulty** in both analyses (Spearman ρ ≈ 0.63).
At ≤60 taxa with modest character counts, TreeSearch converges fully at 30s.
Difficulty increases steadily above ~75 taxa and becomes acute above ~120 taxa.

Character count (nchar) matters only at extremes (e.g. 3660 chars, 2954 chars);
pct_missing and pct_inapp show moderate individual correlations (ρ = 0.49–0.55 in T-265)
but inconsistent signal in the MorphoBank sample — small samples mean these
correlations are unreliable beyond the ntax signal.

---

## TNT comparison gaps (T-265, fitch_mode, 120s, 8 datasets)

These are the only reliable apples-to-apples gaps (Fitch TreeSearch vs TNT Fitch).

| Dataset | ntax | nchar | pct_missing | pct_inapp | median_gap |
|---------|-----:|------:|:-----------:|:---------:|:----------:|
| Zanol2014 | 74 | 213 | 11.7% | 16.6% | **3** |
| Zhu2013 | 75 | 253 | 42.6% | 12.4% | **3** |
| Conrad2008 | 64 | 363 | 23.4% | 5.1% | 2 |
| Giles2015 | 78 | 236 | 41.5% | 11.8% | 2 |
| OMeara2014 | 63 | 317 | 43.4% | 5.4% | 2 |
| Liljeblad2008 | 68 | 308 | 5.2% | 5.6% | 0 |
| Wetterer2000 | 63 | 150 | 21.2% | 7.7% | 0 |
| Wilson2003 | 61 | 165 | 7.7% | 8.6% | 0 |

Spearman correlations with `median_gap`:

| Feature | ρ |
|---------|:-:|
| ntax | 0.63 |
| pct_missing | 0.55 |
| pct_inapp | 0.49 |
| nchar | 0.28 |
| n_patterns | 0.28 |

**Note:** n=8 is too small for reliable multivariate analysis. The pct_missing/pct_inapp
signals may be confounded with ntax (larger datasets often have more missing data).

---

## Convergence gaps (T-252, MorphoBank 25 matrices, 30s → 120s improvement)

Most matrices converge fully at 30s (gap=0). Non-zero gap datasets:

| Dataset | ntax | nchar | pct_miss | pct_inapp | conv_gap |
|---------|-----:|------:|:--------:|:---------:|:--------:|
| project2068 | 86 | **3660** | 20.9% | 24.9% | 238 |
| project4284 | **4062** | 27 | 82.9% | 0% | 75 |
| syab072 | 125 | 2954 | 28.3% | ? | 21 |
| project804 | 173 | 589 | 32.8% | ? | 10 |
| project3938 | 119 | 677 | 52.6% | 4.3% | 8 |
| project2771 | 94 | 124 | 1.0% | 30.0% | 6 |
| (others) | ≤131 | ≤721 | | | ≤2 |

Spearman correlations with `conv_gap` (n=23, excluding 2 extreme outliers):

| Feature | ρ |
|---------|:-:|
| ntax | **0.64** |
| n_patterns | 0.34 |
| pct_inapp | 0.36 |
| nchar | 0.30 |
| pct_missing | −0.04 |

---

## Key findings

1. **ntax is the dominant difficulty predictor** (consistent ρ ≈ 0.63 across two
   independent datasets/metrics). The hard wall is around 75–130 taxa under the
   current strategy presets.

2. **nchar matters only at extremes.** project2068 (86t, 3660c) has the largest
   absolute convergence gap despite modest ntax — the 3660-character search space
   is simply too large per-replicate. syab072 (125t, 2954c) similarly.

3. **Missing data and inapplicable characters** show moderate correlations in T-265
   but not in T-252. This likely reflects a confound with ntax (larger datasets often
   have more missing data in MorphoBank matrices), not an independent effect.

4. **Most datasets are already covered** (≤60 taxa, ≤700 chars): 19 of 25 MorphoBank
   training matrices and all datasets ≤60 taxa converge at 30s. TreeSearch's
   CRAN benchmark suite (14 datasets, ≤88 taxa) is well-covered.

---

## Strategic implications for T-253

| Priority | Action | Targets |
|----------|--------|---------|
| High | **T-245: TBR batching** — reduce per-candidate evaluation cost | ≥75 taxa (nchar moderate) |
| High | **NNI escalation** (already in presets via `nniFirst=TRUE`) | ≥75 taxa |
| Medium | **Character batching / lazy scoring** for high-nchar datasets | ≥1000 chars |
| Low | Missing/inapplicable tuning | Not independently predictive |

The clearest opportunity is the ≥75-taxon regime. T-245 (TBR candidate batching,
estimated ~13% gain) is the highest-value next step for search quality at scale.
