# Benchmarks and Profiling

Load this when: running benchmarks, interpreting benchmark results,
doing VTune profiling, or selecting datasets for strategy validation.

See also: `search-algorithms.md` (NNI, biased Wagner, outer cycles results),
`search_strategy.md` (presets, ratchet tuning).

---

## VTune driver scripts — dry-run first

**Always test a VTune driver script with plain `Rscript` before launching
VTune.** Software-sampling overhead can be 5–20×; if the bare script takes
30s, VTune may need 10 min. Target < 5s bare run for a lite driver.

MaddisonSlatkin is exponential in tip count — even n=20 with k=3 can take
seconds per call. Use small n (≤15 for k=3, ≤12 for k=4, ≤9 for k=5)
and few iterations for VTune drivers.

---

## MorphoBank external benchmark corpus

The neotrans repo (`../neotrans/inst/matrices/`) contains ~800 MorphoBank
NEXUS matrices. Complement to the 14 bundled datasets and 1 large-tree dataset.

**Catalogue:** `dev/benchmarks/mbank_catalogue.csv` (659 usable matrices
after ntax≥20 filter and dedup). Regenerate with
`Rscript dev/benchmarks/build_mbank_catalogue.R`.

**Train/validation split:** Matrices whose MorphoBank project number is
divisible by 5 are **validation** (124 matrices, ~19%). All others are
**training** (535 matrices). The 7 `syab*` files are always training.

**Dedup:** Multi-file projects with ≥95% character identity on shared taxa
(≥80% taxon overlap) are flagged `dedup_drop = TRUE`. 24 near-duplicates excluded.

**IMPORTANT:** Validation results must **never** be used to guide strategy
tuning. They confirm generalization only. This is a one-way door.

**Fixed 25-matrix training sample:** `MBANK_FIXED_SAMPLE` in
`bench_datasets.R` — 7 small, 7 medium, 7 large, 4 xlarge. Selected via
max-min distance on standardized features. **Do not modify.** Used by
`benchmark_mbank_sample()`. Fitch track only.

**Fixed 20-matrix Brazeau-track sample:** `MBANK_BRAZEAU_SAMPLE` in
`bench_datasets.R` — 5 small, 6 medium, 6 large, 3 xlarge. Restricted to
training matrices with **pct_inapp ≥ 4%**. **Do not modify.**

**Key functions** (in `dev/benchmarks/bench_datasets.R`):
- `load_mbank_catalogue()` — loads metadata CSV (excludes dedup by default)
- `load_mbank_sample(cat, n, seed, split)` — stratified random sample
- `load_mbank_datasets(cat, keys)` — load specific matrices by key
- `load_mbank_brazeau_sample(cat)` — fixed 20-matrix Brazeau sample
- `has_meaningful_inapp(cat, threshold)` — filter to pct_inapp ≥ threshold

**Benchmark runners** (in `dev/benchmarks/bench_framework.R`):
- `benchmark_mbank_sample()` — fixed 25-matrix training sample (routine)
- `benchmark_mbank_sweep(split)` — full training or validation sweep
- `benchmark_mbank_validation()` — validation sweep with prominent warning

**Benchmark tracks:**

| Track | Scoring | Datasets | Purpose |
|-------|---------|----------|---------|
| **Fitch** | `fitch_mode()` | 14 bundled + `MBANK_FIXED_SAMPLE` | TNT comparison, core search quality |
| **Brazeau** | Default (Brazeau 2019) | `MBANK_BRAZEAU_SAMPLE` + bundled | NA-algorithm-specific strategy tuning |

TNT comparisons are Fitch track only.

**TNT comparison suite** lives in `../TS-TNT-bench/`. Key files:
- `dev/benchmarks/bench_tnt_compare.R` — runner (smoke/medium/full)
- `dev/benchmarks/tnt_comparison.qmd` — Quarto report
- Requires TNT 1.6 at `C:/Programs/Phylogeny/tnt/TNT-bin/tnt.exe`

Benchmark scripts in `dev/benchmarks/`. Key files:
- `bench_regression.R` — CI regression test (score quality + timing bounds)
- `bench_framework.R` — Dataset × strategy × replicate grid
- `strategies.md` — Strategy space documentation

---

## Benchmarking methodology notes

**Metric:** When comparing strategies with different time costs (e.g.
NNI→TBR vs TBR), use **time-adjusted expected best** (TAEB) — the expected
minimum score from k = budget / time_per_rep independent replicates. Median
per-replicate score is adequate only when comparing parameter changes on a
fixed pipeline (same time-per-rep). Bootstrap estimation: sample k scores
with replacement, take the min, repeat 5000×, take the mean.

**Brazeau vs EW scoring confound (T-265, 2026-03-26):** TreeSearch uses the
Brazeau et al. (2019) inapplicable algorithm by default, which penalizes
inapplicable-to-applicable transitions. TNT treats `-` as `?` (standard EW
Fitch). On 11 gap datasets, the apparent mean gap was +17.8 steps; the
actual EW-vs-EW gap is only +2.2 steps (5 datasets at 0 gap). **All TNT
comparisons MUST use `fitch_mode()` to convert inapplicable to missing**
for apples-to-apples scoring. `fitch_mode()` is defined in
`bench_intra_fuse.R` and `bench_t265_regression.R`.

**Early vs late search:** Early replicates are dominated by initial descent
quality (Wagner → local optimum); late replicates test ratchet/drift escape.
At ≤88 tips, 20s gives 10–40 replicates spanning both regimes. At 180 tips,
20s doesn't complete one replicate.

---

## Phase distribution baselines

**T-290b (2026-03-28, Brazeau-sample datasets, 30s, post-T-255 no-drift presets):**

| Phase | Fitch/EW/default | Fitch/EW/thorough | Brazeau/EW/default | Brazeau/EW/thorough |
|-------|:---:|:---:|:---:|:---:|
| Ratchet | 76% | 65% | 74% | 63% |
| TBR | 8% | 5% | 7% | 4% |
| XSS | 6% | 7% | 5% | 6% |
| RSS | 3% | 10% | 3% | 10% |
| CSS | — | 7% | — | 7% |
| Wagner | 4% | 3% | 9% | 7% |
| Final TBR | 2% | 2% | 2% | 2% |

*(Drift has been 0% in all presets since T-255.)*

**Brazeau / Fitch per-phase cost ratios (T-290b, EW):**

| Phase | default | thorough |
|-------|:-------:|:--------:|
| Wagner | **3.6×** | **3.9×** |
| Ratchet | 1.3× | 1.3× |
| RSS/CSS | 1.3× | 1.3× |
| TBR | 0.9× | 0.9× |

Wagner is the outlier. All other phases are within 0.9–1.4× of Fitch cost.

**wagnerStarts under Brazeau (T-290b/c, 2026-03-28):**
- *Multiple reps/budget*: wagnerStarts=1 and 3 equivalent; w3 marginally better.
- *~1 rep/budget* (60s at 86t/3660c): wagnerStarts=3 better by +564 steps.
- *0 reps/budget* (30s at 86t/3660c): wagnerStarts=1 **better** — Brazeau
  Wagner is expensive (~4×), 3 starts consume budget.
Current presets correct: thorough (w3, gets ≥1 rep at 65–119t) ✓; large (w1) ✓.

Per-candidate indirect scoring is at memory-throughput limit (~23 ns at 75 tips).

---

## Ratchet tuning validation (2026-03-22)

Full 14-dataset comparison, optimized vs original defaults (10s budget, 3 seeds).

| Dataset | Tips | Original | Optimized | Delta |
|---------|:---:|:---:|:---:|:---:|
| Longrich2010 | 20 | 131 | 131 | 0 |
| Vinther2008 | 23 | 79 | 79 | 0 |
| Sansom2010 | 23 | 189 | 189 | 0 |
| DeAssis2011 | 33 | 64 | 64 | 0 |
| Aria2015 | 35 | 143 | 143 | 0 |
| Wortley2006 | 37 | 494 | 491 | +3 |
| Griswold1999 | 43 | 408 | 407 | +1 |
| Schulze2007 | 52 | 165 | 164 | +1 |
| Eklund2004 | 54 | 442 | 441 | +1 |
| Agnarsson2004 | 62 | 778 | 778 | 0 |
| Zanol2014 | 74 | 1338 | 1331 | +7 |
| Zhu2013 | 75 | 649 | 650 | −1 |
| Giles2015 | 78 | 720 | 716 | +4 |
| Dikow2009 | 88 | 1614 | 1614 | 0 |

Zhu2013 marginal regression at 10s resolves at 20s (median 649→644).
At 20s with 5 seeds: Zhu2013 645/643, Giles2015 712/710, Dikow2009
1611/1611 (all improvements).
