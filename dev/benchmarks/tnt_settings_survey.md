# TNT 1.6 Settings Survey — Time-to-Target (TTT)

**Metric:** wall-clock seconds to first reach the best known score B for each dataset  
**Censored:** runs that never reached B within 120 s are marked NA (3 of 420 runs)

**Machine:** DW-CZC429715G · 12th Gen Intel Core i7-12700 · 15.7 GB RAM  
**TNT:** C:/Programs/Phylogeny/tnt/tnt.exe (v1.6, 32-bit)  
**Date:** 2026-06-17  
**Script:** `dev/benchmarks/bench_tnt_settings.R`  
**Data:** `dev/benchmarks/tnt_settings_survey.csv`

---

## Dataset reference

| Dataset     | Tips | B (best score) |
|-------------|-----:|---------------:|
| Wortley2006 |   37 |            479 |
| Eklund2004  |   54 |            440 |
| Zanol2014   |   74 |          1 261 |
| Zhu2013     |   75 |            624 |
| Giles2015   |   78 |            670 |
| Dikow2009   |   88 |          1 606 |

All datasets run in **Fitch mode** (inapplicable `-` tokens replaced with `?`).  
B was established by 10 TNT seeds (Phase 1) with 300 s per seed.

---

## Config definitions

| Config       | xmult options                                  |
|--------------|------------------------------------------------|
| sect-only    | rss css xss nofuse noratchet nodrift           |
| sect+fuse    | rss css xss noratchet nodrift *(fuse=1 default)* |
| sect+ratchet | rss css xss ratchet 10 nodrift                 |
| sect+drift   | rss css xss drift 10 noratchet                 |
| all          | rss css xss ratchet 10 drift 10                |
| ratchet-only | norss nocss noxss ratchet 10 nofuse nodrift    |
| level 0–10   | xmult = level N giveupscore B hits 5 replic 100 |
| default      | xmult = giveupscore B hits 5 replic 100        |

All survey runs used `giveupscore B hits 5 replic 100`.  
TNT 1.6 quirk: `fuse` keyword inside `xmult =` triggers an interactive prompt; configs
wanting fuse simply omit `nofuse` (TNT default is fuse=1).

---

## Main results: median TTT (seconds, 5 seeds)

Configs ranked by median across all six datasets.

| Config       | Wortley | Eklund | Zanol | Zhu  | Giles | Dikow | **Median** |
|:-------------|--------:|-------:|------:|-----:|------:|------:|-----------:|
| sect+ratchet |    0.45 |   0.24 |  1.21 | 0.99 |  0.34 |  2.20 |   **0.56** |
| level3       |    0.34 |   0.25 |  1.68 | 1.01 |  0.45 |  1.98 |   **0.57** |
| level5       |    0.46 |   0.25 |  2.20 | 0.80 |  0.46 |  2.16 |   **0.57** |
| sect+drift   |    0.34 |   0.26 |  1.56 | 1.01 |  0.36 |  2.66 |   **0.63** |
| all          |    0.24 |   0.25 |  0.90 | 1.45 |  0.35 |  2.78 |   **0.63** |
| level0       |    0.58 |   0.26 |  3.66 | 0.79 |  0.37 |  2.54 |   **0.66** |
| level4       |    0.36 |   0.25 |  1.88 | 0.78 |  0.56 |  2.70 |   **0.67** |
| level10      |    0.35 |   0.24 |  2.32 | 0.90 |  0.56 |  3.62 |   **0.72** |
| level2       |    0.56 |   0.26 |  3.38 | 0.89 |  0.55 |  2.51 |   **0.73** |
| ratchet-only |    0.26 |   0.27 |  1.34 | 1.34 |  0.36 |  4.52 |   **0.74** |
| sect-only    |    0.56 |   0.26 |  2.65 | 1.10 |  0.45 |  3.30 |   **0.88** |
| sect+fuse    |    0.55 |   0.26 |  2.65 | 1.00 |  0.45 |  3.31 |   **0.90** |
| default      |    0.67 |   0.25 |  3.52 | 1.02 |  0.47 |  2.97 |   **0.95** |
| level1       |    0.79 |   0.25 |  3.96 | 1.00 |  0.46 |  3.30 |   **0.99** |

Cells are median of 5 seeds. 3 censored runs (level0/4/5 × Dikow2009, one seed each)
excluded from medians.

---

## Fastest config per dataset

| Dataset     | Tips | Fastest config | Median TTT |
|-------------|-----:|:---------------|------------|
| Wortley2006 |   37 | all            | 0.245 s    |
| Eklund2004  |   54 | sect+ratchet   | 0.238 s    |
| Zanol2014   |   74 | all            | 0.905 s    |
| Zhu2013     |   75 | level4         | 0.779 s    |
| Giles2015   |   78 | sect+ratchet   | 0.341 s    |
| Dikow2009   |   88 | level3         | 1.975 s    |

---

## Ratchet / drift verdict

Starting from the same sectorial baseline (rss+css+xss):

| Config       | Median TTT | vs sect-only |
|:-------------|------------|:-------------|
| sect-only    |  0.883 s   | baseline     |
| sect+fuse    |  0.903 s   | −2% (noise)  |
| sect+ratchet |  0.559 s   | **−37%**     |
| sect+drift   |  0.627 s   | **−29%**     |
| all          |  0.628 s   | **−29%**     |
| ratchet-only |  0.735 s   | −17%         |

**Adding ratchet or drift both accelerate convergence substantially.**  
Ratchet alone (+37%) outperforms drift alone (+29%) on average.  
Combining them (`all`) does not compound further — it matches drift alone.  
Fuse contributes nothing beyond what sectorial alone achieves.  
Ratchet-without-sectors works but is slower than either sectorial variant.

---

## Level series

The `level N` controls TNT's new-technology search effort. The response is non-monotonic:

| Level   | Median TTT | Notes                             |
|:--------|------------|:----------------------------------|
| level0  |  0.664 s   | Low effort; slow on Zanol         |
| level1  |  0.994 s   | **Worst in level series**         |
| level2  |  0.725 s   |                                   |
| level3  |  0.568 s   | **Best in level series**          |
| level4  |  0.672 s   |                                   |
| level5  |  0.574 s   | 2nd best                          |
| level10 |  0.715 s   |                                   |
| default |  0.949 s   | No level specified; 2nd worst     |

Levels 3 and 5 are the sweet spot. `default` (TNT's xmult without a level flag) performs
poorly — equivalent to level1 behaviour in this size range.

---

## Summary for emulation

1. **Best single config overall:** `sect+ratchet` (rss+css+xss, ratchet 10, no drift) —
   37% faster than plain sectorial, 40% faster than TNT default.

2. **Best for larger taxa (≥85t):** `level3` wins on Dikow2009 (88t); further testing on
   larger matrices needed to see if this scales.

3. **Fuse is inert** at these matrix sizes — neither helps nor hurts.

4. **Ratchet > drift**, but both help; combining both does not compound.

5. **Do NOT emulate TNT default** (`xmult;` with no extra options) — it ranks 13th of 14.

6. **Priority lever for TreeSearch emulation:** our sectorial search already implements
   rss+css+xss; adding a ratchet perturbation loop (10 iterations) is the single change
   most likely to close the remaining score gap observed in bench_sectorial_shared.R.
