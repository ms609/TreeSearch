# TNT 1.6 Scaling Survey — Time-to-Target (TTT), MorphoBank datasets

**Companion to:** `tnt_settings_survey.md` (6 gap datasets, 37–88 taxa)  
**Datasets:** 4 MorphoBank projects, 103–205 taxa (above the n=90 sector-size inflection)  
**Configs:** same 14 as the gap-dataset survey  
**Seeds:** 3 per (config, dataset); **timeout:** 300 s/run  
**Phase-1:** 5 seeds × 600 s to establish B

**Machine:** DW-CZC429715G · 12th Gen Intel Core i7-12700 · 15.7 GB RAM  
**TNT:** C:/Programs/Phylogeny/tnt/tnt.exe (v1.6, 32-bit)  
**Date:** 2026-06-17  
**Script:** `dev/benchmarks/bench_tnt_settings.R` → `tnt_scaling_full()`  
**Data:** `dev/benchmarks/tnt_scaling_survey.csv`

---

## Dataset reference

| Dataset    | Tips | Chars | B (Phase-1) | Note                           |
|------------|-----:|------:|------------:|:-------------------------------|
| project691 |  103 |   529 |        2169 | just above n=90 inflection     |
| project4230|  125 |   307 |        1149 |                                |
| project4103|  144 |   169 |         671 |                                |
| project3763|  205 |   109 |        1292 | Phase-2 found 1290 → B tight  |

B from Phase-1 (5 seeds × 600 s). For project3763 several Phase-2 configs
found 1290–1291, so the Phase-1 B is slightly conservative; all comparisons
treat ≤ 1292 as "reached".

---

## Main results: median TTT (seconds, 3 seeds)

`NA` = all 3 seeds censored (never reached B within 300 s).

| Config       | 103t (p691) | 125t (p4230) | 144t (p4103) | 205t (p3763) |
|:-------------|------------:|-------------:|-------------:|-------------:|
| sect-only    |        0.48 |         1.77 |         0.27 |       **NA** |
| sect+fuse    |        0.50 |         1.00 |         0.25 |       300.4  |
| sect+ratchet |        0.81 |         1.46 |         0.29 |      **86.6**|
| sect+drift   |        1.00 |         1.11 |         0.35 |      **57.5**|
| all          |        0.98 |         1.84 |         0.46 |      112.8   |
| ratchet-only |        1.22 |         4.97 |         0.47 |       **NA** |
| level0       |        1.11 |         2.15 |         0.49 |       **NA** |
| level1       |        0.80 |         1.57 |         0.33 |       **NA** |
| level2       |        0.70 |         1.45 |         0.37 |       **NA** |
| level3       |        0.81 |         1.58 |         0.35 |       **NA** |
| level4       |        0.72 |         1.69 |         0.27 |       **NA** |
| level5       |        0.68 |         1.59 |         0.46 |       **NA** |
| level10      |        1.00 |         1.67 |         0.46 |       **NA** |
| default      |        0.70 |         1.56 |         0.37 |      300.5   |

Censored totals: project691 0/42, project4230 0/42, project4103 0/42,
project3763 **29/42**.

---

## The 205-taxon failure

At 205 taxa, TNT's sector size is pinned at 45 (the n=90–450 plateau). With
only ~4–5 sectors tiling the tree, the within-sector RAS restarts cover a
small fraction of tree space. Without global perturbation, the search stalls.

| Config at 205t | Median score | Seeds reached B |
|:---------------|-------------:|----------------:|
| sect+ratchet   |         1290 |             3/3 |
| sect+drift     |         1291 |             3/3 |
| sect+fuse      |         1291 |             3/3 |
| all            |         1291 |             3/3 |
| default        |         1296 |             1/3 |
| level0–10      |     1294–NA  |             0/3 |
| ratchet-only   |           NA |             0/3 |
| sect-only      |           NA |             0/3 |

sect+ratchet actually finds 1290 (better than Phase-1 B=1292), confirming the
Phase-1 search was not exhaustive enough for this matrix.

---

## Scaling verdict

### At ≤ 144 taxa

All configs reach B quickly (median < 5 s). The ordering broadly matches the
gap-dataset survey: sect+fuse and sect-only are fastest, ratchet-only is
slowest. Ratchet and drift add overhead with little benefit at this scale.

### At 205 taxa — the picture inverts

Pure sectorial (sect-only, all level variants, ratchet-only) **fail
completely**. Perturbation is no longer optional — it is load-bearing:

- **sect+drift wins** (57 s) — drift perturbation rescues stalled sectorial
- **sect+ratchet** (87 s) — ratchet also rescues, slightly slower
- **ratchet-only fails** — perturbation alone without sectorial is insufficient
- **level configs fail** — TNT's own level-based heuristic cannot find B

This reversal happens somewhere between 144 t and 205 t. The 29/42 censor
rate on project3763 (a 109-char matrix, not unusually complex) suggests the
transition is taxon-driven rather than character-driven.

---

## Combined verdict (37–205 taxa)

| Config       | ≤ 88t rank | 103–144t rank | 205t outcome |
|:-------------|:----------:|:-------------:|:-------------|
| sect+ratchet |  **1st**   |    middle      | ✓ reaches B  |
| sect+drift   |    4th     |    middle      | ✓ fastest    |
| all          |    5th     |    slowest     | ✓ reaches B  |
| sect+fuse    |   11th     |    fastest     | barely (300s)|
| sect-only    |   11th     |    fastest     | ✗ fails      |
| level3       |  **2nd**   |    middle      | ✗ fails      |
| ratchet-only |   10th     |    slowest     | ✗ fails      |
| default      |   13th     |    middle      | barely (300s)|

**sect+ratchet is the only config that is fast at small scale AND reliable at
large scale.** This makes it the unambiguous emulation target.

---

## Implications for TreeSearch emulation

1. **Ratchet is mandatory, not optional** — at 205t, no-perturbation configs
   all fail. The ratchet loop must be part of the core search, not an
   add-on.

2. **Sectorial without perturbation stalls beyond ~150t** — even `sect-only`
   with TNT's native 3-RAS-per-sector fails. The sector resampling alone is
   insufficient once the tree is large enough that sectors cover little of
   the tree each pass.

3. **Level configs scale poorly** — TNT's own `level N` parameter does not
   help beyond 144t in our tests. This is consistent with level being a
   within-pass intensity setting, not a global perturbation mechanism.

4. **Priority fix for `search_sector` ([src/ts_sector.cpp:501](../../src/ts_sector.cpp)):**
   add k RAS+TBR restarts per sector (as per the shared-start finding), AND
   add a ratchet outer loop. Both are needed for large-matrix performance.
