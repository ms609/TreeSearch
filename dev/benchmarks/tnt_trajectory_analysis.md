# T-251: TNT vs TreeSearch Trajectory Analysis

Date: 2026-03-26

## Executive Summary

TreeSearch's score gap with TNT (3–21 steps on gap datasets) arises from two
compounding factors:

1. **Per-evaluation overhead**: TNT evaluates 1.5–3.6× more rearrangements
   per second than TreeSearch, despite TreeSearch having wider SIMD (SSE2
   128-bit vs TNT's 32-bit scalar on Windows). The overhead is in data
   structure manipulation, not the Fitch kernel.

2. **Phase allocation**: TreeSearch spends 16–23% of wall time on drift,
   which has extremely poor return (405–1498 ms per step gained). TNT's
   `xmult` is dominated by sectorial search, which is far more cost-effective.

## Methodology

Three datasets with the largest persistent score gaps (from T-249) were
compared at 30-second budgets, 3 seeds each, EW scoring, inapplicable
tokens treated as missing:

| Dataset | Tips | Chars | Gap (TS − TNT) |
|---------|:----:|:-----:|:---:|
| Geisler2001 | 68 | 186 | 5–9 |
| Zhu2013 | 75 | 253 | 4–6 |
| Wortley2006 | 37 | 105 | 3–4 |

TNT: console-mode Windows 32-bit (v1.6, 2026-02-20), `xmult=hits 10
replic 100`. TreeSearch: cpp-search HEAD, `ts_driven_search()` with
default strategy parameters, `verbosity=2`.

**Caveat:** TNT on Windows is 32-bit; Hamilton benchmarks will use the
64-bit Linux build which may have different throughput characteristics.
The per-evaluation throughput ratios below may not hold on Linux.

## Per-Evaluation Throughput

TNT's total rearrangements are reported directly. TreeSearch's
per-evaluation rate was measured via `ts_tbr_search()` on a single Wagner
→ TBR convergence.

| Dataset | TNT M evals/s | TS M evals/s | TNT/TS ratio |
|---------|:---:|:---:|:---:|
| Geisler2001 (68t) | 16.5 | 10.9 | 1.5× |
| Zhu2013 (75t) | 27.9 | 13.9 | 2.0× |
| Wortley2006 (37t) | 12.2 | 3.4 | 3.6× |

The gap is larger at smaller tree sizes, where the Fitch kernel is a
smaller fraction of per-evaluation cost and overhead dominates.

T-250 showed TreeSearch's Fitch kernel processes 128 bits per SIMD
iteration vs TNT's 32 bits — a ~4× raw throughput advantage. Yet TNT
evaluates more total rearrangements per second. This means TreeSearch's
**per-evaluation overhead** (undo stack management, data structure
traversal, incremental scoring setup) exceeds TNT's by 6–14×, completely
negating the SIMD advantage.

## Total Rearrangements (30s budget)

| Dataset | TNT total evals | TS est. total evals | TNT/TS ratio |
|---------|:---:|:---:|:---:|
| Geisler2001 | 499M | ~210M (est.) | ~2.4× |
| Zhu2013 | 796M | ~280M (est.) | ~2.8× |
| Wortley2006 | 104M | ~54M (est.) | ~1.9× |

TS estimates based on TBR throughput × phase time allocation. TNT
examines roughly twice as many candidates in the same wall time.

## Phase Cost Efficiency

TreeSearch phase efficiency = ms of wall time per step of score
improvement. Lower is better. Averaged over 3 seeds per dataset.

### Geisler2001 (68 taxa)

| Phase | Time (ms) | Steps gained | ms/step | % of time |
|-------|:---------:|:---:|:---:|:---:|
| TBR | 1773 | 2397 | 0.8 | 3% |
| CSS | 1408 | 154 | 9.1 | 3% |
| RSS | 863 | 49 | 18 | 2% |
| XSS | 1502 | 97 | 20 | 3% |
| Ratchet | 34616 | 1070 | 34 | 63% |
| **Drift** | **11843** | **8** | **1498** | **22%** |

### Zhu2013 (75 taxa)

| Phase | Time (ms) | Steps gained | ms/step | % of time |
|-------|:---------:|:---:|:---:|:---:|
| TBR | 1574 | 3321 | 0.5 | 3% |
| XSS | 2028 | 367 | 5.7 | 4% |
| CSS | 1372 | 107 | 14 | 3% |
| RSS | 833 | 46 | 18 | 2% |
| Ratchet | 33710 | 765 | 44 | 62% |
| **Drift** | **12695** | **10** | **1270** | **23%** |

### Wortley2006 (37 taxa)

| Phase | Time (ms) | Steps gained | ms/step | % of time |
|-------|:---------:|:---:|:---:|:---:|
| TBR | 1100 | 2655 | 0.4 | 2% |
| XSS | 1652 | 376 | 4.5 | 3% |
| CSS | 1332 | 226 | 6.1 | 3% |
| RSS | 883 | 83 | 11 | 2% |
| Ratchet | 35945 | 2058 | 18 | 72% |
| **Drift** | **7989** | **22** | **405** | **16%** |

**Pattern:** Drift is 30–170× less efficient than the next-worst phase
(ratchet) across all three datasets.

## TNT's Search Structure

TNT's `xmult` trajectory reveals a fundamentally different phase
composition from TreeSearch's pipeline:

**Geisler2001 (30s, seed 1):** TNT reports 30 sub-replicate results
across 7 replicates. Algorithm breakdown:
- SECT (sectorial search): ~20 entries
- TBR: ~8 entries
- FUSE: ~2 entries

TNT hits score 1293 within replicate 0 (3 seconds, 56M rearrangements)
via TBR following sectorial search. Subsequent replicates hover around
1293–1303, with sectorial search and fusing maintaining the best score.

TreeSearch hits 1298 as its best single-replicate score (replicate 10,
after ~14s of cumulative search time). No replicate reaches 1293.

**Key structural differences:**

1. **TNT does extensive sectorial search within each replicate.** Each TNT
   replicate includes multiple rounds of sectorial search + TBR before
   moving to the next Wagner start. TreeSearch does one pass of
   XSS+RSS+CSS per outer cycle.

2. **TNT's replicates are longer and more productive.** TNT completes ~7
   replicates in 30s on Geisler2001 (~4.3s each), with each replicate
   including intensive sectorial + TBR + fuse. TreeSearch completes 14–19
   replicates (~1.5–2s each), but each is shallower.

3. **TNT fuses frequently within the search.** The FUSE entries in TNT's
   trajectory show tree fusing as an integrated part of the search cycle,
   not a separate post-search step.

## Per-Replicate Score Quality

Median per-replicate score (the typical quality of a single search from
a random Wagner start):

| Dataset | TNT median rep | TS median rep | TNT advantage |
|---------|:---:|:---:|:---:|
| Geisler2001 | ~1297 | 1313 | 16 steps |
| Zhu2013 | ~626 | 636 | 10 steps |
| Wortley2006 | ~487 | 488 | 1 step |

TNT achieves better per-replicate scores, which means its intra-replicate
search (sectorial + TBR) is more thorough.

## TreeSearch Per-Replicate Trajectory

**Geisler2001 (seed 1):** 15 replicates
- Rep 1: 1349 (Wagner 1678 → TBR → Ratchet → Drift)
- Rep 2: 1308 (improvement)
- Rep 5: 1304
- Rep 10: 1298 (best found)
- Rep 15: 1327 (no improvement in last 5 reps)

Score improves from 1349 → 1298 over 15 replicates (51 steps). TNT
improves from ~1298 → 1293 within a single replicate.

## Recommendations

### High priority: Eliminate or drastically reduce drift

Drift consumes 16–23% of search time but contributes <1% of score
improvement. At 405–1498 ms per step gained, it is 30–170× less
efficient than the next-worst phase.

**Proposed change:** Set `driftCycles = 0` in the default preset.
Reallocate the saved time to additional ratchet cycles or sectorial
search rounds. The `thorough` preset (with many more base cycles) could
retain 1–2 drift cycles as a diversity mechanism.

Expected impact: ~20% wall-time savings with negligible score loss.
Equivalent to adding ~4 more replicates per 30s budget.

### Medium priority: Increase sectorial search intensity

TNT's dominance of sectorial search (SECT appears in ~67% of trajectory
entries) suggests TreeSearch's single-pass XSS+RSS+CSS is insufficient.
Currently sectorial search takes only 6–10% of wall time but has
respectable efficiency (5–20 ms/step).

**Proposed change:** Increase sectorial search rounds. Options:
- Double `xssRounds` and `rssRounds` within each outer cycle
- Add a second sectorial search pass after ratchet (currently
  sectorial → ratchet → drift → TBR; change to
  sectorial → ratchet → sectorial → TBR)
- Increase `sectorMaxSize` to capture more of the tree in each sector

### Medium priority: Reduce per-evaluation overhead

The 1.5–3.6× per-evaluation throughput gap means every search phase is
penalized. Likely targets:
- Undo stack management in TBR (PreallocUndo grow/shrink)
- Incremental scoring setup cost (even when not finding improvements)
- Collapsed-flag recomputation (O(n) per move, even when 0% collapsed)

This is a deeper engineering effort (T-245/T-246 overlap) but has the
broadest impact since it accelerates every phase.

### Low priority: Ratchet tuning

Ratchet is the most time-consuming phase (62–72%) and mid-tier in
efficiency. The current 12 cycles at 25% perturbation may be too many;
diminishing returns likely set in after 6–8 cycles. The adaptive level
mechanism already scales this down when hit rates are high, but the
base count could be reduced for the default preset.

## Data Files

- `bench_trajectory.R` — comparison script
- `trajectory_results.rds` — raw results (3 datasets × 3 seeds)
- `tnt_trajectory_analysis.md` — this document
