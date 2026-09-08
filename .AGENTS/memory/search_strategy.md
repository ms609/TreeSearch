### Driven search pipeline per replicate

1. Random Wagner tree → NNI warmup → TBR to local optimum
2. XSS sectorial search (if tree large enough)
3. RSS random sectorial search
4. CSS constrained sectorial search
5. Ratchet perturbation to escape local optima
5a. Post-ratchet XSS+RSS+CSS (if `postRatchetSectorial = TRUE`)
6. NNI-perturbation (topology-space escape, if `nniPerturbCycles > 0`)
7. Drift search (accept suboptimal moves)
8. PCSA perturbation (if `annealCycles > 0`)
9. Final TBR polish
10. Add to pool
11. Fuse against pool (every `fuse_interval` replicates)

Steps 2–9 are wrapped in the `outerCycles` loop (default 1).

Post-search: TBR plateau enumeration from all pool seeds to find MPTs.

### Strategy presets (auto-selected by `NTip` and signal density)

| Preset | Condition | Key settings |
|--------|-----------|-------------|
| sprint | ≤30 tips | 3 ratchet (4%), 0 drift, XSS only, NNI-first |
| default | 31–64 tips; or ≥65 tips with <100 char patterns | 12 ratchet (25%, 5 moves), 0 drift, XSS+RSS, Wagner×3, NNI-first, adaptive level |
| thorough | 65–119 tips with ≥100 char patterns | 20 ratchet (25%, 5 moves, adaptive), 0 NNI-perturb (T-274), 0 drift, XSS+RSS+CSS, Wagner×3, NNI-first, outerCycles=2 |
| large | ≥120 tips with ≥100 char patterns | 12 ratchet (25%, 5 moves, adaptive), 0 NNI-perturb, 0 drift, 1 SA cycle (T=20→0, 5 phases), XSS(3)+RSS(2)+CSS(1), Wagner×1 biased (Goloboff 2014), NNI-first, outerCycles=1, tbrMaxHits=1, sectorMaxSize=100, pruneReinsert=5 cycles NNI-polish (T-289f Stage 5: NNI polish fixes 0-rep failure at 206t; improves 131–180t) |

**T-264 (2026-03-26):** `consensusStableReps` removed from all presets
(disabled, 0). The previous setting of 3 caused catastrophic early
termination — the search stopped after 3 replicates with unchanged
consensus, using only 7–20% of the time budget on most datasets.

**Large preset design rationale (T-179, 2026-03-24):** At 180 tips, each TBR
convergence takes ~5–7s, making phases like NNI-perturbation (~5.5s/cycle) and
drift (~4s/cycle) extremely expensive. Systematic benchmarking on mbank_X30754
(180t, 418p) showed that reducing cycle counts (12 ratchet, 4 drift, no NNI-perturb)
with outerCycles=1 and a single biased Wagner start outperforms the thorough
preset by 4–7 steps (median) at 30–60s budgets and ties at 120s, while
consistently completing more replicates.

**T-289 Stage 4 (2026-03-28, EPYC 7702, 10 seeds, 5 datasets 131–206t):**
PR (c=5, d=5%, MISSING) vs baseline. 60s: mean Δ=+0.5 steps (neutral);
project3701 146t regresses −12 steps; syab07205 206t: 0 replicates complete
(per-rep cost ~60s, budget exceeded). 120s: mean Δ=−9.1 steps but driven
by project3701 (−37 steps); others ≤6 steps. Replicate ratio 0.82 at 60s,
0.68 at 120s. Decision: disable PR (TBR polish) — 0-rep failure at 206t/60s
is a showstopper.

**T-289f Stage 5 (2026-03-29, EPYC 7702, 10 seeds, 5 datasets 131–206t):**
PR (c=5, NNI full-tree polish) vs pr_tbr (TBR polish, Stage 4 reference) vs
baseline. pr_tbr at 206t/60s: still 0 reps (confirmed). pr_nni fixes the
0-rep failure (2 reps at 206t/60s). Score deltas vs baseline: project4133
(131t) ≈0; project3701 (146t) **−178 steps** at 60s, −128 at 120s; project804
(173t) −9/−2; mbank_X30754 (180t) −4/−7; syab07205 (206t) +17.5 at 60s
(neutral at 120s). **Decision: enable pruneReinsertCycles=5, pruneReinsertNni=TRUE
in large preset.** Note G-006: NNI polish ignores ConstraintData — irrelevant
since large preset does not use topological constraints.

**Post-T-206 Hamilton HPC baselines (2026-03-26, EPYC 7702, 5 seeds):**
30s median=1202 (range 1189–1214), 60s median=1190 (1190–1202), 120s
median=1185 (1171–1189). Per-replicate median 17.3s (cf. ~60s pre-T-206).
The 65–74 step improvement over pre-T-206 Intel baselines is primarily
from the outer cycle reset cap (maxOuterResets=0), not hardware.
Phase distribution: TBR 43.6%, Ratchet 32.2%, SA 7.4% (14% hit rate,
0.8 steps/s — least productive phase). T-248 benchmarked annealCycles
0/1/3: AC=1 (400ms/rep, 40% hit rate) is most cost-effective; AC=3
(1370ms/rep, 21% hit rate) showed no significant score gain (p>0.5,
n=5 seeds). Large preset reduced to annealCycles=1.

All presets set `nniFirst = TRUE` (NNI warmup before TBR) and
`sprFirst = FALSE` (SPR is counterproductive when NNI is active —
empirically NNI→TBR outperforms NNI→SPR→TBR). With `nniFirst`, each
Wagner start is NNI-optimized before selection (best of 3 NNI-local optima
rather than 3 raw Wagner scores). `default` also enables `adaptiveLevel =
TRUE` (scale ratchet/drift by hit rate); `thorough` omits it because high
base cycle counts already cover hard landscapes.

**Ratchet perturbation tuning (2026-03-22)**: Systematic profiling across
all 14 benchmark datasets showed the previous 4% perturbation probability
was far too gentle. With 253 characters (Zhu2013), 4% zeroes only ~10
characters — insufficient to reshape the landscape. Increasing to 25%
with fewer perturbed TBR moves (5 instead of auto=20) improves median
scores by 3–7 steps on hard datasets while completing fewer but more
productive replicates. 9/14 datasets improved, 4 unchanged, 1 marginal at
10s budget (resolves at 20s). The key insight: the perturbed-phase TBR
should be short (the landscape is warped, so extensive search on it is
wasteful), but the perturbation itself should be aggressive enough to
meaningfully displace the tree from its current basin of attraction.

Signal-density gate: datasets with few character patterns (<100) have flat
parsimony landscapes where intensive search adds no benefit.

### Adaptive sectorial search

XSS and CSS use **adaptive early-exit**: after each round of sector searches
+ global TBR polish, if the overall best score did not improve, remaining
rounds are skipped. This avoids wasting ~7% of replicate time on datasets
where sectorial search is unproductive (e.g. Dikow2009). On productive
datasets (e.g. Zhu2013), the early exit never fires.

### Conflict-guided RSS

RSS uses **conflict-guided sector selection**: before each replicate's RSS
phase, `driven_search()` computes a `SplitFrequencyTable` from the pool's
best-score trees. Within `rss_search()`, each internal node's "conflict
score" is `1 − (fraction of pool trees containing that split)`.
Max-descendant conflict is propagated upward, and eligible sector roots
are sampled via `std::discrete_distribution` with weight `1 + 3 × conflict`.
Falls back to uniform selection when the pool has <2 best-score trees or
when conflict variation is negligible.

### Consensus-stability stopping

After each replicate, if `consensus_stable_reps > 0` (disabled in all
presets since T-264; available via `SearchControl(consensusStableReps=N)`),
the pool's strict consensus hash is compared to the previous replicate's.
If unchanged for `consensus_stable_reps` consecutive replicates, the
search terminates early. `compute_consensus_hash()` uses
XOR of per-split FNV-1a hashes for O(pool × splits) cost.

### Adaptive search level

When `adaptive_level = true`, ratchet and drift cycle counts are scaled
each replicate based on the cumulative hit rate:
- hit_rate > 0.7 → 0.5× (easy landscape)
- hit_rate > 0.4 → 0.75×
- hit_rate < 0.15 → 1.5× (hard landscape)
- else → 1.0×

### TBR zero-length clip skipping + regraft merging (collapsed flags)

`compute_collapsed_flags()` (`ts_collapsed.h/.cpp`) identifies edges where
clipping provably cannot improve score. Checks 5 conditions: (1) zero
standard-block cost at parent, (2) zero NA-block cost at parent, (3) prelim
preservation (`prelim[sibling] == prelim[parent]`), (4) down2 preservation
(NA), (5) subtree_actives preservation (NA). Works for EW, IW, Profile,
and NA-aware scoring. Integrated into TBR, SPR, and drift search.
Disabled during MPT enumeration (equal-score topologies may exist).
Recomputed after every accepted move.

**Regraft merging** (Goloboff 1996): within a collapsed region (connected
set of nodes linked by zero-length edges), all regraft positions yield the
same full score. Only boundary edges (entering the region) are evaluated;
interior collapsed edges are skipped via `if (collapsed[below]) continue`.
TBR, SPR, and drift all use this. The `CollapsedRegions` struct exists in
the header but callers use `compute_collapsed_flags()` directly (the
`region_id` field is unused — only the boolean flag array matters).

**Collapsed-topology pool dedup**: `compute_collapsed_splits()` in
`ts_splits.cpp` produces the split set excluding collapsed edges. Two
binary trees differing only in zero-length resolutions produce the same
collapsed split set → treated as duplicates by `TreePool::add_collapsed()`.
Both serial (`driven_search`) and parallel (`ThreadSafePool`) paths use
collapsed dedup.

**Benchmark results** (2026-03-22, 4 standard datasets, 3 seeds each):
Skip rate = 0% on all datasets (Vinther2008 23t, Agnarsson2004 62t,
Zhu2013 75t, Dikow2009 88t). Near-optimal trees in these morphological
datasets have negligible zero-length edges. Overhead from flag computation
is negligible. Score equivalence confirmed (enabled vs disabled produce
identical best scores). Benefit expected on sparse/synthetic data.
