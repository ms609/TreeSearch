# Search Algorithm Design Notes

Load this when: researching NNI warmup, biased Wagner, outer cycles,
large-tree behaviour, or reviewing the search optimization history.

See also: `search_strategy.md` (pipeline structure, strategy presets),
`benchmarking.md` (corpus, methodology, benchmark tables).

---

## NNI in the driven pipeline

`nni_search()` in `ts_search.cpp` is implemented. At ≤88 tips, NNI is
redundant — TBR subsumes it. At 180 tips, NNI becomes essential: TBR
evaluates O(n²) candidates per pass (millions of evaluations, many minutes
to converge from Wagner); NNI evaluates O(n) candidates (~1000× cheaper).

**All presets set `nniFirst = TRUE`** (NNI warmup before TBR). Each Wagner
start is NNI-optimized before selection. SPR is counterproductive when NNI
is active — NNI→TBR outperforms NNI→SPR→TBR empirically.

**Empirical comparison at 180 tips** (mbank_X30754, 3 seeds, EW):

| Strategy | Median score | Median time |
|----------|:-----------:|:-----------:|
| TBR alone | 1427 | 13.6s |
| SPR→TBR | 1360 | 13.1s |
| **NNI→TBR** | **1326** | **6.8s** |
| NNI→SPR→TBR | 1369 | 8.8s |

NNI→TBR wins on both score AND time (~2× faster, ~100 steps better).

**Time-adjusted expected best (5 seeds, EW):**

| Budget | 88t: TBR | 88t: NNI→SPR→TBR | 180t: TBR | 180t: NNI→SPR→TBR |
|--------|:--------:|:-----------------:|:---------:|:-----------------:|
| 20s | 1617 | 1619 (+2) | 1388 | 1278 (−110) |
| 60s | 1617 | 1619 (+2) | 1348 | 1253 (−95) |
| 120s | 1617 | 1619 (+2) | 1337 | 1247 (−90) |

At ≤88 tips: NNI has a consistent but negligible 2-step penalty. At 180
tips: NNI saves 90–110 steps. No reactive per-run switching needed — always-on
NNI warmup is optimal.

---

## Stochastic NNI-perturbation (T-186)

`ts_nni_perturb.h/.cpp` implements topology-space escape inspired by
IQ-TREE's `doRandomNNIs()`. Complementary to the weight-perturbation
ratchet: ratchet reshapes the objective function; NNI-perturbation directly
displaces the tree topology.

**Algorithm:** Collect all internal NNI edges. For each edge (with probability
`perturb_fraction`, default 0.5), apply a random NNI swap — skip edges
adjacent to already-swapped edges. Track touched nodes in a hash set.
After all compatible swaps, rebuild postorder and full rescore, then TBR
to a new local optimum. Repeat for `n_cycles`.

**Pipeline placement:** Between ratchet and drift. **Disabled by default
(`nniPerturbCycles = 0`)** and in all presets since T-274 (2026-03-27).

**R API:** `SearchControl(nniPerturbCycles, nniPerturbFraction)`.

**T-274 benchmark (2026-03-27):** 20 seeds, Zhu2013/Giles2015/Dikow2009
(75–88t). NNI-perturb adds 59–69% per-replicate overhead with ≤0.1-step
expected-best benefit at all budgets — within bootstrap noise. Set
`nniPerturbCycles = 0` in thorough preset. Available via `SearchControl()`
for manual use.

---

## Biased Wagner addition (T-188, 2026-03-23)

`biased_wagner_tree()` (`ts_wagner.h/.cpp`) samples the taxon-addition order
from a softmax distribution weighted by informativeness score.

Two criteria:
- **GOLOBOFF** (bias=1): `score[t]` = number of non-ambiguous characters for
  taxon t. Ref: Goloboff 2014 (*Extended implied weighting*) §3.3.
- **ENTROPY** (bias=2): `score[t]` = Σ_c (n_states_c − |state set for t|).

**R API:** `SearchControl(wagnerBias = 0L, wagnerBiasTemp = 0.3)`.
Applied only to the first of `wagnerStarts` starts; remaining starts use
random order for basin diversity.

**Benchmark results** (2026-03-23, 14 standard + crico-174):
- Wagner→TBR gap reduction: ~80% at 174t (random: 1356 steps, Goloboff: 244)
- Score improvement after TBR convergence: ~22 steps at 174t; 1–2 steps at ≤88t
- Anomalous slight regression at 75–100t; T=0.3 stochastic is safer than T=0

---

## Outer search cycle loop (T-189, 2026-03-23)

`outer_cycles` in `SearchParams` / `outerCycles` in `SearchControl()`.
Wraps steps 3–6 of `run_single_replicate()` in a configurable outer loop:
`[XSS+RSS+CSS → Ratchet → NNI-perturb → Drift → TBR] × N`.
Ratchet/NNI-perturb/drift cycles are divided evenly among N outer cycles.

`outerCycles = 1` (default) is bit-for-bit identical to the previous
linear pipeline. `thorough` preset defaults to `outerCycles = 2`.

Matches TNT's `xmult` interleaving (Goloboff 1999 §2.3): after each
ratchet/drift escape, a fresh XSS pass exploits the new topology.

---

## Large-tree scaling issues (discovered 2026-03-23)

The 180-taxon `mbank_X30754` dataset (425 chars, 374 informative patterns,
40% missing, 20% inapplicable) exposed:

1. **C++ TBR convergence at 180 tips takes ~13s** (Wagner ~2560 → local
   optimum ~1420). NNI warmup (~1.5s) followed by TBR reduces this to
   ~7s while finding better scores. T-178 filed.
2. **Strategy presets assume replicate time O(seconds).** At 180 tips,
   a single replicate takes ~60-100s. Cycle counts need recalibration.

**180-taxon baseline (C++ driven search, EW, single replicate):**
- Wagner (best of 3): ~2560 steps, 16ms
- NNI convergence: ~1600 steps, 1.5s
- TBR convergence: ~1330 steps, 7s (from NNI-optimal start)
- XSS: additional ~60 steps improvement, 5s
- Total single replicate: ~25s (before ratchet/drift)

---

## Search optimization roadmap

Items completed as of 2026-03-29. Numbered by original priority.

1. ~~Consensus-guided sector targeting~~ — **Done**: RSS weighted by pool split conflict scores
2. ~~Diverse pool maintenance~~ — **Done**: evict most-similar entry on ties
3. ~~Cross-replicate constraint tightening~~ — **Done**: opt-in via `consensusConstrain = TRUE`
4. ~~Collapsed-tree clip skipping~~ — **Done**: zero-length edges skipped in TBR, SPR, drift. Skip rate 0% on standard morphological datasets (benefit expected on sparse/synthetic data).
5. ~~Collapsed-region regraft merging + pool dedup~~ — **Done**: boundary-only regraft evaluation; collapsed-topology pool dedup.
6. ~~Strategy preset tuning~~ — **Done**: `default` uses `wagnerStarts=3`, `sprFirst=TRUE`, `adaptiveLevel=TRUE`; `thorough` uses `sprFirst=TRUE`.
7. ~~Ratchet perturbation tuning~~ — **Done**: perturbation probability 4%→25%, perturbed TBR moves 20→5, ratchet cycles 5→10 (default), 10→20 (thorough). Drift cycles 2→4, AFD 5, RFD 0.15. Validated on 14 datasets.
8. ~~Biased Wagner addition~~ — **Done** (T-188): see above.
9. ~~Outer search cycle loop~~ — **Done** (T-189): see above.
10. ~~Drift MPT diversity experiment~~ — **Done** (T-254): drift provides zero score benefit, zero MPT enumeration benefit. Delays consensus stability. `driftCycles=0` in all presets (T-255).
11. ~~NNI-perturb cycle count at thorough-preset scale~~ — **Done** (T-274): see above.
12. ~~Size-weighted TBR clip ordering~~ — **Closed** (2026-03-29): Hypothesis FALSIFIED. Tip clips (~51% of all clips) account for only 22–38% of accepted moves (enrichment 0.43–0.76×). Medium-small clips (size 2..√n) are most productive. All three variants (INV_WEIGHT, TIPS_FIRST, BUCKET) favour tips — wrong direction. Diagnostic code preserved in `feature/weighted-clip-order` branch.
13. ~~XSS↔TBR cycling under IW~~ — **Closed** (2026-03-29): IW3 XSS improvement rate ~30% vs EW ~25%; below 2× threshold. Key finding: XSS cycling benefit scales with tree size, not scoring mode. At 180t: XSS adds 12–19% overhead, TAEB Δ = −6.8 to −9.8 EW steps at 30–120s.
14. ~~Targeted post-clip sector search~~ — **Closed** (2026-03-29): Hit rate ~35% but net HARMFUL — local sector refinement after each TBR move changes global trajectory, steering into worse basins. Validates existing design: XSS should run as a separate phase AFTER TBR convergence.
