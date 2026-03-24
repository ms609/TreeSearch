# T-199 Findings: PT Diagnostic Profiling & PCSA Evaluation

## Executive Summary

**Standard parallel tempering with Boltzmann acceptance is fundamentally
incompatible with discrete parsimony.** Pair 0–1 (cold↔warm) swap
acceptance is exactly 0% across all 8 test matrices (40–385 tips), all
temperature ladders, and all moves-per-round settings.

However, **post-convergence simulated annealing (PCSA) with best-tree
restart is highly effective.** Multi-cycle PCSA outperforms cold-only TBR
by 40–110 steps at 125–205 tips (EW) with dramatically lower variance
(SD 10 vs 62 at 205 tips). PCSA functions as a "smooth ratchet" —
same principle (perturb, reconverge, keep best), different perturbation
mechanism.

---

## Part 1: Standard PT is broken for parsimony

### Phase 1: Default ladder {0, 3, 9, 27} across 8 matrices

| Matrix | Tips | Chars | Best | Pair 0-1 | Swaps | Cold improve | Cold ms | Hot ms |
|--------|------|-------|------|----------|-------|-------------|---------|--------|
| project3419 | 40 | 368 | 2062 | 0% | 2/30 | 0 | 78 | 16 |
| project2604 | 43 | 307 | 1115 | 0% | 12/30 | 0 | 320 | 74 |
| project4531 | 71 | 256 | 885 | 0% | 4/30 | 0 | 240 | 37 |
| project3741 | 86 | 110 | 849 | 0% | 6/30 | 0 | 491 | 268 |
| project3253 | 125 | 394 | 5435 | 0% | 10/30 | 0 | 8013 | 1219 |
| project1221 | 150 | 252 | 874 | 0% | 6/30 | 0 | 4160 | 1131 |
| project3763 | 205 | 105 | 1550 | 0% | 7/30 | 0 | 8910 | 1186 |
| project2722 | 385 | 520 | 5584 | 0% | 5/30 | 0 | 220179 | 9601 |

Swaps succeed between hot chain pairs (1-2, 2-3) at 10–90%. The
bottleneck is exclusively at the cold-warm boundary.

### Phase 2–3: Tighter ladders and fewer moves

- Ladders {0, 0.5, 1, 2} and {0, 1, 3, 9}: pair 0-1 still 0%.
- Moves per round = 5, 2, 1: still 0%.
- Ultra-low T (0.01–0.1): warm chain has 72% acceptance (improving moves)
  but diverges by 400–900 steps because stochastic TBR samples random
  moves while cold chain does exhaustive best-improvement sweeps.

### Phase 4: Warm cold chain (T_cold > 0)

Tested T_cold = 0.1, 0.3, 0.5, 1.0:
- T_cold = 0.3: Metropolis prob ≈ 10^-256
- T_cold = 1.0: Metropolis prob ≈ 10^-5
- Chain 0 always runs deterministic TBR regardless of temperature;
  T_cold only affects the swap criterion.

### Phase 5: Cold-only vs 4-chain PT

- 4-chain PT (20 rounds, 205 tips): best = 1550, cold_improve = 0
- 1-chain cold only (20 rounds): best = 1550, 80% of wall time
- Hot chains are pure overhead.

### Root cause

In discrete parsimony with TBR, the cold chain (T=0) converges almost
instantly via deterministic TBR. Any T>0 chain diverges because each
stochastic TBR move changes many edges at once, incurring 10–100+ step
penalties. The score gap at the cold-warm boundary (89–1400+ steps) makes
Metropolis acceptance probability effectively zero. Parsimony lacks the
"gradually overlapping energy distributions" property that PT requires.

---

## Part 2: IW landscape investigation

Under implied weighting (IW, k=3), scores are fractional. The pair 0-1
score gap shrinks from ~760 (EW) to ~22 (IW units).

### IW with high T_cold

With T_cold ≥ 5, pair 0-1 swaps succeed at 5–60%:

| Ladder (86 tips, IW k=3) | Pair 0-1 rate | Best score |
|---------------------------|---------------|------------|
| T=0,3,9,27 | 0% | 59.14 |
| T=3,6,12,24 | 0% (prob=0.02) | 59.14 |
| T=5,10,20,40 | 5% | 58.70 |
| T=10,20,40,80 | 40% | 59.48 |
| T=20,40,80,160 | 60% | 59.85 |

At T=5, a successful swap displaced the cold chain from its 59.16 local
optimum; TBR descent found a different basin reaching 58.70 (a real
improvement). But `cold_improve=0` because the check is immediate
post-swap (score 82.36), not post-TBR-recovery.

### IW multi-seed comparison (86 tips, 10 seeds)

| Approach | Mean | SD | Time |
|----------|------|----|------|
| Cold-only | **58.76** | 0.68 | 1× |
| PT(T=0) | 59.00 | 0.68 | 1.2× |
| PT(T=5) | 58.88 | **0.41** | 3× |

Under IW, cold-only still wins on average. PT(T=5) has lower variance
(forced restarts provide diversification) but at 3× the cost. The IW
landscape is smooth enough that TBR navigates it well without perturbation.

---

## Part 3: Alternative exchange mechanisms

### Score-transfer PT

Replaces Metropolis swaps with direct injection: when any hot chain's
current score beats the cold chain, inject that tree immediately. Simple,
no temperature dependence.

**Result:** Marginal. At 205 tips, one notable win (IW seed 9042: 66.77
vs cold 68.71) but on average equivalent to cold-only. Hot chains at
T≥3 almost never produce trees better than the cold chain's converged
score.

### Post-convergence SA (PCSA) — from scratch

Wagner tree → SA (t_start → 0) → TBR polish. Compared to cold-only
TBR (Wagner → TBR convergence):

| 205 tips EW (5 seeds) | Mean | Time |
|-----------------------|------|------|
| Cold-only (TBR×20) | 1456 | 1× |
| PCSA (t=20→0, 5 phases) | **1399** | ~2× |

**PCSA beats cold-only by ~57 steps** at 205 tips EW. The SA perturbation
displaces the tree from its local optimum; TBR polish finds a different,
often better basin.

**Under IW: SA from scratch is worse than cold-only** because the SA
phases destroy tree structure and TBR can't fully recover in the
smoother IW landscape.

---

## Part 4: Multi-cycle PCSA with best-tree restart

The key innovation: repeat SA+TBR polish cycles, restarting each cycle
from the best tree found so far (like ratchet's perturbation pattern).

### Algorithm

```
tree ← Wagner + TBR convergence
best ← tree
for cycle in 1..N:
    tree ← copy(best)
    tree ← SA(tree, t_start → 0)   # stochastic perturbation
    tree ← TBR_polish(tree)         # reconverge
    if score(tree) < score(best):
        best ← tree
return best
```

### Results: 125 tips EW (project3253, 10 seeds)

| Approach | Mean | SD | Wins/10 |
|----------|------|----|---------|
| Cold-only (TBR×20) | 5412 | 119 | 1 |
| PCSA×3 | 5331 | 63 | 3 |
| **PCSA×5** | **5304** | **63** | **6** |
| Cold+SA×3 (converge then SA) | 5322 | 46 | 1 |

### Results: 205 tips EW (project3763, 10 seeds)

| Approach | Mean | SD | Wins/10 |
|----------|------|----|---------|
| Cold-only (TBR×20) | 1413 | 62 | 1 |
| PCSA×3 | 1373 | 12 | 4 |
| **PCSA×5** | **1367** | **10** | **4** |
| Cold+SA×3 | 1376 | 14 | 2 |

### Equal-time control

Cold-only with 40 rounds (double TBR budget) barely improves over 20
rounds — the chain is stuck at its local optimum. Extra TBR rounds
cannot escape a converged basin.

| 205 tips EW (5 seeds) | Mean |
|-----------------------|------|
| Cold-only (TBR×20) | 1456 |
| Cold-only (TBR×40) | 1451 |
| PCSA×3 | **1399** |

### Key insights

1. **PCSA is a smooth ratchet.** Same principle as ratchet (perturb +
   reconverge + keep best), different perturbation mechanism (SA-style
   stochastic acceptance vs character reweighting).

2. **Variance reduction is the strongest signal.** SD drops from 62 to 10
   at 205 tips. Multiple cycles systematically explore different basins.

3. **More cycles = better**, with diminishing returns. 5 cycles is
   meaningfully better than 3.

4. **Cold+SA combo is competitive but PCSA alone is slightly better.**
   The pre-convergence doesn't hurt but the SA cycles dominate.

5. **Benefit is strongest under EW at 125+ tips.** Under IW, TBR already
   navigates the smoother landscape effectively.

---

## Recommendations

### PR #215 (feature/parallel-temper)

Keep disabled by default (`ptChains = 0`). The Boltzmann swap mechanism
has no path to effectiveness under EW parsimony. Under IW with high
T_cold it marginally works but doesn't beat cold-only.

### Pipeline integration of PCSA

Add SA as a perturbation phase in `run_single_replicate()`, after drift
and before final TBR polish. Proposed design:

- **Trigger:** EW mode, ≥100 tips (where TBR plateaus are common)
- **Parameters:** `t_start` = 10–20, `t_end` = 0, `n_phases` = 5,
  `moves_per_phase` = n_tip
- **Cycles:** 3–5 per replicate (with best-tree restart)
- **Budget:** SA cycles replace some ratchet cycles rather than adding
  to the total budget — the two mechanisms are complementary perturbation
  strategies

### Future work

- Test whether alternating SA and ratchet perturbation (character
  reweighting) within the same replicate discovers trees that neither
  mechanism finds alone
- Evaluate at 300+ tips where both EW and IW may plateau
- Tune SA temperature schedule adaptively based on acceptance rate
