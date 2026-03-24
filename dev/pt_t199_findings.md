# T-199 Findings: PT Diagnostic Profiling

## Executive Summary

**Standard parallel tempering with Boltzmann acceptance is fundamentally
incompatible with discrete parsimony's greedy cold chain.** Pair 0–1
(cold↔warm) swap acceptance is exactly 0% across all 8 test matrices
(40–385 tips), all temperature ladders tested, and all moves-per-round
settings — including 1 move per round at T=0.5.

A single cold-chain-only run matches the 4-chain PT result at 80% of the
wall time. Hot chains contribute zero improvements.

## Detailed Results

### Phase 1: Default ladder {0, 3, 9, 27} across 8 matrices

| Matrix | Tips | Chars | Best | Pair 0-1 | Overall swaps | Cold improve | Cold ms | Hot ms | Total ms |
|--------|------|-------|------|----------|---------------|-------------|---------|--------|----------|
| project3419 | 40 | 368 | 2062 | 0% | 2/30 | 0 | 78 | 16 | 94 |
| project2604 | 43 | 307 | 1115 | 0% | 12/30 | 0 | 320 | 74 | 395 |
| project4531 | 71 | 256 | 885 | 0% | 4/30 | 0 | 240 | 37 | 277 |
| project3741 | 86 | 110 | 849 | 0% | 6/30 | 0 | 491 | 268 | 759 |
| project3253 | 125 | 394 | 5435 | 0% | 10/30 | 0 | 8013 | 1219 | 9233 |
| project1221 | 150 | 252 | 874 | 0% | 6/30 | 0 | 4160 | 1131 | 5292 |
| project3763 | 205 | 105 | 1550 | 0% | 7/30 | 0 | 8910 | 1186 | 10096 |
| project2722 | 385 | 520 | 5584 | 0% | 5/30 | 0 | 220179 | 9601 | 229782 |

Note: swaps succeed between hot chain pairs (1-2, 2-3) at 10–90% rates.
The bottleneck is exclusively at the cold-warm boundary.

### Phase 2: Tighter temperature ladders

Tested {0, 0.5, 1, 2} and {0, 1, 3, 9} on 86, 125, 205 tips.
**Pair 0-1 acceptance: 0% in all cases.**

### Phase 3: Fewer moves per round

With tight ladder {0, 0.5, 1, 2}, tested moves_per_round = 5, 2, 1.
**Pair 0-1 acceptance: 0% even at 1 move per round.**

### Phase 4: Warm cold chain

Tested T_cold = 0.1 and 0.3 (with T_warm = 0.5, 1, 2).
- T_cold = 0.3: Metropolis prob ≈ 10^-256 (still effectively 0)
- T_cold = 0.1: Metropolis prob = 0.0

### Phase 5: Cold-only vs 4-chain PT (205 tips)

- 4-chain PT (20 rounds): best = 1550, time = 13659 ms, cold_improve = 0
- 1-chain cold only (20 rounds): best = 1550, time = 10914 ms
- **Hot chains are pure overhead.**

## Root Cause Analysis

In standard MCMC-based PT (e.g., Bayesian phylogenetics with continuous
branch lengths), adjacent temperatures have overlapping energy distributions
because small parameter perturbations cause small energy changes. This
enables successful swaps.

In discrete parsimony with TBR:
1. **Cold chain (T=0) converges almost instantly** to a local optimum via
   deterministic TBR (typically in round 0–1).
2. **Any T > 0 chain rapidly diverges** because each accepted stochastic TBR
   move changes many edges simultaneously, incurring a score penalty of
   tens to hundreds of steps.
3. **Score gap at pair 0-1 boundary**: 89–1400+ steps after just 1 round.
4. **Metropolis acceptance probability**: exp(-gap × Δβ) ≈ 0 (underflows to
   zero or ≈ 10^-256).

The fundamental incompatibility is that parsimony's discrete landscape with
TBR moves does not have the "gradually overlapping energy distributions"
property that PT requires.

## Implications for the PT Feature

The current PT implementation (T-190–T-193, PR #215) provides no benefit
and should either:

1. **Be removed / kept disabled by default** — the feature is harmless
   when `ptChains = 0` (default), but the code and parameter surface add
   complexity for zero benefit.

2. **Be redesigned with a parsimony-compatible exchange mechanism**:
   - **Score-transfer**: If any hot chain finds score ≤ cold_best, inject
     that tree into the cold chain (no Boltzmann criterion).
   - **Tree fusion**: Extract good subtree features from hot chains and
     graft them onto the cold chain (sector search across chains).
   - **Population diversity**: Hot chains as topology generators for a shared
     pool (essentially multi-start with stochastic restarts).

3. **Pivot to simulated annealing** — the SA infrastructure in ts_temper.h
   already exists and doesn't suffer from the swap problem.

## Recommendation

Close PR #215 or merge it with PT disabled by default, then redirect the
evaluation effort (T-200, T-201) toward:
- SA with adaptive cooling schedule (already implemented, needs tuning)
- Multi-start with topology transfer (closer to what TNT's "xmult" does)
- Cross-chain sector search (use hot chains to generate subtree alternatives
  for the cold chain without Metropolis swaps)
